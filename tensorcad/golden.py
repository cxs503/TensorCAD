"""Golden geometry regression suite: lock what the geometry pipeline produces.

A change that looks local (a sampler, a cap, an orientation fix) moves numbers nobody
re-checks, because analytic tests only probe the shapes they were written for. This
module is the counterweight: a fixed *golden library* of four hulls, an independent
reference for every asserted quantity, and a committed baseline of what the pipeline
measured the last time the baseline was blessed.

===================  ===========================================  =====================
case                 geometry path                                reference
===================  ===========================================  =====================
``wigley``           analytic surface -> loft -> cap -> integrals  closed forms
``series60``         Series-60 family -> loft -> cap -> integrals  closed forms
``sphere_ellipsoid`` SDF primitives -> marching tetrahedra         closed forms
``dtmb5415``         shipped reduced model of the 5415 dataset     dataset arrays (npz)
===================  ===========================================  =====================

Three disciplines are structural here, not conventions:

1. **No hand-copied criteria.** A reference is never a literal: every reference metric
   is produced by a callable named in ``reference_source``, and the checker re-invokes
   it on every run, refusing to pass when the stored value no longer reproduces
   (``reference_not_reproducible``). Editing a number in the baseline cannot buy a pass.
2. **Tolerances are derived, not chosen.** Acceptance tolerance = twice the observed
   coarse-to-fine change of the metric (a two-grid bound), with a 1e-4 relative floor
   for metrics that do not move. The derivation is recomputed and printed every run.
3. **Fail closed.** A missing baseline, a case the baseline does not know, a case the
   baseline knows and the code forgot, a missing npz, a baseline entry that does not lock
   the mesh (no ``geometry_sha256``) — all are errors, never "skipped".

Usage::

    python -m tensorcad.golden --check            # CI / git hook gate
    python -m tensorcad.golden --check --quick    # coarse level only (fast)
    python -m tensorcad.golden --bless            # record a baseline (refuses on errors)
    python -m tensorcad.golden --list
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import platform
import subprocess
import sys
import time
from collections.abc import Iterable, Mapping
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, ClassVar

SCHEMA = "tensorcad.golden/1"
ROOT = Path(__file__).resolve().parents[1]
BASELINE = ROOT / "docs" / "golden" / "tensorcad_geometry_golden.json"

#: Relative floor on an acceptance tolerance: keeps ``2 * |fine - coarse|`` from
#: collapsing to ~0 on an exactly-integrated metric.
TOLERANCE_FLOOR = 1e-4
#: Relative tolerance of the drift comparison against the blessed measurements. Geometry
#: is deterministic; the only legitimate movement is floating point (libm/BLAS).
DRIFT_TOLERANCE = 1e-9
#: Relative tolerance for "the stored reference equals what its source computes now".
REFERENCE_TOLERANCE = 1e-9
#: Hard ceiling on a tolerance a case may *declare* instead of deriving it from the two-grid
#: change.  When a reference is a published table rather than a closed form, the measured
#: value has no discretisation signal to bound it: the band has to be quoted from the
#: contract that owns the reference (the dataset's own 3 % particulars standard, say).  The
#: ceiling is what keeps that escape hatch from becoming "declare 40 % and pass": it sits
#: well below the platform's 15 % reject threshold, where an error is no longer a tolerance
#: question, and above the few-percent band a four-significant-figure reference table and a
#: reconstructed hull can honestly support.
DECLARED_TOLERANCE_CAP = 0.05
#: Prefix of the ``basis`` a *declared* tolerance carries in the baseline.  It is the
#: marker the checker reads to know which derivation a recorded band came from.
DECLARED_BASIS_PREFIX = "declared contract"

#: Safety factor of the grid-convergence index (GCI).  Roache's three-grid recommendation
#: (Verification and Validation in Computational Fluid Dynamics, 1998/2002): with an
#: observed order fitted from three levels a factor of 1.25 covers the order uncertainty.
#: The band it produces is a *derived* criterion -- the extrapolated measurement must agree
#: with the analytic reference within the band that the measurements' own convergence
#: justifies -- so it is not a chosen number.
GCI_SAFETY_FACTOR = 1.25
ERROR, WARNING, NOTICE = "error", "warning", "notice"
LEVELS: tuple[str, ...] = ("coarse", "fine")


# --------------------------------------------------------------------------------------
# records
# --------------------------------------------------------------------------------------
@dataclass(frozen=True)
class Finding:
    severity: str
    code: str
    case: str
    detail: str
    metric: str | None = None

    def to_dict(self) -> dict[str, Any]:
        data: dict[str, Any] = {
            "severity": self.severity,
            "code": self.code,
            "case": self.case,
            "detail": self.detail,
        }
        if self.metric is not None:
            data["metric"] = self.metric
        return data


@dataclass(frozen=True)
class MetricCheck:
    """One asserted quantity: value, reference, relative error and derived tolerance."""

    name: str
    value: float
    kind: str  # "analytic" | "invariant" | "drift"
    reference: float | None = None
    relative_error: float | None = None
    tolerance: float | None = None
    basis: str = ""
    #: Where the tolerance came from: "two-grid" (derived from the coarse-to-fine change on
    #: this run), "declared" (quoted from the contract that owns the reference, see
    #: :data:`DECLARED_TOLERANCE_CAP`) or "invariant" (written by the case).
    origin: str = "two-grid"

    @property
    def passed(self) -> bool:
        if self.kind == "drift":
            return True
        if self.relative_error is None or self.tolerance is None:
            return False
        return abs(self.relative_error) <= self.tolerance

    def to_dict(self) -> dict[str, Any]:
        data: dict[str, Any] = {
            "name": self.name,
            "value": _number(self.value),
            "kind": self.kind,
            "origin": self.origin,
            "passed": self.passed,
        }
        for key, value in (
            ("reference", self.reference),
            ("relative_error", self.relative_error),
            ("tolerance", self.tolerance),
        ):
            if value is not None:
                data[key] = _number(value)
        if self.basis:
            data["basis"] = self.basis
        return data


@dataclass
class CaseResult:
    name: str
    title: str
    reference_source: str
    metrics: list[MetricCheck] = field(default_factory=list)
    findings: list[Finding] = field(default_factory=list)
    measurements: dict[str, dict[str, float]] = field(default_factory=dict)
    reference: dict[str, float] = field(default_factory=dict)
    invariants: dict[str, tuple[float, float, str]] = field(default_factory=dict)
    geometry_sha256: str = ""
    known_limitations: list[dict[str, Any]] = field(default_factory=list)
    seconds: float = 0.0

    @property
    def passed(self) -> bool:
        return not any(finding.severity == ERROR for finding in self.findings)

    def to_dict(self) -> dict[str, Any]:
        return {
            "name": self.name,
            "title": self.title,
            "reference_source": self.reference_source,
            "passed": self.passed,
            "seconds": round(self.seconds, 3),
            "geometry_sha256": self.geometry_sha256,
            "measurements": {
                level: {key: _number(value) for key, value in values.items()}
                for level, values in self.measurements.items()
            },
            "reference": {key: _number(value) for key, value in self.reference.items()},
            "metrics": [metric.to_dict() for metric in self.metrics],
            "known_limitations": self.known_limitations,
            "findings": [finding.to_dict() for finding in self.findings],
        }


@dataclass
class SuiteResult:
    cases: list[CaseResult] = field(default_factory=list)
    findings: list[Finding] = field(default_factory=list)
    baseline: str = ""
    quick: bool = False
    seconds: float = 0.0

    @property
    def all_findings(self) -> list[Finding]:
        return list(self.findings) + [f for case in self.cases for f in case.findings]

    @property
    def passed(self) -> bool:
        return not any(finding.severity == ERROR for finding in self.all_findings)

    @property
    def errors(self) -> list[Finding]:
        return [finding for finding in self.all_findings if finding.severity == ERROR]

    def to_dict(self) -> dict[str, Any]:
        return {
            "schema": SCHEMA,
            "passed": self.passed,
            "quick": self.quick,
            "seconds": round(self.seconds, 3),
            "baseline": self.baseline,
            "cases": [case.to_dict() for case in self.cases],
            "findings": [finding.to_dict() for finding in self.findings],
        }


# --------------------------------------------------------------------------------------
# the golden library
# --------------------------------------------------------------------------------------
class GoldenCase:
    """A geometry plus an independent reference for the quantities it is asserted on."""

    name: str = ""
    title: str = ""
    reference_source: str = ""

    def measure(self, level: str) -> dict[str, float]:
        raise NotImplementedError

    def reference(self) -> dict[str, float]:
        raise NotImplementedError

    def invariants(self, values: Mapping[str, float]) -> dict[str, tuple[float, float, str]]:
        """Identities over the measurements: name -> (residual, tolerance, basis)."""
        return {}

    def reference_tolerance(
        self, values: Mapping[str, float], reference: Mapping[str, float]
    ) -> dict[str, tuple[float, str]]:
        """Declared acceptance band for reference metrics that carry no two-grid signal.

        A reference metric is bounded by ``2 x |fine - coarse|`` when refining a grid is what
        moves it.  A metric that is *not* refreshed by a finer grid -- because the reference
        is a published table and the measurement is the same reconstruction at both levels --
        has no such signal, so the derived bound collapses to the 1e-4 floor and would reject
        the reference's own precision.  A case that owns such a metric declares the band it
        is held to instead, as ``metric -> (tolerance, basis)``, where ``basis`` must name
        the contract the band is quoted from.  The engine refuses an unjustified or too-loose
        declaration (:data:`DECLARED_TOLERANCE_CAP`); everything not declared here keeps the
        derived two-grid bound.
        """
        return {}

    def known_limitations(self, values: Mapping[str, float]) -> list[dict[str, Any]]:
        return []

    def geometry_sha256(self) -> str:
        return ""

    def geometry_inputs(self) -> tuple[Path, ...]:
        """Data files this case reads, beyond the library code.

        A case whose geometry is generated from code alone reads nothing, and the gate is
        triggered by ``tensorcad/`` itself.  A case that measures a *shipped dataset* has
        inputs a change to which moves the measured geometry without touching a line of
        ``tensorcad/``: those must be declared here, because they are also what the trigger
        set (:func:`geometry_files`, :func:`staged_requires_golden`) and the source
        fingerprint are derived from.  Listing them by hand anywhere else would let the
        gate go quiet on a dataset edit.
        """
        return ()

    def levels(self) -> tuple[str, ...]:
        """The measurement ladder; the last level is the acceptance level.

        Two levels (``coarse``/``fine``) is the base discipline: the acceptance band is
        twice the *last* coarse-to-fine change.  A case whose two-grid band is too loose to
        see a small coherent regression declares a third level here; that level does not
        widen the acceptance band (which stays a last-step two-grid bound), it is what the
        grid-convergence assertion is extrapolated over.
        """
        return LEVELS

    def refinement_factor(self, coarse_level: str, fine_level: str) -> float | None:
        """Cell-size ratio ``h(coarse) / h(fine)``, or ``None`` when not declared.

        This is the ratio the extrapolation is raised to, so it is computed from the
        ladder, never typed: an invented ratio would silently rescale every band.
        """
        return None


class WigleyGolden(GoldenCase):
    """Classic Wigley hull: every hydrostatic quantity has a closed form."""

    name = "wigley"
    title = "Wigley hull (L=10, B=1, T=0.625, n=m=2): lofted and capped surface"
    reference_source = "tensorcad.wigley.wigley_analytic_hydrostatics (closed forms)"
    METRICS: ClassVar[dict[str, str]] = {
        "displacement_volume_m3": "displacement_volume",
        "waterplane_area_m2": "waterplane_area",
        "wetted_surface_area_m2": "wetted_surface_area",
        "block_coefficient": "block_coefficient",
        "prismatic_coefficient": "prismatic_coefficient",
        "longitudinal_centre_buoyancy_m": "longitudinal_centre_buoyancy",
        "vertical_centre_buoyancy_m": "vertical_centre_buoyancy",
        "maximum_section_area_m2": "maximum_section_area",
    }

    #: Ladder of grids.  The third level exists because the two-grid band on the analytic
    #: metrics was 4.3e-3 while the extrapolated value reproduces the closed forms to ~6e-6
    #: (tools/golden/convergence_study.py): the measured error was discretisation, not model
    #: error, and only a third level can say so.
    GRIDS: ClassVar[dict[str, tuple[int, int]]] = {
        "coarse": (41, 13),
        "fine": (161, 41),
        "finer": (321, 81),
    }

    def _parameters(self, level: str) -> Any:
        from . import wigley

        nx, nz = self.GRIDS[level]
        return wigley.reference_parameters(nx=nx, nz=nz)

    def measure(self, level: str) -> dict[str, float]:
        from . import wigley

        hydro = wigley.generate_wigley(self._parameters(level)).hydrostatics
        assert hydro is not None
        return _hydrostatics_values(hydro, self.METRICS)

    def reference(self) -> dict[str, float]:
        from . import wigley

        analytic = wigley.wigley_analytic_hydrostatics(wigley.reference_parameters()).as_dict()
        return {key: float(analytic[key]) for key in self.METRICS}

    def geometry_sha256(self) -> str:
        """Fingerprint of the acceptance grid the measurements are taken on."""
        from . import wigley

        return _mesh_sha256(wigley.generate_wigley(self._parameters(self.levels()[-1])).mesh)

    def levels(self) -> tuple[str, ...]:
        return ("coarse", "fine", "finer")

    def refinement_factor(self, coarse_level: str, fine_level: str) -> float:
        coarse = self.GRIDS[coarse_level][1] - 1
        fine = self.GRIDS[fine_level][1] - 1
        return fine / coarse


class Series60Golden(GoldenCase):
    """Series-60 family (Cb = 0.60): the exponents are re-derived, not typed."""

    name = "series60"
    title = "Series 60 hull (L=10, B=1.5, T=0.6, Cb=0.60): fitted exponent family"
    reference_source = "tensorcad.series60.analytic_* over the exponents fitted from Cb"
    METRICS: ClassVar[dict[str, str]] = {
        "displacement_volume_m3": "displacement_volume",
        "waterplane_area_m2": "waterplane_area",
        "block_coefficient": "block_coefficient",
        "prismatic_coefficient": "prismatic_coefficient",
        "maximum_section_area_m2": "maximum_section_area",
        "longitudinal_centre_buoyancy_m": "longitudinal_centre_buoyancy",
    }

    #: Ladder of grids.  This case carries the one documented blind spot of the two-grid
    #: criterion (``maximum_section_area_m2`` at 2.72e-2, docs/golden/mutation_audit.md):
    #: volume / Cb / mid-section converge at order ~1.1 (the keel exponent singularity of
    #: the section law), so a third level plus the extrapolation is the only way to hold
    #: them below 1 % without paying for a fourth grid.
    GRIDS: ClassVar[dict[str, tuple[int, int]]] = {
        "coarse": (81, 25),
        "fine": (161, 49),
        "finer": (321, 97),
    }

    def _parameters(self, level: str) -> Any:
        from . import series60

        nx, nz = self.GRIDS[level]
        return series60.Series60Parameters(nx=nx, nz=nz)

    def _exponents(self) -> tuple[float, float]:
        from . import series60

        parameters = self._parameters("coarse")
        return series60._interpolate_exponents(parameters.cb)

    def measure(self, level: str) -> dict[str, float]:
        from . import series60

        hydro = series60.generate_series60(self._parameters(level)).hydrostatics
        assert hydro is not None
        return _hydrostatics_values(hydro, self.METRICS)

    def reference(self) -> dict[str, float]:
        from . import series60

        parameters = self._parameters("coarse")
        n, m = self._exponents()
        return {
            "displacement_volume_m3": float(series60.analytic_volume(  # type: ignore[no-untyped-call]
                parameters.length, parameters.beam, parameters.draft, n, m
            )),
            "waterplane_area_m2": float(series60.analytic_waterplane_area(  # type: ignore[no-untyped-call]
                parameters.length, parameters.beam, n
            )),
            "block_coefficient": float(series60.analytic_cb(n, m)),  # type: ignore[no-untyped-call]
            "prismatic_coefficient": float(series60.analytic_cp(n)),  # type: ignore[no-untyped-call]
            "maximum_section_area_m2": float(
                series60.analytic_maximum_section_area(  # type: ignore[no-untyped-call]
                    parameters.beam, parameters.draft, m
                )
            ),
            "longitudinal_centre_buoyancy_m": parameters.length / 2.0,
        }

    def invariants(self, values: Mapping[str, float]) -> dict[str, tuple[float, float, str]]:
        parameters = self._parameters("coarse")
        cb = values["block_coefficient"]
        midship = values["maximum_section_area_m2"] / (parameters.beam * parameters.draft)
        return {
            "identity_cb_equals_cp_times_cm": (
                cb - values["prismatic_coefficient"] * midship,
                1e-9,
                "Cb = Cp·Cm with Cm = A_midship / (B·T)",
            ),
            "identity_volume_equals_cb_lbt": (
                values["displacement_volume_m3"] - cb * parameters.length * parameters.beam
                * parameters.draft,
                1e-9,
                "V = ∫dV must equal Cb·L·B·T on the same mesh",
            ),
        }

    def geometry_sha256(self) -> str:
        """Fingerprint of the acceptance grid the measurements are taken on."""
        from . import series60

        return _mesh_sha256(series60.generate_series60(self._parameters(self.levels()[-1])).mesh)

    def levels(self) -> tuple[str, ...]:
        return ("coarse", "fine", "finer")

    def refinement_factor(self, coarse_level: str, fine_level: str) -> float:
        coarse = self.GRIDS[coarse_level][1] - 1
        fine = self.GRIDS[fine_level][1] - 1
        return fine / coarse


class SphereEllipsoidGolden(GoldenCase):
    """SDF primitives through marching tetrahedra: the closed forms are exact."""

    name = "sphere_ellipsoid"
    title = "Sphere (r=1) and ellipsoid (2.0, 1.5, 0.5): SDF -> marching tetrahedra"
    reference_source = (
        "(4/3)πr³ and (4/3)πabc over the exact primitive fields; the areas from the exact "
        "semi-axes: 4πr² for the sphere and "
        "tensorcad.analytic_solids.ellipsoid_surface_area for the triaxial ellipsoid "
        "(converged Gauss-Legendre quadrature of the exact area element, cross-checked "
        "against the spheroid closed form)"
    )
    AXES = (2.0, 1.5, 0.5)
    RADIUS = 1.0

    def _resolution(self, level: str) -> int:
        return {"coarse": 24, "fine": 48}[level]

    def measure(self, level: str) -> dict[str, float]:
        from .geometry_field import EllipsoidField, SphereField, to_mesh
        from .hydrostatics import compute_hydrostatics

        resolution = self._resolution(level)
        values: dict[str, float] = {}
        for name, geometry in (
            ("sphere", SphereField(radius=self.RADIUS, center=(0.0, 0.0, 0.0))),
            ("ellipsoid", EllipsoidField(semi_axes=self.AXES, center=(0.0, 0.0, 0.0))),
        ):
            mesh = to_mesh(geometry, resolution=resolution).mesh
            hydro = compute_hydrostatics(mesh, waterline_z=-1e3, closed=True)
            values[f"{name}_volume_m3"] = abs(hydro.displacement_volume)
            values[f"{name}_surface_area_m2"] = _closed_surface_area(mesh)
            for axis, extent in zip("xyz", _extents(mesh), strict=True):
                values[f"{name}_extent_{axis}_m"] = extent
        return values

    def reference(self) -> dict[str, float]:
        from .analytic_solids import ellipsoid_surface_area, ellipsoid_volume

        reference = {
            "sphere_volume_m3": 4.0 / 3.0 * math.pi * self.RADIUS**3,
            "sphere_surface_area_m2": 4.0 * math.pi * self.RADIUS**2,
            "ellipsoid_volume_m3": ellipsoid_volume(self.AXES),
            "ellipsoid_surface_area_m2": ellipsoid_surface_area(self.AXES),
        }
        for axis, extent in zip("xyz", (2 * self.RADIUS,) * 3, strict=True):
            reference[f"sphere_extent_{axis}_m"] = extent
        for axis, extent in zip("xyz", (2.0 * axis for axis in self.AXES), strict=True):
            reference[f"ellipsoid_extent_{axis}_m"] = extent
        return reference


    def geometry_sha256(self) -> str:
        """Fingerprint of both sampled fields (sphere + triaxial ellipsoid, fine grid).

        Two meshes, one case: each mesh is hashed and the pair is folded into a single
        digest so a change to either body drifts the case.
        """
        from .geometry_field import EllipsoidField, SphereField, to_mesh

        digest = hashlib.sha256()
        for name, geometry in (
            ("sphere", SphereField(radius=self.RADIUS, center=(0.0, 0.0, 0.0))),
            ("ellipsoid", EllipsoidField(semi_axes=self.AXES, center=(0.0, 0.0, 0.0))),
        ):
            mesh = to_mesh(geometry, resolution=self._resolution("fine")).mesh
            digest.update(name.encode())
            digest.update(b"\0")
            digest.update(_mesh_sha256(mesh).encode())
        return digest.hexdigest()


class DTMB5415Golden(GoldenCase):
    """The shipped reduced model of the DTMB 5415 design dataset.

    The reference is the dataset itself (``docs/data/dtmb5415/holdout_fixture.npz``): the
    holdout residual the dataset records for one configuration. Absolute displacement has
    no published counterpart inside the repository, so it is drift-locked only, and the
    mesh is asserted to be what it is — an open hull sheet — so that its tetrahedron sum
    is never mistaken for a displacement.
    """

    name = "dtmb5415"
    title = "DTMB 5415 baseline hull: shipped rank-40 reduced model (npz, no 1 GB source)"
    #: Metric -> (measured key, published key, what the comparison is worth).
    #: Each published key is derived by :func:`tensorcad.dtmb5415.readme_reference` from the
    #: dataset README's model-scale table, never typed here.
    PARTICULARS: ClassVar[dict[str, tuple[str, str]]] = {
        "length_stations_m": ("Lpp_m", "station span vs published Lpp: the mesh spans the "
                                       "station grid, which runs ~0.9 % past Lpp"),
        "waterline_beam_m": ("beam_m", "beam at the design waterline vs published moulded "
                                       "B_max: the two differ by the flare above the WL"),
        "displacement_volume_m3": ("displacement_volume_m3", "displacement volume at the "
                                                            "design waterline"),
        "block_coefficient": ("cb", "Cb over the published L.B.T vs Cb over the measured "
                                    "station span/WL beam/moulded draft"),
    }
    reference_source = (
        "docs/data/dtmb5415/holdout_fixture.npz:per_config_rms_m[0]; principal particulars "
        "from tensorcad.dtmb5415.readme_reference() (dataset README model-scale table x lambda)"
    )
    CONFIG = 0
    MIRROR_AXIS = 1

    #: The shipped reduced model: ``design_space.npz`` is the rank-40 basis the hull is
    #: rebuilt from, ``design_space_metrics.json`` carries the fit that produced it, and
    #: ``holdout_fixture.npz`` supplies the reference residual of :meth:`measure`.  All
    #: three are inputs of this case: editing any of them changes what the gate asserts.
    DATA_DIR = ROOT / "docs" / "data" / "dtmb5415"

    def geometry_inputs(self) -> tuple[Path, ...]:
        return (
            self.DATA_DIR / "design_space.npz",
            self.DATA_DIR / "design_space_metrics.json",
            self.DATA_DIR / "holdout_fixture.npz",
        )

    def _fixture(self) -> Any:
        import numpy as np

        return np.load(ROOT / "docs" / "data" / "dtmb5415" / "holdout_fixture.npz")

    def _space(self) -> Any:
        from .dtmb5415 import DTMB5415Space

        return DTMB5415Space()

    def measure(self, level: str) -> dict[str, float]:
        import numpy as np

        from . import mesh as meshlib

        space = self._space()
        fixture = self._fixture()
        cloud = np.asarray(space.base_cloud, dtype=float)
        mesh = space.mesh()
        vertices = np.asarray(mesh.vertices, dtype=float)
        hydro = space.hydrostatics()
        values = {
            "displacement_volume_m3": float(hydro["displacement_volume_m3"]),
            "block_coefficient": float(hydro["cb"]),
            "waterplane_area_m2": float(hydro["waterplane_area_m2"]),
            "longitudinal_centre_buoyancy_m": float(hydro["lcb_m"]),
            "waterline_beam_m": float(hydro["beam_waterline_m"]),
            "draft_m": float(hydro["draft_m"]),
            "length_stations_m": float(hydro["length_stations_m"]),
            "mesh_vertex_count": float(len(mesh.vertices)),
            "mesh_face_count": float(len(mesh.faces)),
            # The shell is an open sheet (the hull *surface*, not a solid). Its
            # tetrahedron sum is a signed integral over an open surface: never a volume.
            # It is stored under a name that cannot be mistaken for one.
            "signed_tetra_sum_open_shell_m3": float(meshlib.total_signed_volume(mesh)),
            "mesh_is_watertight": float(meshlib.is_watertight(mesh)),
            "boundary_edge_count": float(len(meshlib.boundary_edges(mesh))),
            "boundary_vertex_count": float(
                len({index for edge in meshlib.boundary_edges(mesh) for index in edge})
            ),
            "boundary_component_count": float(_boundary_component_count(mesh)),
        }
        for axis, index in (("x", 0), ("y", 1), ("z", 2)):
            values[f"cloud_extent_{axis}_m"] = float(np.ptp(cloud[:, index]))
            values[f"mesh_extent_{axis}_m"] = float(np.ptp(vertices[:, index]))
        design = np.asarray(fixture["design_vector"][self.CONFIG], dtype=float)
        truth = np.asarray(fixture["true_cloud"][self.CONFIG], dtype=float)
        predicted = np.asarray(space.cloud(design), dtype=float)
        values["holdout_rms_m"] = float(
            np.sqrt(np.mean(np.sum((predicted - truth) ** 2, axis=1)))
        )
        return values

    def reference(self) -> dict[str, float]:
        from .dtmb5415 import readme_reference

        fixture = self._fixture()
        published = readme_reference()
        reference = {"holdout_rms_m": float(fixture["per_config_rms_m"][self.CONFIG])}
        for measured, (published_key, _) in self.PARTICULARS.items():
            reference[measured] = float(published[published_key])
        return reference

    def reference_tolerance(
        self, values: Mapping[str, float], reference: Mapping[str, float]
    ) -> dict[str, tuple[float, str]]:
        """The dataset's own particulars band, quoted rather than re-derived.

        ``draft_m`` is deliberately absent: it is an *input* of the hydrostatics (the
        waterline the hull is cut at), so asserting it against the README draft would be a
        constant compared with itself -- no independent information.  It stays drift-locked
        like every other measurement.
        """
        from .dtmb5415 import README_PARTICULARS_TOLERANCE, DRAFT_M

        band = README_PARTICULARS_TOLERANCE
        declared: dict[str, tuple[float, str]] = {}
        for measured, (_, caveat) in self.PARTICULARS.items():
            declared[measured] = (
                band,
                f"{band:.0%} band of the hull's own pre-declared geometry contract "
                f"(tools/dtmb5415_param/verify_geometry.py: 3 % accept, 15 % reject), "
                f"reference = dataset README table x lambda; {caveat}; draft_m "
                f"(= {DRAFT_M:.4f} m) is an input and is not asserted",
            )
        return declared

    def invariants(self, values: Mapping[str, float]) -> dict[str, tuple[float, float, str]]:
        axis = "xyz"[self.MIRROR_AXIS]
        invariants: dict[str, tuple[float, float, str]] = {
            "identity_cb_equals_volume_over_lbt": (
                values["block_coefficient"]
                - values["displacement_volume_m3"]
                / (values["length_stations_m"] * values["waterline_beam_m"] * values["draft_m"]),
                1e-9,
                "Cb = ∇ / (L·B·T) on the reported particulars",
            ),
            f"identity_mirror_doubles_{axis}_extent": (
                values[f"mesh_extent_{axis}_m"] - 2.0 * values[f"cloud_extent_{axis}_m"],
                1e-9,
                f"the full hull is the mirror of the half cloud about {axis}",
            ),
        }
        for other in (name for name in "xyz" if name != axis):
            invariants[f"identity_mirror_preserves_{other}_extent"] = (
                values[f"mesh_extent_{other}_m"] - values[f"cloud_extent_{other}_m"],
                1e-9,
                f"mirroring about {axis} leaves the {other} extent unchanged",
            )
        # Topology, not numerics: the shipped shell is the hull *surface*. If a future
        # generator closes it, this fires and the displacement can be asserted instead.
        invariants["identity_open_shell_is_not_a_solid"] = (
            values["mesh_is_watertight"],
            0.0,
            "the shipped 5415 mesh is an open hull sheet, so no closed-surface volume exists",
        )
        invariants["identity_open_shell_has_one_boundary_component"] = (
            values["boundary_component_count"] - 1.0,
            0.0,
            "the open sheet carries a single boundary loop (deck edge + stem seam)",
        )
        return invariants

    def geometry_sha256(self) -> str:
        """Fingerprint of the shipped shell: vertices, faces and their order."""
        return _mesh_sha256(self._space().mesh())

    def known_limitations(self, values: Mapping[str, float]) -> list[dict[str, Any]]:
        return [
            {
                "id": "displacement-withheld-open-shell",
                "detail": (
                    "the shipped 5415 shell is an open sheet, not a solid: one boundary loop "
                    f"of {int(values['boundary_edge_count'])} edges (deck edge plus the stem "
                    "seam). Surface-integral displacement is therefore withheld: the "
                    "displacement asserted on this case is the dataset/station-area one, "
                    "cross-checked against the README particulars table (see reference), not "
                    "an integral over this shell. Passing closed=True to "
                    "compute_hydrostatics on this shell returns the signed tetrahedron sum of "
                    "an open surface (negative here), which is not a volume; the previous "
                    "'relative_discrepancy' of 2.6 was measuring that artifact, not an "
                    "integration error. Closing the sheet is not a one-liner: after welding, "
                    "mesh_repair.boundary_loops refuses this boundary (a vertex is walked by "
                    "two walks at the stem pinch), so a watertight repair needs the seam "
                    "stitched first."
                ),
                "displacement_volume_m3": _number(values["displacement_volume_m3"]),
                "signed_tetra_sum_open_shell_m3": _number(values["signed_tetra_sum_open_shell_m3"]),
                "boundary_edge_count": _number(values["boundary_edge_count"]),
                "boundary_component_count": _number(values["boundary_component_count"]),
                "asserted": False,
            }
        ]


def golden_cases() -> list[GoldenCase]:
    return [WigleyGolden(), Series60Golden(), SphereEllipsoidGolden(), DTMB5415Golden()]


# --------------------------------------------------------------------------------------
# helpers
# --------------------------------------------------------------------------------------
def _boundary_component_count(mesh: Any) -> int:
    """Connected components of the mesh boundary edges (a topological fact).

    ``mesh_repair.boundary_loops`` walks loops and rejects a boundary that pinches at a
    vertex, so this counts components instead: one component means one hole, however the
    walk happens to be oriented.
    """
    from . import mesh as meshlib

    parent: dict[int, int] = {}

    def find(item: int) -> int:
        parent.setdefault(item, item)
        while parent[item] != item:
            parent[item] = parent[parent[item]]
            item = parent[item]
        return item

    for a, b in meshlib.boundary_edges(mesh):
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[ra] = rb
    return len({find(node) for node in parent})


def _hydrostatics_values(hydro: Any, metrics: Mapping[str, str]) -> dict[str, float]:
    """Map reported metric names to the :class:`Hydrostatics` attributes behind them."""
    return {name: float(getattr(hydro, attribute)) for name, attribute in metrics.items()}


def _number(value: float) -> float | int:
    if isinstance(value, bool):  # pragma: no cover - defensive
        return int(value)
    if math.isfinite(value) and float(value).is_integer() and abs(value) < 1e15:
        return int(value)
    return float(value)


def _relative(delta: float, reference: float) -> float:
    return abs(delta) if reference == 0.0 else abs(delta) / abs(reference)


def _closed_surface_area(mesh: Any) -> float:
    vertices = mesh.vertices
    total = 0.0
    for a, b, c in mesh.faces:
        pa, pb, pc = vertices[a], vertices[b], vertices[c]
        ab = (pb[0] - pa[0], pb[1] - pa[1], pb[2] - pa[2])
        ac = (pc[0] - pa[0], pc[1] - pa[1], pc[2] - pa[2])
        cx = ab[1] * ac[2] - ab[2] * ac[1]
        cy = ab[2] * ac[0] - ab[0] * ac[2]
        cz = ab[0] * ac[1] - ab[1] * ac[0]
        total += 0.5 * math.sqrt(cx * cx + cy * cy + cz * cz)
    return total


def _extents(mesh: Any) -> tuple[float, float, float]:
    xs = [vertex[0] for vertex in mesh.vertices]
    ys = [vertex[1] for vertex in mesh.vertices]
    zs = [vertex[2] for vertex in mesh.vertices]
    return (max(xs) - min(xs), max(ys) - min(ys), max(zs) - min(zs))


def _mesh_sha256(mesh: Any) -> str:
    """Quantised mesh hash: equal geometry hashes equal even across machines."""
    digest = hashlib.sha256()
    for vertex in mesh.vertices:
        digest.update(",".join(f"{coordinate:.9f}" for coordinate in vertex).encode())
        digest.update(b";")
    for face in mesh.faces:
        digest.update(",".join(str(index) for index in face).encode())
        digest.update(b";")
    return digest.hexdigest()


def _declared_inputs(root: Path = ROOT) -> list[Path]:
    """The data inputs every case declares, restricted to files inside ``root``."""
    inputs: list[Path] = []
    for case in golden_cases():
        for path in case.geometry_inputs():
            path = Path(path)
            if path in inputs:
                continue
            if root != ROOT:
                try:
                    path = root / path.relative_to(ROOT)
                except ValueError:
                    continue
            inputs.append(path)
    return sorted(inputs)


def geometry_files(root: Path = ROOT) -> list[Path]:
    """Every file whose change obliges a golden run.

    The library code, the suite's own tests, and the data inputs the cases declare
    (:meth:`GoldenCase.geometry_inputs`).  The last group is what makes the fingerprint a
    statement about the *inputs* of the measurement, not only about the code that reads
    them: a dataset edit moves neither file name nor module, so a code-only fingerprint
    would call it unchanged.
    """
    files = sorted((root / "tensorcad").rglob("*.py"))
    files.append(root / "tests" / "test_golden_regression.py")
    files.extend(_declared_inputs(root))
    return [path for path in files if path.is_file()]


def geometry_fingerprint(root: Path = ROOT) -> str:
    digest = hashlib.sha256()
    for path in geometry_files(root):
        digest.update(path.relative_to(root).as_posix().encode())
        digest.update(b"\0")
        digest.update(path.read_bytes())
        digest.update(b"\0")
    return digest.hexdigest()


def staged_requires_golden(paths: Iterable[str], root: Path = ROOT) -> bool:
    """Does this set of staged paths touch geometry, the suite, or the baseline?"""
    for raw in paths:
        path = raw.strip()
        if not path:
            continue
        # anything under tensorcad/ (including a module that does not exist yet), the
        # suite's own tests, and the baseline itself
        if path.startswith("tensorcad/") and path.endswith(".py"):
            return True
        if path in ("tests/test_golden_regression.py", BASELINE.relative_to(ROOT).as_posix()):
            return True
        # a declared data input of any case (see GoldenCase.geometry_inputs): staging one of
        # these edits the shipped dataset the gate measures, which is a geometry change
        # that touches no file under tensorcad/ and so would otherwise pass unnoticed
        for declared in _declared_inputs(root):
            try:
                if declared.relative_to(root).as_posix() == path:
                    return True
            except ValueError:
                continue
        if path.startswith("docs/golden/"):
            return True
    return False


def git_commit(root: Path = ROOT) -> str:
    completed = subprocess.run(
        ["git", "rev-parse", "HEAD"], cwd=root, capture_output=True, text=True, check=False
    )
    return completed.stdout.strip() if completed.returncode == 0 and completed.stdout.strip() else "unknown"


def _staged_paths(explicit: str | None) -> list[str]:
    if explicit is not None:
        return Path(explicit).read_text(encoding="utf-8").splitlines()
    completed = subprocess.run(
        ["git", "diff", "--cached", "--name-only", "--diff-filter=ACMR"],
        cwd=ROOT,
        capture_output=True,
        text=True,
        check=False,
    )
    return completed.stdout.splitlines()


# --------------------------------------------------------------------------------------
# evaluation
# --------------------------------------------------------------------------------------
def _tolerance(
    reference: float, coarse: float | None, fine: float
) -> tuple[float, str]:
    """Two-grid acceptance tolerance: twice the observed coarse-to-fine change.

    A discretisation error that halves between two grids is over-estimated by this bound;
    the factor two is the safety factor of a two-grid estimate. A metric that does not
    move between the grids falls back to the relative floor.
    """
    if coarse is None:
        return TOLERANCE_FLOOR, (
            f"reproduction floor {TOLERANCE_FLOOR:g}: no coarse level in this run "
            "(quick mode or no coarse/fine axis)"
        )
    change = _relative(fine - coarse, reference)
    if change > TOLERANCE_FLOOR / 2.0:
        return 2.0 * change, (
            f"two-grid: 2 x |fine - coarse| / |reference| = 2 x {change:.3e} = {2.0 * change:.3e}"
        )
    return TOLERANCE_FLOOR, (
        f"two-grid change {change:.3e} is below the floor, so the floor "
        f"{TOLERANCE_FLOOR:g} applies"
    )


def _observed_order(change1: float, change2: float, r1: float, r2: float) -> float | None:
    """Solve ``(r1^p - 1) / (1 - r2^-p) = |change1| / |change2|`` for the order ``p``.

    The two refinement ratios are not assumed equal (the ladders are declared per case), so
    the order is solved rather than read off a ratio of logs.  ``None`` means the two
    consecutive changes cannot be explained by any positive order -- a metric that does not
    move, or moves the wrong way -- and no extrapolation is claimed for it.
    """
    if abs(change1) <= 0.0 or abs(change2) <= 0.0:
        return None
    target = abs(change1) / abs(change2)

    def residual(p: float) -> float:
        return (r1**p - 1.0) / (1.0 - r2 ** (-p)) - target

    low, high = 1e-3, 12.0
    if not (residual(low) < 0.0 < residual(high)):
        return None
    for _ in range(200):
        mid = 0.5 * (low + high)
        if residual(mid) < 0.0:
            low = mid
        else:
            high = mid
    return 0.5 * (low + high)


def _gci_estimates(
    case: "GoldenCase",
    levels: tuple[str, ...],
    measurements: Mapping[str, Mapping[str, float]],
    reference: Mapping[str, float],
    findings: list["Finding"],
) -> dict[str, tuple[float, float, str]]:
    """Richardson extrapolant and its grid-convergence band, per reference metric.

    With three levels ``v1, v2, v3`` and refinement ratios ``r1, r2`` the observed order
    ``p`` makes two extrapolants: ``X12`` from the coarse pair and ``X23`` from the fine
    pair.  The asserted value is ``X23`` (the one built on the finest information) and the
    band is Roache's three-grid index ``Fs * |X23 - v3| / |reference|``: the extrapolated
    value must agree with the analytic reference at least as well as the finest measurement
    is estimated to be converged.  Both quantities are functions of the measurements, so
    the criterion cannot be hand-picked -- and it is checked against its own derivation on
    every run (:func:`_check_recorded_extrapolation`).
    """
    ladder = levels[-3:]
    r1 = case.refinement_factor(ladder[0], ladder[1])
    r2 = case.refinement_factor(ladder[1], ladder[2])
    if r1 is None or r2 is None or r1 <= 1.0 or r2 <= 1.0:
        findings.append(
            Finding(
                ERROR,
                "richardson_ratio_missing",
                case.name,
                f"the ladder {ladder} has no declared refinement factor, so no "
                "grid-convergence estimate can be taken over it",
            )
        )
        return {}
    estimates: dict[str, tuple[float, float, str]] = {}
    skipped: list[str] = []
    for key in sorted(reference):
        try:
            v1 = float(measurements[ladder[0]][key])
            v2 = float(measurements[ladder[1]][key])
            v3 = float(measurements[ladder[2]][key])
        except KeyError:
            skipped.append(key)
            continue
        p = _observed_order(v2 - v1, v3 - v2, r1, r2)
        if p is None:
            skipped.append(key)
            continue
        x12 = v2 - (v2 - v1) / (r1**p - 1.0)
        x23 = v3 - (v3 - v2) / (r2**p - 1.0)
        scale = abs(reference[key]) if reference[key] else 1.0
        # Two independent error terms, both functions of the measurements: the estimated
        # remaining discretisation error of the finest level (|X23 - v3|) and the observed
        # instability of the extrapolant between the two refinement steps (|X23 - X12|).
        # The second is what stops a non-asymptotic sequence from claiming convergence it
        # does not have: when the order is drifting, the two extrapolants disagree and the
        # band widens until it is honest again.
        remaining = abs(x23 - v3)
        instability = max(remaining, abs(x23 - x12))
        band = max(GCI_SAFETY_FACTOR * instability / scale, TOLERANCE_FLOOR)
        basis = (
            f"GCI: {GCI_SAFETY_FACTOR:g} x max(|extrapolated - {ladder[2]}|, "
            f"|X23 - X12|) / |reference| = {GCI_SAFETY_FACTOR:g} x "
            f"{instability / scale:.3e} (remaining {remaining / scale:.3e}; observed order "
            f"p = {p:.2f} over r = {r1:.4f}, {r2:.4f}: {x12:.9g} then {x23:.9g})"
        )
        estimates[key] = (x23, band, basis)
    if skipped:
        findings.append(
            Finding(
                WARNING,
                "richardson_unavailable",
                case.name,
                f"no grid-convergence estimate for {sorted(skipped)}: those metrics do not "
                "move consistently between the three levels (exactly integrated, or moving "
                "the wrong way), so their two-grid band is the only criterion in force",
            )
        )
    return estimates


def _check_recorded_extrapolation(
    case: "GoldenCase",
    blessed: Mapping[str, Any],
    reference: Mapping[str, float],
    findings: list["Finding"],
) -> None:
    """The recorded grid-convergence bands must be the ones this code derives.

    Same discipline as :func:`_check_recorded_criterion`: the baseline stores the bands a
    ``--bless`` derived, and every run -- quick included, since the derivation only needs
    the blessed measurements -- recomputes them.  A hand-widened GCI band is refused.
    """
    levels = case.levels()
    recorded = blessed.get("extrapolated") or {}
    if len(levels) < 3:
        if recorded:
            findings.append(
                Finding(
                    ERROR,
                    "extrapolation_unknown_criterion",
                    case.name,
                    "the baseline records grid-convergence bands for a case whose ladder has "
                    "no third level",
                )
            )
        return
    if not blessed:
        return
    if not isinstance(recorded, Mapping) or not recorded:
        findings.append(
            Finding(
                ERROR,
                "extrapolation_not_blessed",
                case.name,
                f"the ladder {levels} supports a grid-convergence assertion but the baseline "
                "records none for this case; re-run '--bless' so the band is recorded",
            )
        )
        return
    measurements = blessed.get("measurements") or {}
    if not isinstance(measurements, Mapping):
        measurements = {}
    derived = _gci_estimates(case, levels, measurements, reference, [])
    if set(recorded) != set(derived):
        findings.append(
            Finding(
                ERROR,
                "extrapolation_not_derived",
                case.name,
                f"the baseline records grid-convergence bands for {sorted(recorded)} but the "
                f"measurements derive them for {sorted(derived)}",
            )
        )
    for key in sorted(set(recorded) & set(derived)):
        entry = recorded[key] if isinstance(recorded[key], Mapping) else {}
        value, band, basis = derived[key]
        for field, expected in (("value", value), ("tolerance", band)):
            stored = entry.get(field)
            if not isinstance(stored, (int, float)) or isinstance(stored, bool):
                findings.append(
                    Finding(
                        ERROR,
                        "extrapolation_not_derived",
                        case.name,
                        f"'{key}' records {field} {stored!r}, which is not a number",
                        metric=key,
                    )
                )
                continue
            difference = _relative(float(stored) - float(expected), float(expected))
            if abs(difference) > REFERENCE_TOLERANCE:
                findings.append(
                    Finding(
                        ERROR,
                        "extrapolation_not_derived",
                        case.name,
                        (
                            f"'{key}' records {field} {float(stored):.6g} but its own "
                            f"derivation gives {float(expected):.6g} (relative difference "
                            f"{difference:.3e}; {basis}). The recorded grid-convergence band "
                            "is not the one the checker derives; a criterion may not be chosen"
                        ),
                        metric=key,
                    )
                )
    if not findings and isinstance(recorded, Mapping):
        # provenance only: the derivation is reproduced, so nothing to report
        pass


def _check_recorded_criterion(
    case: "GoldenCase",
    blessed: Mapping[str, Any],
    reference: Mapping[str, float],
    declared: Mapping[str, tuple[float, str]],
    findings: list[Finding],
) -> None:
    """The tolerance stored in the baseline must be its own derivation, in every mode.

    ``--quick`` -- the mode the git hook runs -- has no coarse level, so it asserts the
    acceptance against the band the baseline recorded.  That band is therefore an *input a
    hand can edit*: inflate it and the local gate becomes a rubber stamp whose table still
    reads green (reproduced: a 100x Wigley band passes ``--check --quick``).  The stored
    value is never taken on trust.  It is re-derived from the blessed inputs -- the two
    measurements the two-grid bound was taken on, or the band the case declares -- and a
    value that does not reproduce is an error, in quick mode and in full mode alike.  The
    band is derived, not chosen, and the derivation is re-run on every check.
    """
    recorded = blessed.get("tolerances") or {}
    if not isinstance(recorded, Mapping) or not recorded:
        return
    measurements = blessed.get("measurements") or {}
    if not isinstance(measurements, Mapping):
        measurements = {}
    ladder = case.levels()
    if len(ladder) >= 3:
        acceptance, partner = ladder[-2], ladder[-3]
    else:
        acceptance, partner = ladder[-1], (ladder[0] if len(ladder) > 1 else None)
    coarse = (measurements.get(partner) if partner is not None else None) or {}
    fine = measurements.get(acceptance) or {}
    if not isinstance(coarse, Mapping) or not isinstance(fine, Mapping):
        coarse, fine = {}, {}

    for key in sorted(recorded):
        entry = recorded[key]
        entry = entry if isinstance(entry, Mapping) else {}
        stored = entry.get("value")
        basis = str(entry.get("basis") or "")
        if key not in reference:
            findings.append(
                Finding(
                    ERROR,
                    "tolerance_unknown_metric",
                    case.name,
                    (
                        f"the baseline records a tolerance for '{key}', which no reference "
                        "asserts; a band without a criterion to bound is a chosen band"
                    ),
                    metric=key,
                )
            )
            continue
        is_declared = basis.startswith(DECLARED_BASIS_PREFIX)
        if is_declared != (key in declared):
            findings.append(
                Finding(
                    ERROR,
                    "tolerance_origin_mismatch",
                    case.name,
                    (
                        f"'{key}' is recorded with basis {basis!r} but the case "
                        f"{'declares' if key in declared else 'does not declare'} a band for "
                        "it; quick mode and full mode would then assert different criteria"
                    ),
                    metric=key,
                )
            )
            continue
        if is_declared:
            expected, expected_basis = declared[key]
        elif key in coarse and key in fine:
            expected, expected_basis = _tolerance(reference[key], coarse[key], fine[key])
        else:
            findings.append(
                Finding(
                    ERROR,
                    "tolerance_not_derived",
                    case.name,
                    (
                        f"'{key}' records tolerance {stored!r} but the baseline carries no "
                        "coarse/fine measurement to derive it from; the recorded band cannot "
                        "be checked against a two-grid bound"
                    ),
                    metric=key,
                )
            )
            continue
        if not isinstance(stored, (int, float)) or isinstance(stored, bool):
            findings.append(
                Finding(
                    ERROR,
                    "tolerance_not_derived",
                    case.name,
                    f"'{key}' records tolerance {stored!r}, which is not a number",
                    metric=key,
                )
            )
            continue
        difference = _relative(float(stored) - float(expected), float(expected))
        if difference > REFERENCE_TOLERANCE:
            findings.append(
                Finding(
                    ERROR,
                    "tolerance_not_derived",
                    case.name,
                    (
                        f"'{key}' records tolerance {float(stored):.6g} but its own "
                        f"derivation gives {float(expected):.6g} (relative difference "
                        f"{difference:.3e}; {expected_basis}). The recorded band is not the "
                        "one the checker derives, so a hand-edited criterion would decide "
                        "the gate; a criterion may not be chosen"
                    ),
                    metric=key,
                )
            )


def _checked_declarations(
    case: GoldenCase,
    values: Mapping[str, float],
    reference: Mapping[str, float],
    findings: list[Finding],
) -> dict[str, tuple[float, str]]:
    """Validate the bands a case declares; an unjustified declaration is an error.

    Three ways a declaration can be worthless, all fail closed: it names a metric no
    reference asserts, it carries no basis, or it is so wide that the metric is not being
    held to anything (:data:`DECLARED_TOLERANCE_CAP`).
    """
    declared: dict[str, tuple[float, str]] = {}
    for key, entry in case.reference_tolerance(values, reference).items():
        tolerance, basis = entry
        tolerance = float(tolerance)
        basis = str(basis).strip()
        if key not in reference:
            findings.append(
                Finding(
                    ERROR,
                    "declared_tolerance_unknown_metric",
                    case.name,
                    f"a tolerance is declared for '{key}', which no reference asserts",
                    metric=key,
                )
            )
            continue
        if not basis:
            findings.append(
                Finding(
                    ERROR,
                    "declared_tolerance_without_basis",
                    case.name,
                    f"'{key}' declares tolerance {tolerance:g} without naming the contract "
                    "it is quoted from; an unjustified band is a chosen one",
                    metric=key,
                )
            )
            continue
        if not (0.0 < tolerance <= DECLARED_TOLERANCE_CAP):
            findings.append(
                Finding(
                    ERROR,
                    "declared_tolerance_too_loose",
                    case.name,
                    f"'{key}' declares tolerance {tolerance:g}, outside (0, "
                    f"{DECLARED_TOLERANCE_CAP:g}]: a band that wide cannot tell a "
                    "regression from the reference's own precision",
                    metric=key,
                )
            )
            continue
        declared[key] = (tolerance, basis)
    return declared


def evaluate_case(
    case: GoldenCase,
    blessed: Mapping[str, Any] | None,
    *,
    quick: bool,
    findings: list[Finding],
    baseline_optional: bool = False,
) -> CaseResult:
    started = time.perf_counter()
    result = CaseResult(name=case.name, title=case.title, reference_source=case.reference_source)
    ladder = case.levels()
    levels = ladder[:1] if quick else ladder
    try:
        for level in levels:
            result.measurements[level] = {
                key: float(value) for key, value in case.measure(level).items()
            }
        result.reference = {key: float(value) for key, value in case.reference().items()}
        result.geometry_sha256 = case.geometry_sha256()
    except Exception as failure:
        findings.append(
            Finding(ERROR, "case_crashed", case.name, f"{type(failure).__name__}: {failure}")
        )
        result.seconds = time.perf_counter() - started
        return result

    reference = result.reference
    # The acceptance level and its two-grid partner.  A third level in the ladder does not
    # move the acceptance level: it is the extra information the grid-convergence assertion
    # is extrapolated over.  Moving acceptance to the third level would make the cheap hook
    # (coarse measurement only) assert a band its own grid cannot meet, which is how a gate
    # starts lying in the other direction.
    if len(levels) >= 3:
        acceptance, partner = levels[-2], levels[-3]
    else:
        acceptance, partner = levels[-1], (levels[0] if len(levels) > 1 else None)
    finest = result.measurements[acceptance]
    coarsest = result.measurements[partner] if partner is not None else None
    result.invariants = case.invariants(finest)
    result.known_limitations = case.known_limitations(finest)
    declared_tolerances = _checked_declarations(case, finest, reference, findings)

    blessed = dict(blessed or {})
    _check_recorded_criterion(case, blessed, reference, declared_tolerances, findings)
    _check_recorded_extrapolation(case, blessed, reference, findings)
    if not blessed and not baseline_optional:
        findings.append(
            Finding(
                ERROR,
                "case_not_in_baseline",
                case.name,
                "the golden baseline does not know this case; run '--bless' to record it",
            )
        )

    # ---- the stored references must still be what their source computes ---------------
    blessed_reference = {
        str(key): float(value) for key, value in (blessed.get("reference") or {}).items()
    }
    for key, value in reference.items():
        stored = blessed_reference.get(key)
        if stored is None:
            if blessed and not baseline_optional:
                findings.append(
                    Finding(
                        ERROR,
                        "reference_not_in_baseline",
                        case.name,
                        f"reference metric '{key}' is not recorded in the baseline",
                        metric=key,
                    )
                )
            continue
        difference = _relative(stored - value, value)
        if difference > REFERENCE_TOLERANCE:
            findings.append(
                Finding(
                    ERROR,
                    "reference_not_reproducible",
                    case.name,
                    (
                        f"reference '{key}' stored as {stored!r} but {case.reference_source} "
                        f"now yields {value!r} (relative difference {difference:.3e}); "
                        "a reference may not be hand-edited"
                    ),
                    metric=key,
                )
            )
    for key in blessed_reference:
        if key not in reference:
            findings.append(
                Finding(
                    WARNING,
                    "reference_metric_removed",
                    case.name,
                    f"baseline reference '{key}' is no longer produced by the case",
                    metric=key,
                )
            )

    # ---- acceptance: measured fine value against the reference ------------------------
    for key in sorted(reference):
        expected = reference[key]
        blessed_tolerance = (blessed.get("tolerances") or {}).get(key) or {}
        origin = "two-grid"
        if coarsest is None and blessed_tolerance.get("value") is not None:
            # quick mode has no coarse level to re-derive the bound from, so it runs on
            # the tolerance the baseline recorded: same criterion, one measurement.
            tolerance = float(blessed_tolerance["value"])
            basis = (
                f"blessed tolerance {tolerance:.3e} (quick mode: "
                f"{blessed_tolerance.get('basis', 'derived at bless time')})"
            )
            if str(blessed_tolerance.get("basis", "")).startswith(DECLARED_BASIS_PREFIX):
                origin = "declared"
        elif key in declared_tolerances:
            tolerance, declared_basis = declared_tolerances[key]
            basis = f"declared contract: {declared_basis}"
            origin = "declared"
        else:
            tolerance, basis = _tolerance(
                expected, coarsest.get(key) if coarsest else None, finest[key]
            )
            origin = "two-grid"
        if key in declared_tolerances and coarsest is None and blessed_tolerance.get("value"):
            # quick mode runs on the blessed band; keep its recorded provenance, not a
            # re-declaration, so a blessed run and a quick run cannot disagree.
            tolerance = float(blessed_tolerance["value"])
            basis = (
                f"blessed tolerance {tolerance:.3e} (quick mode: "
                f"{blessed_tolerance.get('basis', 'declared at bless time')})"
            )
            origin = "declared"
        error = (finest[key] - expected) / expected if expected else finest[key] - expected
        metric = MetricCheck(
            name=key,
            value=finest[key],
            kind="analytic",
            reference=expected,
            relative_error=error,
            tolerance=tolerance,
            basis=basis,
            origin=origin,
        )
        result.metrics.append(metric)
        if not metric.passed:
            findings.append(
                Finding(
                    ERROR,
                    "reference_mismatch",
                    case.name,
                    (
                        f"{key} = {metric.value:.9g} vs reference {expected:.9g} = "
                        f"{error:+.3%} (tolerance {tolerance:.3%}; {basis})"
                    ),
                    metric=key,
                )
            )

    # ---- grid convergence: the extrapolated value against the same reference ----------
    if not quick and len(levels) >= 3:
        for key, (value, band, basis) in _gci_estimates(
            case, levels, result.measurements, reference, findings
        ).items():
            error = (value - reference[key]) / reference[key] if reference[key] else value - reference[key]
            metric = MetricCheck(
                name=key,
                value=value,
                kind="richardson",
                reference=reference[key],
                relative_error=error,
                tolerance=band,
                basis=basis,
                origin="richardson",
            )
            result.metrics.append(metric)
            if not metric.passed:
                findings.append(
                    Finding(
                        ERROR,
                        "richardson_mismatch",
                        case.name,
                        (
                            f"{key} extrapolated to {value:.9g} vs reference "
                            f"{reference[key]:.9g} = {error:+.3%}, outside the "
                            f"grid-convergence band {band:.3%} ({basis}): the three levels "
                            "converge consistently but not to the analytic value, so a "
                            "discretisation error does not explain the gap"
                        ),
                        metric=key,
                    )
                )

    # ---- invariants: identities over the measured numbers -----------------------------
    for key, (residual, tolerance, basis) in result.invariants.items():
        result.metrics.append(
            MetricCheck(
                name=key,
                value=residual,
                kind="invariant",
                reference=0.0,
                relative_error=abs(residual),
                tolerance=tolerance,
                basis=basis,
                origin="invariant",
            )
        )
        if abs(residual) > tolerance:
            findings.append(
                Finding(
                    ERROR,
                    "invariant_violated",
                    case.name,
                    f"{key}: residual {residual:.3e} > {tolerance:g} ({basis})",
                    metric=key,
                )
            )

    # ---- convergence: refining the grid must move towards the reference ---------------
    if coarsest is not None:
        for key, expected in reference.items():
            coarse_error = abs(coarsest[key] - expected)
            fine_error = abs(finest[key] - expected)
            if fine_error > coarse_error * (1.0 + 1e-9):
                findings.append(
                    Finding(
                        WARNING,
                        "convergence_not_monotone",
                        case.name,
                        (
                            f"{key}: refining the grid moved the value away from the reference "
                            f"({coarse_error:.6g} -> {fine_error:.6g})"
                        ),
                        metric=key,
                    )
                )

    # ---- drift against the blessed measurements --------------------------------------
    blessed_measurements = blessed.get("measurements") or {}
    for level, values in result.measurements.items():
        stored_values = blessed_measurements.get(level)
        if stored_values is None:
            continue
        for key, value in sorted(values.items()):
            if key not in stored_values:
                if blessed:
                    findings.append(
                        Finding(
                            NOTICE,
                            "metric_added",
                            case.name,
                            f"metric '{key}' is new since the baseline was blessed",
                            metric=key,
                        )
                    )
                continue
            stored = float(stored_values[key])
            drift = _relative(value - stored, stored)
            if drift > DRIFT_TOLERANCE:
                findings.append(
                    Finding(
                        ERROR,
                        "drift",
                        case.name,
                        (
                            f"{level}/{key} drifted from {stored!r} to {value!r} "
                            f"({drift:.3e} > {DRIFT_TOLERANCE:g}); blessing the new value "
                            "is a deliberate act, not a side effect of a refactor"
                        ),
                        metric=key,
                    )
                )

    # ---- geometry identity: the mesh itself, not only the numbers integrated on it ----
    # ``geometry_sha256`` is written by ``--bless``; until this block existed it was never
    # read back, so a change that keeps every measured quantity bit-identical
    # (re-triangulation, vertex renumbering, inverted winding -- the surface area takes
    # the norm and the sphere/ellipsoid volume the absolute value) passed silently.
    blessed_fingerprint = str(blessed.get("geometry_sha256") or "")
    measured_fingerprint = result.geometry_sha256
    if not measured_fingerprint:
        findings.append(
            Finding(
                WARNING,
                "geometry_fingerprint_missing",
                case.name,
                "this case does not fingerprint its geometry: only the integrated "
                "quantities are drift-locked",
            )
        )
    elif not blessed_fingerprint:
        if blessed and not baseline_optional:
            # Fail closed: a baseline entry without a fingerprint locks the integrated
            # numbers but not the mesh they were integrated on, so a re-tessellation that
            # keeps every metric would pass. That is an unguarded case, not a warning: a
            # merge that drops the field (observed: the committed baseline carried
            # geometry_sha256 = "" for three of four cases) must turn the gate red.
            findings.append(
                Finding(
                    ERROR,
                    "geometry_fingerprint_not_blessed",
                    case.name,
                    "the baseline records no geometry fingerprint for this case, so its "
                    "mesh is not locked (only its integrated quantities are); "
                    "re-run '--bless' and commit the fingerprint",
                )
            )
    elif measured_fingerprint != blessed_fingerprint:
        findings.append(
            Finding(
                ERROR,
                "geometry_changed",
                case.name,
                (
                    f"geometry fingerprint {blessed_fingerprint[:16]}... -> "
                    f"{measured_fingerprint[:16]}... (the mesh changed; blessing the new "
                    "mesh is a deliberate act, not a side effect of a refactor)"
                ),
            )
        )

    result.seconds = time.perf_counter() - started
    return result


def run_suite(
    *,
    quick: bool = False,
    baseline_path: Path | str = BASELINE,
    cases: Iterable[GoldenCase] | None = None,
    baseline_optional: bool = False,
) -> SuiteResult:
    started = time.perf_counter()
    path = Path(baseline_path)
    result = SuiteResult(baseline=str(path), quick=quick)
    findings: list[Finding] = []
    data: dict[str, Any] = {}
    if not path.is_file():
        if not baseline_optional:
            findings.append(
                Finding(ERROR, "baseline_missing", "*", f"no golden baseline at {path}")
            )
    else:
        try:
            data = json.loads(path.read_text(encoding="utf-8"))
        except json.JSONDecodeError as failure:
            findings.append(
                Finding(ERROR, "baseline_unreadable", "*", f"{path}: {failure}")
            )
            data = {}
        else:
            if data.get("schema") != SCHEMA:
                findings.append(
                    Finding(
                        ERROR,
                        "baseline_schema",
                        "*",
                        f"baseline schema {data.get('schema')!r} != {SCHEMA!r}",
                    )
                )
            forced = [
                entry
                for entry in ((data.get("provenance") or {}).get("forced_bless") or [])
                if isinstance(entry, dict)
            ]
            if forced:
                overridden = sorted(
                    {str(code) for entry in forced for code in entry.get("exceeded") or []}
                )
                commits = ", ".join(
                    str(entry.get("git_commit") or "?")[:12] for entry in forced
                )
                findings.append(
                    Finding(
                        WARNING,
                        "forced_bless_recorded",
                        "*",
                        f"this baseline was written by {len(forced)} forced bless(es) "
                        f"overriding {overridden or 'no boundary'} at {commits}: a forced "
                        "bless is a reviewed decision and stays visible in every run",
                    )
                )
    blessed_cases = {str(key): value for key, value in (data.get("cases") or {}).items()}
    registered = list(cases) if cases is not None else golden_cases()
    registered_names = {case.name for case in registered}
    for stale in sorted(set(blessed_cases) - registered_names):
        findings.append(
            Finding(
                ERROR,
                "baseline_case_unknown",
                stale,
                "the baseline records a case the code no longer registers; "
                "either restore the case or re-bless deliberately",
            )
        )
    for case in registered:
        # Each case keeps the findings it raised. Without this, ``CaseResult.passed`` is
        # vacuously true and the report tells a reader that a failing case passed -- the
        # suite-level flag was right, the per-case column was noise. The same Finding
        # objects stay in the single suite list, so ``SuiteResult.all_findings``
        # (suite-level + per-case) still reports each error exactly once.
        case_findings: list[Finding] = []
        evaluated = evaluate_case(
            case,
            blessed_cases.get(case.name),
            quick=quick,
            findings=case_findings,
            baseline_optional=baseline_optional,
        )
        evaluated.findings = case_findings
        result.cases.append(evaluated)
    result.findings = findings
    result.seconds = time.perf_counter() - started
    return result


# --------------------------------------------------------------------------------------
# blessing and reporting
# --------------------------------------------------------------------------------------
def bless(
    baseline_path: Path | str = BASELINE,
    *,
    force: bool = False,
    cases: Iterable[GoldenCase] | None = None,
) -> tuple[dict[str, Any], SuiteResult]:
    """Measure every case on this tree and write the baseline. Refuses on errors."""
    result = run_suite(
        quick=False, baseline_path=baseline_path, baseline_optional=True, cases=cases
    )
    if not result.passed and not force:
        raise RuntimeError(
            "refusing to bless a tree that fails its own golden gate: "
            + "; ".join(f"[{finding.case}/{finding.code}] {finding.detail}" for finding in result.errors[:5])
        )
    history: list[dict[str, Any]] = []
    try:
        previous = json.loads(Path(baseline_path).read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError, TypeError):
        previous = {}
    for entry in ((previous.get("provenance") or {}).get("forced_bless") or []):
        if isinstance(entry, dict):
            history.append(dict(entry))
    if not result.passed:
        # ``--force`` is the one path that can write a baseline the gate itself refuses.
        # An unrecorded override would make a laundered baseline indistinguishable from a
        # reviewed one, so the override is written into the baseline and re-announced on
        # every later run.  The list is append-only: a later blessing carries it forward,
        # which is what stops a forced bless from being forgotten by a second bless.
        history.append(
            {
                "git_commit": git_commit(),
                "python": platform.python_version(),
                "exceeded": sorted({finding.code for finding in result.errors}),
                "details": [
                    f"[{finding.case}/{finding.code}] {finding.detail}"
                    for finding in result.errors[:10]
                ],
            }
        )

    payloads: dict[str, Any] = {}
    for case in result.cases:
        payloads[case.name] = {
            "title": case.title,
            "reference_source": case.reference_source,
            "geometry_sha256": case.geometry_sha256,
            "reference": {key: _number(value) for key, value in case.reference.items()},
            "measurements": {
                level: {key: _number(value) for key, value in values.items()}
                for level, values in case.measurements.items()
            },
            "tolerances": {
                metric.name: {
                    "value": _number(metric.tolerance) if metric.tolerance is not None else None,
                    "basis": metric.basis,
                }
                for metric in case.metrics
                if metric.kind == "analytic"
            },
            "extrapolated": {
                metric.name: {
                    "value": _number(metric.value),
                    "tolerance": _number(metric.tolerance) if metric.tolerance is not None else None,
                    "basis": metric.basis,
                }
                for metric in case.metrics
                if metric.kind == "richardson"
            },
            "invariants": {
                key: {
                    "residual": _number(residual),
                    "tolerance": tolerance,
                    "basis": basis,
                }
                for key, (residual, tolerance, basis) in case.invariants.items()
            },
            "known_limitations": case.known_limitations,
        }
    baseline = {
        "schema": SCHEMA,
        "provenance": {
            "generated_by": "python -m tensorcad.golden --bless",
            "git_commit": git_commit(),
            "python": platform.python_version(),
            "platform": platform.platform(),
            "geometry_fingerprint": geometry_fingerprint(),
        },
        "policy": {
            "acceptance_tolerance": (
                "two-grid: 2 x |fine - coarse| / |reference|, floor "
                f"{TOLERANCE_FLOOR:g} (relative)"
            ),
            "richardson_tolerance": (
                f"GCI: {GCI_SAFETY_FACTOR:g} x |extrapolated - finest| / |reference|, floor "
                f"{TOLERANCE_FLOOR:g} (relative): on a three-level ladder the extrapolated "
                "value must agree with the analytic reference within the band the "
                "measurements' own convergence justifies"
            ),
            "drift_tolerance": f"{DRIFT_TOLERANCE:g} (relative, vs the blessed measurement)",
            "reference_tolerance": (
                f"{REFERENCE_TOLERANCE:g} (relative): a stored reference must equal what its "
                "declared source computes on every run"
            ),
        },
        "cases": payloads,
    }
    if history:
        baseline["provenance"]["forced_bless"] = history
    Path(baseline_path).parent.mkdir(parents=True, exist_ok=True)
    Path(baseline_path).write_text(json.dumps(baseline, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return baseline, result


def format_table(result: SuiteResult) -> str:
    lines = [
        f"golden geometry regression — {'PASS' if result.passed else 'FAIL'} "
        f"({len(result.cases)} cases, {result.seconds:.2f} s, "
        f"{'quick' if result.quick else 'full'})",
    ]
    header = f"{'case':<17}{'metric':<40}{'measured':>16}{'reference':>16}{'error':>11}{'tol':>10}"
    lines.append(header)
    lines.append("-" * len(header))
    for case in result.cases:
        for metric in case.metrics:
            error = "" if metric.relative_error is None else f"{metric.relative_error:+.2e}"
            tolerance = "" if metric.tolerance is None else f"{metric.tolerance:.2e}"
            reference = "" if metric.reference is None else f"{metric.reference:.9g}"
            flag = "" if metric.passed else "  <-- FAIL"
            label = f"{metric.name} [GCI]" if metric.kind == "richardson" else metric.name
            lines.append(
                f"{case.name:<17}{label:<40}{metric.value:>16.9g}{reference:>16}"
                f"{error:>11}{tolerance:>10}{flag}"
            )
    for limitation in (entry for case in result.cases for entry in case.known_limitations):
        lines.append(f"known limitation ({limitation['id']}): {limitation['detail']}")
    for finding in result.all_findings:
        if finding.severity == NOTICE:
            continue
        lines.append(f"{finding.severity.upper():<8}{finding.case:<17}{finding.code:<30}{finding.detail}")
    return "\n".join(lines)


# --------------------------------------------------------------------------------------
# CLI
# --------------------------------------------------------------------------------------
def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        prog="python -m tensorcad.golden",
        description="Golden geometry regression suite for TensorCAD.",
    )
    action = parser.add_mutually_exclusive_group()
    action.add_argument("--check", action="store_true", help="compare against the baseline (default)")
    action.add_argument("--bless", action="store_true", help="record a new baseline")
    action.add_argument("--list", action="store_true", help="list the golden cases")
    action.add_argument(
        "--staged",
        nargs="?",
        const="",
        default=None,
        metavar="LIST",
        help=(
            "run the check only if the staged paths touch geometry; LIST is a file of "
            "paths (git hook mode), default 'git diff --cached'"
        ),
    )
    parser.add_argument("--quick", action="store_true", help="coarse level only")
    parser.add_argument("--json", metavar="PATH", help="write the machine-readable report here")
    parser.add_argument("--baseline", default=str(BASELINE), help="baseline path")
    parser.add_argument("--force", action="store_true", help="bless even if the gate fails")
    parser.add_argument("--verbose", action="store_true", help="print notice-level findings")
    args = parser.parse_args(argv)

    if args.list:
        for case in golden_cases():
            print(f"{case.name:<17}{case.title}\n{'':<17}reference: {case.reference_source}")
        return 0

    if args.bless:
        baseline, result = bless(args.baseline, force=args.force)
        print(
            f"blessed {len(baseline['cases'])} cases -> {args.baseline} "
            f"({result.seconds:.2f} s)"
        )
        return 0

    if args.staged is not None:
        staged = _staged_paths(args.staged or None)
        if not staged_requires_golden(staged):
            print("golden: no staged geometry change, check skipped")
            return 0
        print(f"golden: {len(staged)} staged path(s) touch geometry, running the gate")

    result = run_suite(quick=args.quick, baseline_path=args.baseline)
    print(format_table(result))
    if args.verbose:
        for finding in result.all_findings:
            if finding.severity == NOTICE:
                print(f"NOTICE  {finding.case:<17}{finding.code:<30}{finding.detail}")
    if args.json:
        Path(args.json).write_text(
            json.dumps(result.to_dict(), indent=2, sort_keys=True) + "\n", encoding="utf-8"
        )
        print(f"report -> {args.json}")
    return 0 if result.passed else 1


if __name__ == "__main__":  # pragma: no cover
    sys.exit(main())
