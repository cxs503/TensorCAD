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
   baseline knows and the code forgot, a missing npz — all are errors, never "skipped".

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

    def known_limitations(self, values: Mapping[str, float]) -> list[dict[str, Any]]:
        return []

    def geometry_sha256(self) -> str:
        return ""


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

    def _parameters(self, level: str) -> Any:
        from . import wigley

        nx, nz = {"coarse": (41, 13), "fine": (161, 41)}[level]
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
        from . import wigley

        return _mesh_sha256(wigley.generate_wigley(self._parameters("fine")).mesh)


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

    def _parameters(self, level: str) -> Any:
        from . import series60

        nx, nz = {"coarse": (81, 25), "fine": (161, 49)}[level]
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


class SphereEllipsoidGolden(GoldenCase):
    """SDF primitives through marching tetrahedra: the closed forms are exact."""

    name = "sphere_ellipsoid"
    title = "Sphere (r=1) and ellipsoid (2.0, 1.5, 0.5): SDF -> marching tetrahedra"
    reference_source = "(4/3)πr³, 4πr² and (4/3)πabc over the exact primitive fields"
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
        reference = {
            "sphere_volume_m3": 4.0 / 3.0 * math.pi * self.RADIUS**3,
            "sphere_surface_area_m2": 4.0 * math.pi * self.RADIUS**2,
            "ellipsoid_volume_m3": 4.0 / 3.0 * math.pi * self.AXES[0] * self.AXES[1] * self.AXES[2],
        }
        for axis, extent in zip("xyz", (2 * self.RADIUS,) * 3, strict=True):
            reference[f"sphere_extent_{axis}_m"] = extent
        for axis, extent in zip("xyz", (2.0 * axis for axis in self.AXES), strict=True):
            reference[f"ellipsoid_extent_{axis}_m"] = extent
        return reference

    def known_limitations(self, values: Mapping[str, float]) -> list[dict[str, Any]]:
        return [
            {
                "id": "ellipsoid-area-not-asserted",
                "detail": (
                    "the ellipsoid surface area has no elementary closed form (elliptic "
                    "integrals), so it is drift-locked but not asserted against a reference"
                ),
                "measured_surface_area_m2": _number(values["ellipsoid_surface_area_m2"]),
                "asserted": False,
            }
        ]


class DTMB5415Golden(GoldenCase):
    """The shipped reduced model of the DTMB 5415 design dataset.

    The reference is the dataset itself (``docs/data/dtmb5415/holdout_fixture.npz``): the
    holdout residual the dataset records for one configuration. Absolute displacement has
    no published counterpart inside the repository, so it is drift-locked only, and the
    surface-integral/cloud-integral disagreement is recorded instead of dressed up as an
    acceptance.
    """

    name = "dtmb5415"
    title = "DTMB 5415 baseline hull: shipped rank-40 reduced model (npz, no 1 GB source)"
    reference_source = "docs/data/dtmb5415/holdout_fixture.npz:per_config_rms_m[0]"
    CONFIG = 0
    MIRROR_AXIS = 1

    def _fixture(self) -> Any:
        import numpy as np

        return np.load(ROOT / "docs" / "data" / "dtmb5415" / "holdout_fixture.npz")

    def _space(self) -> Any:
        from .dtmb5415 import DTMB5415Space

        return DTMB5415Space()

    def measure(self, level: str) -> dict[str, float]:
        import numpy as np

        from .dtmb5415 import WATERLINE_Z
        from .hydrostatics import compute_hydrostatics

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
            "surface_integral_volume_m3": float(
                compute_hydrostatics(mesh, waterline_z=WATERLINE_Z, closed=True).displacement_volume
            ),
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
        fixture = self._fixture()
        return {"holdout_rms_m": float(fixture["per_config_rms_m"][self.CONFIG])}

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
        return invariants

    def known_limitations(self, values: Mapping[str, float]) -> list[dict[str, Any]]:
        cloud_volume = values["displacement_volume_m3"]
        surface_volume = values["surface_integral_volume_m3"]
        return [
            {
                "id": "surface-integral-volume-disagrees-with-cloud-integral",
                "detail": (
                    "compute_hydrostatics over space.mesh() does not reproduce the cloud "
                    "integral reported by space.hydrostatics(); recorded as a diagnostic, "
                    "not asserted against the dataset"
                ),
                "cloud_volume_m3": _number(cloud_volume),
                "surface_integral_volume_m3": _number(surface_volume),
                "relative_discrepancy": _number(_relative(surface_volume - cloud_volume, cloud_volume)),
                "asserted": False,
            }
        ]


def golden_cases() -> list[GoldenCase]:
    return [WigleyGolden(), Series60Golden(), SphereEllipsoidGolden(), DTMB5415Golden()]


# --------------------------------------------------------------------------------------
# helpers
# --------------------------------------------------------------------------------------
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


def geometry_files(root: Path = ROOT) -> list[Path]:
    """Every file whose change obliges a golden run."""
    files = sorted((root / "tensorcad").rglob("*.py"))
    files.append(root / "tests" / "test_golden_regression.py")
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
    levels = LEVELS[:1] if quick else LEVELS
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
    finest = result.measurements[levels[-1]]
    coarsest = result.measurements[levels[0]] if len(levels) > 1 else None
    result.invariants = case.invariants(finest)
    result.known_limitations = case.known_limitations(finest)

    blessed = dict(blessed or {})
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
        if coarsest is None and blessed_tolerance.get("value") is not None:
            # quick mode has no coarse level to re-derive the bound from, so it runs on
            # the tolerance the baseline recorded: same criterion, one measurement.
            tolerance = float(blessed_tolerance["value"])
            basis = (
                f"blessed tolerance {tolerance:.3e} (quick mode: "
                f"{blessed_tolerance.get('basis', 'derived at bless time')})"
            )
        else:
            tolerance, basis = _tolerance(
                expected, coarsest.get(key) if coarsest else None, finest[key]
            )
        error = (finest[key] - expected) / expected if expected else finest[key] - expected
        metric = MetricCheck(
            name=key,
            value=finest[key],
            kind="analytic",
            reference=expected,
            relative_error=error,
            tolerance=tolerance,
            basis=basis,
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
        result.cases.append(
            evaluate_case(
                case,
                blessed_cases.get(case.name),
                quick=quick,
                findings=findings,
                baseline_optional=baseline_optional,
            )
        )
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
            "drift_tolerance": f"{DRIFT_TOLERANCE:g} (relative, vs the blessed measurement)",
            "reference_tolerance": (
                f"{REFERENCE_TOLERANCE:g} (relative): a stored reference must equal what its "
                "declared source computes on every run"
            ),
        },
        "cases": payloads,
    }
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
            lines.append(
                f"{case.name:<17}{metric.name:<40}{metric.value:>16.9g}{reference:>16}"
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
