"""Fifth golden case: the attitude module, measured at combined heel and trim.

Wigley, Series 60, the SDF primitives and DTMB 5415 are all measured *upright and level*.
A regression in the attitude path -- the waterline normal, the sign of the heel, the
clipping of the box faces, the in-plane factor of the waterplane -- moves nothing that
those four cases integrate, so before this file the golden gate was blind to it.

The case is deliberately small: a barge ``L x B x depth`` floating at ``V = L·B·T`` with
``heel = 15 deg`` and ``trim = 3 deg``.  Its measured side is
:func:`tensorcad.box_inclined.box_inclined_hydrostatics` (clip the six faces, close the
cut with the cap polygon, integrate by the divergence theorem).  Its reference side is
:mod:`tools.golden.inclined_box_oracle`, which integrates the submerged column
``0 <= z <= z_w(x, y)`` analytically in ``z`` and by Simpson in the two horizontal
directions: no face, no cap polygon, no triangle moment.  Two structurally different
derivations of one number, so the criterion is dual-derived and no number is transcribed.

The level axis is the solve tolerance of the waterline-offset bisection, not a grid: both
routes are exact to round-off at every waterline, so what moves between ``coarse`` and
``fine`` is only how tightly the floating position is held, and the acceptance band is the
two-grid change of that (2.7e-4).  This is a *modelling* lock -- it pins the attitude
conventions and the clip geometry -- not a convergence study; the mesh cases remain the
convergence ones.
"""
from __future__ import annotations

import hashlib
import math
from pathlib import Path
from typing import Any, ClassVar, Mapping

from .golden import GoldenCase

__all__ = ["InclinedBoxGolden"]


class InclinedBoxGolden(GoldenCase):
    """Barge at combined heel + trim: exact polyhedron clip vs independent column integral."""

    name = "inclined_box"
    title = (
        "Barge L=10 B=4 depth=3, V=L·B·T=60 m3, heel=15 deg trim=3 deg: "
        "polyhedron clip vs column integral"
    )
    reference_source = (
        "tools/golden/inclined_box_oracle.oracle_hydrostatics: submerged column z_w(x,y) "
        "integrated analytically in z and by composite Simpson (1200 intervals/axis) in "
        "x and y -- no face clip, no cap polygon"
    )

    LENGTH: ClassVar[float] = 10.0
    BEAM: ClassVar[float] = 4.0
    DEPTH: ClassVar[float] = 3.0
    #: Design draft, so the displaced volume is ``L·B·T`` -- a number the case cannot move
    #: by changing attitude, which is what makes the displacement identity an identity.
    DRAFT: ClassVar[float] = 1.5
    HEEL_DEG: ClassVar[float] = 15.0
    TRIM_DEG: ClassVar[float] = 3.0
    #: Simpson intervals per horizontal axis on the reference route.
    RESOLUTION: ClassVar[int] = 1200
    #: Solve tolerance of the waterline-offset bisection (relative to the model scale) per
    #: level.  This *is* the level axis: the two routes are exact at every waterline, so the
    #: only thing that moves between the levels is how tightly the floating position is
    #: solved.  The gate therefore derives the band from the coarse-to-fine change it
    #: records (2.7e-4 here) instead of the case declaring one.
    SOLVE_TOLERANCE: ClassVar[dict[str, float]] = {"coarse": 3e-5, "fine": 1e-13}

    METRICS: ClassVar[dict[str, str]] = {
        "displacement_volume_m3": "displacement_volume",
        "lcb_x_m": "longitudinal_centre_of_buoyancy",
        "tcb_y_m": "transverse_centre_of_buoyancy",
        "vcb_z_m": "vertical_centre_of_buoyancy",
        "waterplane_area_m2": "waterplane_area",
        "waterline_offset_m": "waterline_offset",
    }

    @property
    def displaced_volume(self) -> float:
        return self.LENGTH * self.BEAM * self.DRAFT

    def levels(self) -> tuple[str, ...]:
        return ("coarse", "fine")

    def measure(self, level: str) -> dict[str, float]:
        from .box_inclined import box_inclined_hydrostatics

        hydro = box_inclined_hydrostatics(
            self.LENGTH,
            self.BEAM,
            self.DEPTH,
            volume_m3=self.displaced_volume,
            heel_deg=self.HEEL_DEG,
            trim_deg=self.TRIM_DEG,
            tolerance=self.SOLVE_TOLERANCE[level],
        )
        x, y, z = hydro.centre_of_buoyancy
        return {
            "displacement_volume_m3": float(hydro.volume_m3),
            "lcb_x_m": float(x),
            "tcb_y_m": float(y),
            "vcb_z_m": float(z),
            "waterplane_area_m2": float(hydro.waterplane_area_m2),
            "waterline_offset_m": float(hydro.waterline_offset),
        }

    def reference(self) -> dict[str, float]:
        """The same six numbers from the column-integration oracle (imported by path)."""
        oracle = _oracle_module()
        values = oracle.oracle_hydrostatics(
            self.LENGTH,
            self.BEAM,
            self.DEPTH,
            volume_m3=self.displaced_volume,
            heel_deg=self.HEEL_DEG,
            trim_deg=self.TRIM_DEG,
            resolution=self.RESOLUTION,
        )
        return {key: float(values[key]) for key in self.METRICS}

    def invariants(self, values: Mapping[str, float]) -> dict[str, tuple[float, float, str]]:
        """Identities that involve no reference at all.

        * the displaced volume must be the requested displacement -- the solve holds the
          displacement and moves the waterline, so an attitude change may not move ``V``;
        * mirroring the heel must mirror the transverse centre and leave the rest alone,
          which is the sign convention of the waterline normal, checked on the module.
        """
        from .box_inclined import box_inclined_hydrostatics

        residual_volume = values["displacement_volume_m3"] - self.displaced_volume
        mirrored = box_inclined_hydrostatics(
            self.LENGTH,
            self.BEAM,
            self.DEPTH,
            volume_m3=self.displaced_volume,
            heel_deg=-self.HEEL_DEG,
            trim_deg=self.TRIM_DEG,
        )
        mx, my, mz = mirrored.centre_of_buoyancy
        residual_mirror = abs(values["tcb_y_m"] + my)
        return {
            "identity_volume_equals_requested_displacement": (
                residual_volume,
                1e-6,
                "V = L·B·T is the constraint the waterline offset is solved for, so it does "
                "not depend on the attitude",
            ),
            "identity_mirrored_heel_mirrors_transverse_centre": (
                residual_mirror,
                1e-8,
                "heel -> -heel must give y_B -> -y_B (body y is starboard, heel>0 is "
                "starboard down)",
            ),
        }

    def geometry_sha256(self) -> str:
        """Fingerprint of the clipped polyhedron itself, not only of the integrals.

        The volumes of a polyhedron are insensitive to a face that is retriangulated, an
        ordering that changes, or a convention that flips a vertex duplication count, so
        the fingerprint is taken over the rounded vertices and the face count of the
        submerged solid.
        """
        from .box_inclined import submerged_box
        from .box_inclined import waterline_normal

        normal = waterline_normal(self.HEEL_DEG, self.TRIM_DEG)
        # The offset is the solved one: read it off the acceptance measurement so the
        # fingerprint is taken on the same solid the acceptance level integrates.
        state = self.measure(self.levels()[-1])
        hydro = submerged_box(
            self.LENGTH,
            self.BEAM,
            self.DEPTH,
            heel_deg=self.HEEL_DEG,
            trim_deg=self.TRIM_DEG,
            waterline_offset=state["waterline_offset_m"],
        )
        points = sorted(
            tuple(round(component, 12) for component in point)
            for point in hydro.waterplane_polygon
        )
        digest = hashlib.sha256()
        digest.update(f"{self.LENGTH!r}|{self.BEAM!r}|{self.DEPTH!r}".encode())
        digest.update(f"|{self.HEEL_DEG!r}|{self.TRIM_DEG!r}".encode())
        digest.update(f"|n={tuple(round(c, 15) for c in normal)!r}".encode())
        digest.update(f"|faces={hydro.faces}|vertices={hydro.vertices}".encode())
        for point in points:
            digest.update(repr(point).encode())
        return digest.hexdigest()


def _oracle_module() -> Any:
    """Import the repository's ``tools/golden/inclined_box_oracle.py`` by path.

    ``tools/`` is not a package and must not become one (nothing in ``tensorcad`` may gain
    an import dependency on the tooling), so the oracle is loaded from its file.  The
    reference route is therefore *outside* the library it judges: a change to
    ``tensorcad`` cannot move the number the case asserts against.
    """
    import importlib.util

    root = Path(__file__).resolve().parents[1]
    path = root / "tools" / "golden" / "inclined_box_oracle.py"
    if not path.is_file():  # pragma: no cover - broken checkout
        raise FileNotFoundError(f"the inclined-box oracle is missing: {path}")
    spec = importlib.util.spec_from_file_location("tensorcad_inclined_box_oracle", path)
    if spec is None or spec.loader is None:  # pragma: no cover - broken checkout
        raise ImportError(f"cannot load the inclined-box oracle from {path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module