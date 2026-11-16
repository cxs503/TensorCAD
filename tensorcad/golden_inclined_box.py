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

    #: Relative round-off floor added to every invariant band: the integrals are sums of a
    #: few hundred triangle moments of magnitude ``V``, so their conditioning floor is far
    #: above one machine epsilon, and a band without it would fail on exact arithmetic.
    ROUNDOFF_FLOOR: ClassVar[float] = 1e-9

    #: Two-sided safety factor: the offset error may push the waterline either way and the
    #: cap polygon's area is not exactly the reference waterplane area.
    SAFETY: ClassVar[float] = 2.0

    @property
    def scale(self) -> float:
        """The bisection's model scale, as :func:`box_inclined.box_inclined_hydrostatics` uses it."""
        return max(self.LENGTH, self.BEAM, self.DEPTH)

    def offset_half_width(self, level: str) -> float:
        """How far the solved waterline offset may sit from the one that floats ``V`` exactly.

        The bisection stops when the bracket is narrower than ``tolerance x scale``, so the
        returned offset is within half that width of the root -- this is the only inexact
        step in the case, which is why it *is* the level axis.
        """
        return 0.5 * self.SOLVE_TOLERANCE[level] * self.scale

    def invariants(
        self, values: Mapping[str, float], *, level: str | None = None
    ) -> dict[str, tuple[float, float, str]]:
        """Identities that involve no reference at all, bounded by the level's solve accuracy.

        * the displaced volume must be the requested displacement -- the solve holds the
          displacement and moves the waterline, so an attitude change may not move ``V``;
        * mirroring the heel must mirror the transverse centre and leave the rest alone,
          which is the sign convention of the waterline normal, checked on the module.

        Both bands are *derived from the level*, never typed: the bisection holds the offset
        to ``offset_half_width(level)`` and ``dV/dc`` is the waterplane area, so the volume
        identity cannot be tighter than ``A_wp x half-width``; the transverse centre moves by
        at most ``(B/2) x dV / V``.  A single band for both rungs was the round-1 defect: the
        1e-6 m3 the fine level meets to 1e-5 of its band is 2600x below what the coarse solve
        can deliver, so ``--quick`` (the git hook) reported two ERRORs on every tree the full
        suite passed, and the hook's only escape was ``TENSORCAD_GOLDEN_SKIP=1``.
        """
        from .box_inclined import box_inclined_hydrostatics

        raised = level or self.levels()[-1]
        half_width = self.offset_half_width(raised)
        waterplane = values["waterplane_area_m2"]
        volume_core = self.SAFETY * waterplane * half_width
        volume_band = volume_core + self.ROUNDOFF_FLOOR * self.displaced_volume
        mirror_core = self.SAFETY * (self.BEAM / 2.0) * volume_core / self.displaced_volume
        mirror_band = mirror_core + self.ROUNDOFF_FLOOR * (self.BEAM / 2.0)

        residual_volume = values["displacement_volume_m3"] - self.displaced_volume
        # The mirror route is solved at the *same* tolerance as the measurement it is
        # compared with: comparing a coarse solve against an exact one reports the coarse
        # solve's offset error, which is what the volume identity already asserts.
        mirrored = box_inclined_hydrostatics(
            self.LENGTH,
            self.BEAM,
            self.DEPTH,
            volume_m3=self.displaced_volume,
            heel_deg=-self.HEEL_DEG,
            trim_deg=self.TRIM_DEG,
            tolerance=self.SOLVE_TOLERANCE[raised],
        )
        mx, my, mz = mirrored.centre_of_buoyancy
        residual_mirror = abs(values["tcb_y_m"] + my)
        criterion = (
            f"bisection half-width at '{raised}' = 0.5 x {self.SOLVE_TOLERANCE[raised]:g} x "
            f"{self.scale:g} m = {half_width:.3e} m"
        )
        return {
            "identity_volume_equals_requested_displacement": (
                residual_volume,
                volume_band,
                "V = L·B·T is the constraint the waterline offset is solved for, so it does "
                "not depend on the attitude: |dV| <= 2 x A_wp x half-width + "
                f"{self.ROUNDOFF_FLOOR:g} x V; {criterion}, A_wp = {waterplane:.6g} m2 "
                f"-> {volume_band:.3e} m3",
            ),
            "identity_mirrored_heel_mirrors_transverse_centre": (
                residual_mirror,
                mirror_band,
                "heel -> -heel must give y_B -> -y_B (body y is starboard, heel>0 is "
                "starboard down), both routes solved to the same tolerance: |dy_B| <= "
                "2 x (B/2) x dV / V + "
                f"{self.ROUNDOFF_FLOOR:g} x B/2; {criterion} -> {mirror_band:.3e} m",
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