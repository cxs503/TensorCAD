"""File-frame evidence and normalisation for *imported* hull geometry.

A benchmark hull is published in whatever frame its originator used, and the two
conventions seen in the fleet's own archives are not the same:

* the KCS IGES (NMRI/G2015) is keel-origin, ``z`` up, ``x`` running aft->fwd from the
  aft perpendicular, so the half-breadth field route reads it as-is;
* the KVLCC2 IGES (same publisher, different export) puts ``z = 0`` at the design
  waterline with ``z`` growing *downwards* (measured extent ``0 ... -28.067 m``) and
  puts ``x = 0`` near the forward end with the stern at ``x ~= -328 m``.

Feeding the second file to a route that keeps samples with ``0 <= z <= draft`` and
places the grid on ``0 <= x <= Lpp`` keeps *nothing*: the waterline slices vanish and
the benchmark reports "fewer than two usable waterline slices".  That is a frame
defect, not a shape defect, and the honest fix is to *measure* the file's frame, apply
an explicit transformation, and record the evidence next to the result.

This module holds only the measurement and the transformation, as pure functions over
a point cloud, so the policy is testable without a 6 MB IGES file and without a mesh:

* :func:`detect` reports the axis ranges, the fraction of samples below the file's own
  ``z = 0`` plane and the fraction with negative ``y``;
* :func:`normalise` applies a *declarative* :class:`FramePolicy` (``z_axis``,
  ``x_origin``) and returns the transformed points together with the evidence and the
  exact ``x``-shift / ``z``-flip it applied, so a build can be re-checked from its
  recorded numbers.

The conventions after normalisation are the ones the mesh routes assume:
``x`` from 0 at the aft-most point of the model, ``y`` to starboard, ``z`` up from the
design waterline, metres.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Iterable, Mapping, Sequence

#: A file counts as ``z``-down when the part of the model *below* its own ``z = 0`` plane
#: is more than this many times the part above it.  The ratio, not the sample fraction, is
#: the criterion: the published KVLCC2 IGES has only ~37 % of its samples below z = 0 at
#: coarse sampling (state of the art: ``z_max = +0.055 m``, ``z_min = -28.067 m``) because
#: the above-water sliver is sampled densely per patch, so a fraction test would call a
#: z-down export "up" and keep nothing.  A keel-origin export has essentially nothing
#: below 0 and is unaffected by the threshold.
Z_DOWN_RATIO = 3.0

#: Tolerance below which an ``x`` translation is *not* recorded (the file already starts
#: at its own origin, and a 1e-9 m shift is round-off, not a policy).
X_TOLERANCE = 1e-6

#: When a file's own ``x = 0`` sits *inside* its model, the share of the extent that lies
#: aft of the origin decides whether that origin is the aft perpendicular or the forward
#: end of a stern-negative export.  The two published exports are unambiguous: KCS keeps
#: 2.5 % of its 243.8 m extent aft of ``x = 0`` (the transom overhang; the origin *is* the
#: AP), KVLCC2 keeps 98.3 % of its 333.6 m extent aft of it (the origin is at the bow).
#: Any split between the two decides the same way; 0.5 is used because it is the point
#: where "the origin is nearer the forward end" becomes literally true.
X_ORIGIN_SPLIT = 0.5

#: Frame string written next to a normalised build.
NORMALISED_FRAME = "x from 0 at the aft-most point, y to starboard, z up from the design waterline (m)"


@dataclass(frozen=True)
class FrameEvidence:
    """What the sampled surface says about the file's frame."""

    x_min: float
    x_max: float
    y_min: float
    y_max: float
    z_min: float
    z_max: float
    samples: int
    z_below_zero_fraction: float
    y_negative_fraction: float

    @property
    def below_zero(self) -> float:
        """Depth of the model below the file's own ``z = 0`` plane (m, >= 0)."""
        return max(0.0, -self.z_min)

    @property
    def above_zero(self) -> float:
        """Height of the model above the file's own ``z = 0`` plane (m, >= 0)."""
        return max(0.0, self.z_max)

    @property
    def z_axis(self) -> str:
        """``"down"`` when the model hangs below its own z = 0 plane, else ``"up"``."""
        if self.below_zero <= 1e-9:
            return "up"
        return "down" if self.below_zero > Z_DOWN_RATIO * self.above_zero else "up"

    @property
    def x_extent(self) -> float:
        return self.x_max - self.x_min

    @property
    def depth(self) -> float:
        """Vertical extent of the model as published (m)."""
        return self.z_max - self.z_min

    @property
    def full_beam(self) -> float:
        """Beam implied by the file: doubled when the file carries one side only."""
        if self.y_negative_fraction < 0.01 or self.y_min >= 0.0:
            return 2.0 * max(abs(self.y_min), abs(self.y_max))
        return self.y_max - self.y_min

    @property
    def half_model(self) -> bool:
        """True when the file carries (essentially) only ``y >= 0``: a mirrored half hull."""
        return self.y_negative_fraction < 0.01

    def as_dict(self) -> dict[str, Any]:
        return {
            "samples": self.samples,
            "x_min": round(self.x_min, 4),
            "x_max": round(self.x_max, 4),
            "y_min": round(self.y_min, 4),
            "y_max": round(self.y_max, 4),
            "z_min": round(self.z_min, 4),
            "z_max": round(self.z_max, 4),
            "x_extent": round(self.x_extent, 4),
            "depth": round(self.depth, 4),
            "z_axis": self.z_axis,
            "below_zero": round(self.below_zero, 4),
            "above_zero": round(self.above_zero, 4),
            "z_below_zero_fraction": round(self.z_below_zero_fraction, 6),
            "y_negative_fraction": round(self.y_negative_fraction, 6),
            "half_model": self.half_model,
        }


@dataclass(frozen=True)
class FramePolicy:
    """Declarative frame transformation, recorded per ship in the manifest.

    ``z_axis``: ``"auto"`` (use :func:`detect`), ``"up"`` or ``"down"`` (force).
    ``x_origin``: ``"auto"`` (translate aft-most point to 0) or ``"as_is"``.
    """

    z_axis: str = "auto"
    x_origin: str = "auto"

    def __post_init__(self) -> None:
        if self.z_axis not in ("auto", "up", "down"):
            raise ValueError(f"frame z_axis must be auto|up|down, got {self.z_axis!r}")
        if self.x_origin not in ("auto", "as_is"):
            raise ValueError(f"frame x_origin must be auto|as_is, got {self.x_origin!r}")

    @classmethod
    def from_mapping(cls, payload: Mapping[str, Any] | "FramePolicy" | None) -> "FramePolicy":
        if payload is None:
            return cls()
        if isinstance(payload, FramePolicy):
            return payload
        return cls(z_axis=str(payload.get("z_axis", "auto")), x_origin=str(payload.get("x_origin", "auto")))

    def as_dict(self) -> dict[str, str]:
        return {"z_axis": self.z_axis, "x_origin": self.x_origin}


def _points(points: Iterable[Sequence[float]]) -> list[tuple[float, float, float]]:
    out = [(float(p[0]), float(p[1]), float(p[2])) for p in points]
    if not out:
        raise ValueError("no points to inspect")
    return out


def detect(points: Iterable[Sequence[float]]) -> FrameEvidence:
    """Measure the frame of a sampled surface."""
    cloud = _points(points)
    n = len(cloud)
    xs = [p[0] for p in cloud]
    ys = [p[1] for p in cloud]
    zs = [p[2] for p in cloud]
    below = sum(1 for z in zs if z < 0.0) / n
    negatives = sum(1 for y in ys if y < 0.0) / n
    return FrameEvidence(
        x_min=min(xs),
        x_max=max(xs),
        y_min=min(ys),
        y_max=max(ys),
        z_min=min(zs),
        z_max=max(zs),
        samples=n,
        z_below_zero_fraction=below,
        y_negative_fraction=negatives,
    )


@dataclass(frozen=True)
class FrameReport:
    """Evidence + the transformation that was actually applied."""

    policy: FramePolicy
    evidence: FrameEvidence
    z_flip: bool
    x_shift: float
    x_origin_reason: str = ""

    @property
    def frame(self) -> str:
        return NORMALISED_FRAME

    def as_dict(self) -> dict[str, Any]:
        return {
            "policy": self.policy.as_dict(),
            "evidence": self.evidence.as_dict(),
            "z_flip": self.z_flip,
            "x_shift": round(self.x_shift, 6),
            "x_origin_reason": self.x_origin_reason,
            "frame": self.frame,
            "note": (
                "detected from the sampled surface, not assumed: a z-down, stern-negative file "
                "is reflected and translated before it reaches the half-breadth field"
            ),
        }


def _x_origin_decision(evidence: FrameEvidence) -> tuple[float, str]:
    """Decide the ``x`` translation for the ``auto`` policy, and say why.

    ``x = 0`` is the fleet frame's aft end, so the question is *where the file's own origin
    is*, not simply *where the model starts*.  Always translating the aft-most point to 0
    is wrong for a file that is already AP-origin: the published KCS model runs from
    ``-6.0`` (transom overhang) to ``237.84`` m (bulb) around an origin that *is* the aft
    perpendicular, so a 6 m shift slides the hull into the Lpp window and the forebody is
    cut off by it -- measured length 227 m against a declared 230 m, which the fleet's
    principal-dimension gate refuses.  The three cases are therefore treated separately.
    """
    extent = evidence.x_extent
    if evidence.x_min >= -X_TOLERANCE:
        # Nothing aft of the file's own origin: the model is authored forward of it, so the
        # aft-most point has to be brought to 0 or the Lpp window would start mid-hull.
        shift = -evidence.x_min
        return shift, (
            f"model starts forward of its own x = 0 (x_min = {evidence.x_min:.3f} m); "
            f"translated by {shift:+.3f} m so the aft-most point is at x = 0"
        )
    if evidence.x_max <= X_TOLERANCE or extent <= 0.0:
        shift = -evidence.x_min
        return shift, (
            f"stern-negative export (x in [{evidence.x_min:.3f}, {evidence.x_max:.3f}] m, fully "
            f"aft of its own origin); translated by {shift:+.3f} m"
        )
    aft_fraction = abs(evidence.x_min) / extent
    if aft_fraction > X_ORIGIN_SPLIT:
        shift = -evidence.x_min
        return shift, (
            f"origin sits near the forward end ({aft_fraction:.1%} of the {extent:.1f} m extent is "
            f"aft of x = 0); translated by {shift:+.3f} m"
        )
    return 0.0, (
        f"origin is the aft perpendicular: only {aft_fraction:.1%} of the {extent:.1f} m extent "
        f"({abs(evidence.x_min):.3f} m) is aft of x = 0, so no translation is applied"
    )


def normalise(
    points: Iterable[Sequence[float]],
    policy: Mapping[str, Any] | FramePolicy | None = None,
) -> tuple[list[tuple[float, float, float]], FrameReport]:
    """Reflect/translate a point cloud into the fleet frame; return points and the report."""
    pol = FramePolicy.from_mapping(policy)
    cloud = _points(points)
    evidence = detect(cloud)
    z_axis = evidence.z_axis if pol.z_axis == "auto" else pol.z_axis
    z_flip = z_axis == "down"
    if pol.x_origin == "auto":
        x_shift, reason = _x_origin_decision(evidence)
        if abs(x_shift) <= X_TOLERANCE:
            x_shift = 0.0
    else:
        x_shift, reason = 0.0, "policy x_origin = 'as_is': the file's own origin is kept"
    out = [(x + x_shift, y, (-z if z_flip else z)) for x, y, z in cloud]
    return out, FrameReport(policy=pol, evidence=evidence, z_flip=z_flip, x_shift=x_shift,
                            x_origin_reason=reason)


__all__ = [
    "NORMALISED_FRAME",
    "X_ORIGIN_SPLIT",
    "X_TOLERANCE",
    "Z_DOWN_RATIO",
    "FrameEvidence",
    "FramePolicy",
    "FrameReport",
    "detect",
    "normalise",
]