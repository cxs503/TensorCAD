"""DARPA SUBOFF bare hull as an analytic geometry field.

Why this hull: it is the one submarine form with published, tabulated geometry and a
long history of CFD comparisons, so it is the natural first *non-trivial* shape for
the Continuous Geometry Core -- axisymmetric, slender (L/D = 8.57), with a real bow,
parallel midbody, stern taper and stern cap rather than a sphere.

Geometry source
---------------
Groves, Huang & Chang (1989), *Geometric Characteristics of DARPA SUBOFF Models*,
DTRC/SHD-1298-01: the four-segment offset-table equations (bow, parallel midbody,
stern taper, stern cap) expressed here in feet, exactly as transcribed by TensorLBM's
``src/tensorlbm/suboff_cad.py::suboff_radius_profile`` (read 2026-10-24). The
transcription is verified against that module in ``tests/test_suboff.py`` by pinning
its published values.

Units and frame
---------------
The field's own length unit is whatever ``units`` says (``"m"`` by default). The
x axis runs **bow (x = 0) to stern (x = length)**; y and z are the transverse and
vertical offsets from the centreline. That is the ship convention used by the rest
of the Tensor ecosystem (and by TensorLBM's mask builder, which additionally places
the hull inside a lattice box -- that translation is a simulation setup, not a
property of the shape, so it is not baked in here).

The signed distance is the true distance to the hull surface: for an axisymmetric
body it equals the distance in the (x, rho) half-plane to the meridian curve, so it
is exact rather than a radial gap scaled by a slope factor. ``tests/test_suboff.py``
checks it against a brute-force dense meridian to ``1e-3 * R``.
"""

from __future__ import annotations

import math
from collections.abc import Callable, Iterable, Sequence
from dataclasses import dataclass
from typing import Any, ClassVar

from .geometry import BoundingBox
from .geometry_field import (
    DEFAULT_FRAME,
    DEFAULT_UNITS,
    SIGNED_DISTANCE,
    FieldError,
    GeometryField,
    Point,
    _points,
)

__all__ = [
    "SUBOFF_LENGTH_M",
    "SUBOFF_RADIUS_M",
    "SuboffHullField",
    "suboff_dimensions",
    "suboff_radius_fraction",
    "suboff_reference_geometry",
]

# --------------------------------------------------------------------------- #
# The DARPA offset table, in the feet its source is published in
# --------------------------------------------------------------------------- #

#: Total hull length, ft (Groves et al. 1989; 14.291667 ft = 4.3561 m).
_LENGTH_FT = 14.291667
#: Maximum hull radius, ft (0.8333333 ft = 0.2540 m).
_RADIUS_FT = 0.8333333
#: Segment boundaries along the hull, ft: bow tip, midbody start, taper start,
#: stern-cap start, stern tip.
_BOW_END_FT = 3.333333
_MID_END_FT = 10.645833
_STERN_END_FT = 13.979167
#: Bow polynomial coefficients and the stern-taper/cap constants, as published.
_BOW_C1 = 1.126395101
_BOW_C2 = 0.442874707
_BOW_EXPONENT = 1.0 / 2.1
_TAPER_R1 = 0.1175
_TAPER_K0 = 10.0
_TAPER_K1 = 44.6244
_TAIL_C = 44.733333
_TAIL_K = 3.2
_TAIL_SCALE = 0.1175

#: International foot, exact by definition (0.3048 m).
FEET_TO_METRES = 0.3048
#: Nominal dimensions in metres, for callers who want the reference model.
SUBOFF_LENGTH_M = _LENGTH_FT * FEET_TO_METRES
SUBOFF_RADIUS_M = _RADIUS_FT * FEET_TO_METRES


def suboff_radius_fraction(xi: float) -> float:
    """Normalised radius ``r(xi) / R_max`` of the DARPA bare hull, ``xi`` in [0, 1].

    Four segments, in the axial coordinate ``x = xi * L`` measured in feet:

    1. **Bow** (``0 <= x < 3.3333``): a quintic polynomial in ``x`` raised to
       ``1/2.1``, which gives the sharp, finely-tapered entrance of the real hull
       (the radius tends to zero with an infinite slope -- see
       :func:`suboff_reference_geometry` for what that does to area quadrature);
    2. **Parallel midbody** (``3.3333 <= x <= 10.6458``): ``r = R_max``;
    3. **Stern taper** (``10.6458 < x <= 13.9792``): a sextic in
       ``ksi = (13.979167 - x) / 3.333333`` fitted to the published offsets,
       ``r1 = 0.1175`` at ``ksi = 0``;
    4. **Stern cap** (``13.9792 < x <= L``): an ellipse of semi-axis
       ``0.1175 R_max`` closing the hull at the stern tip.
    """
    value = float(xi)
    if not math.isfinite(value):
        raise FieldError(f"xi must be finite, got {xi!r}")
    if not 0.0 <= value <= 1.0:
        raise FieldError(f"xi must lie in [0, 1], got {xi!r}")
    x_ft = value * _LENGTH_FT
    if x_ft < _BOW_END_FT:
        tmp = 0.3 * x_ft - 1.0
        square = tmp * tmp
        fourth = square * square
        polynomial = (
            _BOW_C1 * x_ft * fourth
            + _BOW_C2 * x_ft * x_ft * (square * tmp)
            + 1.0
            - fourth * (1.2 * x_ft + 1.0)
        )
        return float(max(polynomial, 0.0) ** _BOW_EXPONENT)
    if x_ft <= _MID_END_FT:
        return 1.0
    if x_ft <= _STERN_END_FT:
        ksi = (_STERN_END_FT - x_ft) / 3.333333
        square = ksi * ksi
        cube = square * ksi
        fourth = cube * ksi
        fifth = fourth * ksi
        sixth = fifth * ksi
        r1 = _TAPER_R1
        polynomial = (
            r1 * r1
            + r1 * _TAPER_K0 * square
            + (20.0 - 20.0 * r1 * r1 - 4.0 * r1 * _TAPER_K0 - _TAPER_K1 / 3.0) * cube
            + (-45.0 + 45.0 * r1 * r1 + 6.0 * r1 * _TAPER_K0 + _TAPER_K1) * fourth
            + (36.0 - 36.0 * r1 * r1 - 4.0 * r1 * _TAPER_K0 - _TAPER_K1) * fifth
            + (-10.0 + 10.0 * r1 * r1 + r1 * _TAPER_K0 + _TAPER_K1 / 3.0) * sixth
        )
        return math.sqrt(max(polynomial, 0.0))
    tail = 1.0 - (_TAIL_K * x_ft - _TAIL_C) ** 2
    return math.sqrt(max(tail, 0.0)) * _TAIL_SCALE


def suboff_radius_fraction_array(xi: Any) -> Any:
    """NumPy mirror of :func:`suboff_radius_fraction` (same formulas, vectorised).

    The sampling path evaluates the profile 10^5--10^6 times per extraction, so it
    needs the vectorised form; ``tests/test_suboff.py`` pins the two together on a
    dense grid, so the fast path cannot drift from the documented one.
    """
    import numpy as np

    values = np.clip(np.asarray(xi, dtype=np.float64), 0.0, 1.0)
    x_ft = values * _LENGTH_FT
    out = np.zeros_like(x_ft)

    bow = x_ft < _BOW_END_FT
    if np.any(bow):
        x = x_ft[bow]
        tmp = 0.3 * x - 1.0
        square = tmp * tmp
        fourth = square * square
        polynomial = (
            _BOW_C1 * x * fourth + _BOW_C2 * x * x * (square * tmp) + 1.0 - fourth * (1.2 * x + 1.0)
        )
        out[bow] = np.maximum(polynomial, 0.0) ** _BOW_EXPONENT

    mid = (~bow) & (x_ft <= _MID_END_FT)
    out[mid] = 1.0

    taper = (x_ft > _MID_END_FT) & (x_ft <= _STERN_END_FT)
    if np.any(taper):
        ksi = (_STERN_END_FT - x_ft[taper]) / 3.333333
        square = ksi * ksi
        cube = square * ksi
        fourth = cube * ksi
        fifth = fourth * ksi
        sixth = fifth * ksi
        r1 = _TAPER_R1
        polynomial = (
            r1 * r1
            + r1 * _TAPER_K0 * square
            + (20.0 - 20.0 * r1 * r1 - 4.0 * r1 * _TAPER_K0 - _TAPER_K1 / 3.0) * cube
            + (-45.0 + 45.0 * r1 * r1 + 6.0 * r1 * _TAPER_K0 + _TAPER_K1) * fourth
            + (36.0 - 36.0 * r1 * r1 - 4.0 * r1 * _TAPER_K0 - _TAPER_K1) * fifth
            + (-10.0 + 10.0 * r1 * r1 + r1 * _TAPER_K0 + _TAPER_K1 / 3.0) * sixth
        )
        out[taper] = np.sqrt(np.maximum(polynomial, 0.0))

    tail = x_ft > _STERN_END_FT
    if np.any(tail):
        polynomial = 1.0 - (_TAIL_K * x_ft[tail] - _TAIL_C) ** 2
        out[tail] = np.sqrt(np.maximum(polynomial, 0.0)) * _TAIL_SCALE

    return np.clip(out, 0.0, 1.0)


def suboff_dimensions(length: float | None = None, radius: float | None = None) -> tuple[float, float]:
    """Resolve ``(length, radius)``, keeping the DARPA L/D when only one is given."""
    nominal_length = SUBOFF_LENGTH_M if length is None else float(length)
    if not math.isfinite(nominal_length) or nominal_length <= 0.0:
        raise FieldError(f"length must be finite and positive, got {length!r}")
    ratio = SUBOFF_RADIUS_M / SUBOFF_LENGTH_M
    resolved = nominal_length * ratio if radius is None else float(radius)
    if not math.isfinite(resolved) or resolved <= 0.0:
        raise FieldError(f"radius must be finite and positive, got {radius!r}")
    return nominal_length, resolved


@dataclass(frozen=True)
class SuboffHullField(GeometryField):
    """The DARPA SUBOFF bare hull (AFF-1) as a signed-distance field.

    ``length`` and ``radius`` are in the field's declared ``units`` (metres by
    default). The hull occupies ``0 <= x <= length`` on the centreline axis, with
    the bow tip at ``x = 0`` and the stern tip at ``x = length``.
    """

    length: float = SUBOFF_LENGTH_M
    radius: float | None = None
    units: str = DEFAULT_UNITS
    coordinate_frame: str = DEFAULT_FRAME

    kind: ClassVar[str] = "suboff_bare_hull"
    quantity: ClassVar[str] = SIGNED_DISTANCE

    def __post_init__(self) -> None:
        length, radius = suboff_dimensions(self.length, self.radius)
        object.__setattr__(self, "length", length)
        object.__setattr__(self, "radius", radius)

    # -- shape ------------------------------------------------------------- #\

    @property
    def max_radius(self) -> float:
        """Hull radius, resolved at construction time (``None`` means "derive")."""
        radius = self.radius
        assert radius is not None, "resolved in __post_init__"
        return radius

    def radius_at(self, x: float) -> float:
        """Hull radius at axial position ``x`` (0 outside ``[0, length]``)."""
        if not 0.0 <= x <= self.length:
            return 0.0
        return self.max_radius * suboff_radius_fraction(x / self.length)

    def _profile(self) -> tuple[Any, Any, float]:
        """Sample grid for the closest-point search, cached on the instance."""
        cached = getattr(self, "_profile_cache", None)
        if cached is None:
            import numpy as np

            samples = 129
            xs = np.linspace(0.0, self.length, samples)
            rs = suboff_radius_fraction_array(xs / self.length) * self.max_radius
            cached = (xs, rs, xs[1] - xs[0])
            object.__setattr__(self, "_profile_cache", cached)
        return cached

    # -- the protocol ------------------------------------------------------ #

    def signed_distance(self, points: Iterable[Sequence[float]]) -> list[float]:
        batch = _points(points)
        if not batch:
            return []
        evaluate = self.numpy_callable()
        if evaluate is not None:
            import numpy as np

            return [float(value) for value in evaluate(np.asarray(batch, dtype=np.float64))]
        return [self._distance_scalar(point) for point in batch]

    def _closest_scalar(self, point: Point) -> tuple[float, float, float]:
        """Nearest meridian point ``(x, r)`` and the distance to it, pure Python."""
        x, rho = point[0], math.hypot(point[1], point[2])
        xs, rs, step = self._profile()
        index = min(range(len(xs)), key=lambda i: (x - xs[i]) ** 2 + (rho - rs[i]) ** 2)
        low = max(float(xs[0]), float(xs[index]) - float(step))
        high = min(float(xs[-1]), float(xs[index]) + float(step))
        golden = (math.sqrt(5.0) - 1.0) / 2.0
        left = high - golden * (high - low)
        right = low + golden * (high - low)

        def gap(position: float) -> float:
            radius = self.max_radius * suboff_radius_fraction(position / self.length)
            return (x - position) ** 2 + (rho - radius) ** 2

        left_value = gap(left)
        right_value = gap(right)
        for _ in range(60):
            if left_value < right_value:
                high, right, right_value = right, left, left_value
                left = high - golden * (high - low)
                left_value = gap(left)
            else:
                low, left, left_value = left, right, right_value
                right = low + golden * (high - low)
                right_value = gap(right)
        best = left if left_value < right_value else right
        return best, self.radius_at(best), math.sqrt(min(left_value, right_value))

    def _distance_scalar(self, point: Point) -> float:
        _, _, distance = self._closest_scalar(point)
        return -distance if self.contains([point])[0] else distance

    def gradient(self, points: Iterable[Sequence[float]]) -> list[Point]:
        """Exact gradient (unit vector to the nearest surface point).

        The gradient of a signed distance points *away from* the surface, so it is
        ``(p - q)/|p - q|`` outside the solid and ``(q - p)/|q - p|`` inside it,
        where ``q`` is the nearest surface point. On the centreline inside the body
        every direction is equally "outward" for an axisymmetric solid; the +y
        direction is returned there so the value is deterministic.
        """
        batch = _points(points)
        if not batch:
            return []
        gradients: list[Point] = []
        inside = self.contains(batch)
        for point, nearest, is_inside in zip(batch, self._nearest(batch), inside, strict=True):
            x, rho = point[0], math.hypot(point[1], point[2])
            sign = -1.0 if is_inside else 1.0
            along = sign * (x - nearest[0])
            radial_gap = sign * (rho - nearest[1])
            if rho > 0.0:
                direction = (along, radial_gap * point[1] / rho, radial_gap * point[2] / rho)
            else:
                direction = (along, radial_gap, 0.0)
            length = math.sqrt(sum(component * component for component in direction))
            if length <= 1e-15:  # pragma: no cover - exactly on a smooth surface
                gradients.append((0.0, 0.0, 0.0))
            else:
                gradients.append(
                    (direction[0] / length, direction[1] / length, direction[2] / length)
                )
        return gradients

    def _nearest(self, batch: list[Point]) -> list[tuple[float, float]]:
        """Nearest point on the meridian curve, as ``(x, r)`` pairs."""
        evaluate = self.numpy_callable()
        if evaluate is None:  # pragma: no cover - numpy is installed for every test
            return [(x, radius) for x, radius, _ in (self._closest_scalar(p) for p in batch)]
        import numpy as np

        arrays = np.asarray(batch, dtype=np.float64)
        x = arrays[:, 0]
        rho = np.hypot(arrays[:, 1], arrays[:, 2])
        xs, rs, step = self._profile()
        coarse = (x[:, None] - xs[None, :]) ** 2 + (rho[:, None] - rs[None, :]) ** 2
        index = np.argmin(coarse, axis=1)
        low = np.maximum(xs[0], xs[index] - step)
        high = np.minimum(xs[-1], xs[index] + step)
        golden = (math.sqrt(5.0) - 1.0) / 2.0
        left = high - golden * (high - low)
        right = low + golden * (high - low)

        def gap(position: Any) -> Any:
            radius = suboff_radius_fraction_array(position / self.length) * self.max_radius
            return (x - position) ** 2 + (rho - radius) ** 2

        left_value = gap(left)
        right_value = gap(right)
        for _ in range(48):
            take_left = left_value < right_value
            high = np.where(take_left, right, high)
            low = np.where(take_left, low, left)
            left = high - golden * (high - low)
            right = low + golden * (high - low)
            left_value = gap(left)
            right_value = gap(right)
        best = np.where(left_value < right_value, left, right)
        profile = suboff_radius_fraction_array(best / self.length) * self.max_radius
        return [
            (float(position), float(radius_value))
            for position, radius_value in zip(best, profile, strict=True)
        ]

    def bounds(self) -> BoundingBox:
        radius = self.max_radius
        return BoundingBox((0.0, -radius, -radius), (self.length, radius, radius))

    def contains(self, points: Iterable[Sequence[float]]) -> list[bool]:
        """Exact containment in closed form, independent of the distance search."""
        out = []
        for point in _points(points):
            rho = math.hypot(point[1], point[2])
            out.append(bool(0.0 <= point[0] <= self.length and rho <= self.radius_at(point[0])))
        return out

    def state(self) -> dict[str, float]:
        return {"length": self.length, "radius": self.max_radius}

    def with_state(self, **values: float) -> SuboffHullField:
        unknown = set(values) - {"length", "radius"}
        if unknown:
            raise FieldError(f"{self.kind} has no control(s) {sorted(unknown)}")
        return SuboffHullField(
            length=float(values.get("length", self.length)),
            radius=float(values.get("radius", self.max_radius)),
            units=self.units,
            coordinate_frame=self.coordinate_frame,
        )

    def parameters(self) -> dict[str, Any]:
        return {
            "length": self.length,
            "radius": self.max_radius,
            "length_over_diameter": self.length / (2.0 * self.max_radius),
            "axial_origin": "bow tip at x = 0, stern tip at x = length",
            "profile": "DARPA SUBOFF bare hull, 4-segment offset table (feet internally)",
            "profile_source": (
                "Groves, Huang & Chang (1989), DTRC/SHD-1298-01, via "
                "TensorLBM src/tensorlbm/suboff_cad.py::suboff_radius_profile"
            ),
        }

    def provenance(self) -> dict[str, Any]:
        return {
            "producer": "TensorCAD",
            "representation": "analytic",
            "kind": self.kind,
            "geometry_source": "DARPA SUBOFF AFF-1 bare hull offset table (Groves et al. 1989)",
            "transcription": "TensorLBM suboff_cad.suboff_radius_profile, read 2026-10-24",
            "units": self.units,
            "coordinate_frame": self.coordinate_frame,
        }

    def numpy_callable(self) -> Callable[[Any], Any] | None:
        try:
            import numpy as np
        except ImportError:  # pragma: no cover - numpy is optional for the core
            return None
        xs, rs, step = self._profile()
        length = self.length
        radius_max = self.max_radius

        def evaluate(points: Any) -> Any:
            arrays = np.asarray(points, dtype=np.float64).reshape(-1, 3)
            x = arrays[:, 0]
            rho = np.hypot(arrays[:, 1], arrays[:, 2])
            # Coarse scan on the cached meridian, then golden-section refinement:
            # the closest point is unique for these profiles and the coarse grid
            # brackets its basin well inside one cell.
            coarse = (x[:, None] - xs[None, :]) ** 2 + (rho[:, None] - rs[None, :]) ** 2
            index = np.argmin(coarse, axis=1)
            low = np.maximum(xs[0], xs[index] - step)
            high = np.minimum(xs[-1], xs[index] + step)
            golden = 0.6180339887498949
            left = high - golden * (high - low)
            right = low + golden * (high - low)

            def gap(position: Any) -> Any:
                profile = suboff_radius_fraction_array(position / length) * radius_max
                return (x - position) ** 2 + (rho - profile) ** 2

            left_value = gap(left)
            right_value = gap(right)
            for _ in range(48):
                take_left = left_value < right_value
                high = np.where(take_left, right, high)
                low = np.where(take_left, low, left)
                left = high - golden * (high - low)
                right = low + golden * (high - low)
                left_value = gap(left)
                right_value = gap(right)
            distance = np.sqrt(np.minimum(left_value, right_value))
            inside = (x >= 0.0) & (x <= length) & (
                rho <= suboff_radius_fraction_array(x / length) * radius_max
            )
            return np.where(inside, -distance, distance)

        return evaluate


def suboff_reference_geometry(
    length: float = SUBOFF_LENGTH_M,
    radius: float | None = None,
) -> dict[str, Any]:
    """Analytic volume, wetted area and form coefficients of the bare hull.

    Computed by quadrature over the meridian, and reported with the resolution it
    was computed at plus a Richardson-style check against a refined one: the bow
    tip has an infinite slope (``r ~ xi**0.476``), so its contribution converges
    slowly and a single number would overstate the precision. See
    ``docs/suboff.md`` for the measured convergence.
    """
    import numpy as np

    resolved_length, resolved_radius = suboff_dimensions(length, radius)
    results: dict[str, Any] = {
        "length": resolved_length,
        "radius": resolved_radius,
        "length_over_diameter": resolved_length / (2.0 * resolved_radius),
    }
    for label, samples in (("coarse", 100_001), ("fine", 1_000_001)):
        xi = np.linspace(0.0, 1.0, samples)
        profile = suboff_radius_fraction_array(xi)
        slope = np.gradient(profile, xi)
        volume = math.pi * resolved_radius**2 * resolved_length * float(np.trapezoid(profile**2, xi))
        metric = np.sqrt(1.0 + ((resolved_radius / resolved_length) * slope) ** 2)
        area = (
            2.0
            * math.pi
            * resolved_radius
            * resolved_length
            * float(np.trapezoid(profile * metric, xi))
        )
        mid_area = math.pi * resolved_radius**2
        results[label] = {
            "samples": samples,
            "volume": volume,
            "wetted_area": area,
            "prismatic_coefficient": volume / (mid_area * resolved_length),
            "midship_area": mid_area,
        }
    results["volume"] = results["fine"]["volume"]
    results["wetted_area"] = results["fine"]["wetted_area"]
    results["prismatic_coefficient"] = results["fine"]["prismatic_coefficient"]
    results["volume_resolution_change"] = (
        results["fine"]["volume"] / results["coarse"]["volume"] - 1.0
    )
    results["wetted_area_resolution_change"] = (
        results["fine"]["wetted_area"] / results["coarse"]["wetted_area"] - 1.0
    )
    return results
