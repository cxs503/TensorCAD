"""DARPA SUBOFF appendages (AFF-8) as analytic geometry fields: sail and stern fins.

The bare hull (:mod:`tensorcad.suboff`) is the easy half of the DARPA model: one
offset table, four segments, an exact axisymmetric distance. The appendages are
the half that looks like it needs solid modelling -- a conning-tower sail whose
planform is three polynomial segments with an elliptical cap, and four swept
NACA-thickness stern fins. It does not: both are closed-form predicates in the
DARPA report, and this module implements those formulas directly, so the AFF-8
model is analytic all the way down and needs no boolean engine.

Geometry source
---------------
Groves, Huang & Chang (1989), *Geometric Characteristics of DARPA SUBOFF Models*,
DTRC/SHD-1298-01, as transcribed by TensorLBM's
``src/tensorlbm/suboff_cad.py`` (``_SAIL_*`` / ``_FIN_*`` / ``_NACA_COEFFS``,
``suboff_sail_contains_points`` and ``suboff_fins_contain_points``, read
2026-10-24). Every constant below is that module's, in the feet it is published
in; ``tests/test_suboff_appendages.py`` pins the two implementations point by
point.

The sail
--------
Three axial segments give the sail's half width ``a(x)`` (ft):

* entrance, ``_SAIL_X1_START < x < _SAIL_X1_END``::

      d = 3.0720 (x - X1_START)
      a = Zmax sqrt(2.094759 A + 0.2071781 B + C)
      A = 2 d (d - 1)^4,  B = d^2 (d - 1)^3 / 3,  C = 1 - (d - 1)^4 (4 d + 1)

* middle, ``X1_END < x <= X2_END``: ``a = Zmax`` (the parallel-sided cap);
* exit, ``X2_END < x <= X3_END``::

      e = (X3_END - x) / 0.6822917
      a = Zmax (2.238361 e f^4 + 3.106529 e^2 f^3 + 1 - f^4 (4 e + 1)),  f = e - 1

with ``Zmax = 0.109375 ft``. The cross-section is a rectangle of half width
``a`` and height ``YTMP = 1.507813 ft`` (measured up from the hull's
centreline, so its lower part is buried in the hull) topped by a half-ellipse
of semi-axes ``(a, a/2)``. Writing the cross-section as the single radial
coordinate

    rho(x, y, z) = |y|                        for z <= YTMP
                 = sqrt(y^2 + (2 (z - YTMP))^2)   for z >  YTMP

turns both pieces into ``rho <= a(x)``, which is exactly the form the report's
predicate has -- and the reason this module needs no booleans: the sail is
``{0 < z, X1_START < x < X3_END, rho < a(x)}``.

The fins
--------
Four fins in a cross, each a swept NACA 4-digit section: chord
``c(r) = K r + C`` at span radius ``r`` from the hull axis (``K < 0`` sweeps the
leading edge forward as the span grows), root at ``R_inner``, tip at
``R_outer``, and a half-thickness ``h(s)`` at chord position ``s`` from
``(0.2969, 0.126, 0.3516, 0.2852, 0.1045)``. In the ``(x, r)`` plane the
planform is a quadrilateral, because ``c`` is linear in ``r``; in three
dimensions it is the same predicate twice -- once with ``r = |y|, t = |z|``
(the two side fins) and once with ``r = |z|, t = |y|`` (the two vertical ones).

Distance semantics
------------------
The sail's cap and the fin's thickness surface have no closed-form nearest
point, so these fields declare ``quantity = "implicit_function"``: the zero set
and the sign are exact, and the magnitude is the first-order normalised gap
``g / |grad g|`` (exact on each face's flat part, first-order accurate near the
surface). ``tests/test_suboff_appendages.py`` measures how close: along the
analytic surface normal the gap is within 0.013 mm (sail) and 0.0002 mm (fins)
of the true offset at a 2 mm probe, and ``|grad|`` is 1.0 to 0.4 % at the median.
The one place the normalisation is genuinely weak is the sail's knife-edge
leading and trailing tips, where ``da/dx`` is infinite; the reported magnitude
there is small, so a sampler sees a slightly smeared tip. The sign is still
exact, which is what decides the surface's topology.

Units and frame
---------------
Same convention as :class:`tensorcad.suboff.SuboffHullField`: ``length`` and the
query coordinates are in the field's ``units`` (metres by default), x runs bow
(x = 0) to stern, y is transverse and z is up. The formulas are evaluated in
feet internally (``feet = coordinate * 14.291667 / length``), so the appendages
scale with the hull's length exactly as TensorLBM's lattice predicates do.
``scale`` reproduces TensorLBM's ``sail_scale`` / ``fin_scale``: the sail grows
by ``scale`` about its own axial centre / centreplane / deck plane, and the fins
by ``scale`` about their trailing edge / root radius, so a scaled appendage
stays attached to the hull.
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
    IMPLICIT_FUNCTION,
    FieldError,
    GeometryField,
    Point,
    _points,
)
from .suboff import SUBOFF_LENGTH_M, SuboffHullField

__all__ = [
    "SuboffAppendagesField",
    "SuboffFinsField",
    "SuboffFullField",
    "SuboffSailField",
    "suboff_fin_chord",
    "suboff_sail_half_width",
]

# --------------------------------------------------------------------------- #
# The AFF-8 constants, in the feet their source is published in
# --------------------------------------------------------------------------- #

#: Total hull length, ft -- the appendages' scale reference (14.291667 ft).
_LENGTH_FT = 14.291667
#: Maximum hull radius / sail deck plane, ft (the sail grows out of it).
_DECK_FT = 0.8333333

#: Sail maximum half width (transverse), ft.
_SAIL_HALF_WIDTH_FT = 0.109375
#: Sail height of the parallel-sided body above the hull axis, ft.
_SAIL_BODY_TOP_FT = 1.507813
#: Sail segment boundaries along the hull axis, ft (entrance / middle / exit).
_SAIL_X_START_FT = 3.032986
_SAIL_X1_END_FT = 3.358507
_SAIL_X2_END_FT = 3.559028
_SAIL_X_END_FT = 4.241319
#: Sail entrance rate and exit length, ft, as published.
_SAIL_ENTRANCE_RATE = 3.0720
_SAIL_EXIT_LENGTH_FT = 0.6822917
#: Sail similarity anchors: axial centre of the footprint, and the deck plane.
_SAIL_X_CENTER_FT = 0.5 * (_SAIL_X_START_FT + _SAIL_X_END_FT)

#: Fin trailing-edge axial position along the hull, ft.
_FIN_TRAILING_FT = 13.146284
#: Fin span (from the hull axis) and its root/tip radii, ft.
_FIN_R_INNER_FT = 0.075
_FIN_R_OUTER_FT = 0.825
#: Fin chord law ``c(r) = K r + C``, ft.
_FIN_SWEEP_K = -0.466308
_FIN_SWEEP_C = 0.88859
#: NACA 4-digit half-thickness coefficients, as published for this section.
_NACA_COEFFS = (0.2969, 0.126, 0.3516, 0.2852, 0.1045)

# --------------------------------------------------------------------------- #
# Sail: half width a(x) and da/dx, both in feet
# --------------------------------------------------------------------------- #


def _sail_half_width_ft(x_ft: float) -> float:
    """Sail half width ``a(x)`` in ft, 0 outside the three segments."""
    if x_ft <= _SAIL_X_START_FT or x_ft > _SAIL_X_END_FT:
        return 0.0
    if x_ft < _SAIL_X1_END_FT:
        d = _SAIL_ENTRANCE_RATE * (x_ft - _SAIL_X_START_FT)
        c = 1.0 - (d - 1.0) ** 4 * (4.0 * d + 1.0)
        b = (d * d * (d - 1.0) ** 3) / 3.0
        a = 2.0 * d * (d - 1.0) ** 4
        return _SAIL_HALF_WIDTH_FT * math.sqrt(max(2.094759 * a + 0.2071781 * b + c, 0.0))
    if x_ft <= _SAIL_X2_END_FT:
        return _SAIL_HALF_WIDTH_FT
    e = (_SAIL_X_END_FT - x_ft) / _SAIL_EXIT_LENGTH_FT
    f = e - 1.0
    return _SAIL_HALF_WIDTH_FT * max(
        2.238361 * e * f**4 + 3.106529 * e * e * f**3 + 1.0 - f**4 * (4.0 * e + 1.0),
        0.0,
    )


def _sail_half_width_slope_ft(x_ft: float) -> float:
    """``da/dx`` of :func:`_sail_half_width_ft`, in ft per ft.

    Written out rather than differenced: the entrance segment has an infinite
    slope at its start (``a ~ sqrt(x - X_START)``), which is a property of the
    published section, and a finite difference would only hide it.
    """
    if x_ft <= _SAIL_X_START_FT or x_ft > _SAIL_X_END_FT:
        return 0.0
    if x_ft < _SAIL_X1_END_FT:
        d = _SAIL_ENTRANCE_RATE * (x_ft - _SAIL_X_START_FT)
        c = 1.0 - (d - 1.0) ** 4 * (4.0 * d + 1.0)
        b = (d * d * (d - 1.0) ** 3) / 3.0
        a = 2.0 * d * (d - 1.0) ** 4
        dc = -(4.0 * (d - 1.0) ** 3 * (4.0 * d + 1.0) + 4.0 * (d - 1.0) ** 4)
        db = (2.0 * d * (d - 1.0) ** 3 + 3.0 * d * d * (d - 1.0) ** 2) / 3.0
        da = 2.0 * (d - 1.0) ** 4 + 8.0 * d * (d - 1.0) ** 3
        value = 2.094759 * a + 0.2071781 * b + c
        if value <= 0.0:  # pragma: no cover - the tip itself is excluded above
            return math.inf
        slope = 2.094759 * da + 0.2071781 * db + dc
        return _SAIL_HALF_WIDTH_FT * slope * _SAIL_ENTRANCE_RATE / (2.0 * math.sqrt(value))
    if x_ft <= _SAIL_X2_END_FT:
        return 0.0
    e = (_SAIL_X_END_FT - x_ft) / _SAIL_EXIT_LENGTH_FT
    f = e - 1.0
    slope = (
        2.238361 * (f**4 + 4.0 * e * f**3)
        + 3.106529 * (2.0 * e * f**3 + 3.0 * e * e * f * f)
        - (4.0 * f**3 * (4.0 * e + 1.0) + 4.0 * f**4)
    )
    return _SAIL_HALF_WIDTH_FT * slope * (-1.0 / _SAIL_EXIT_LENGTH_FT)


def _sail_half_width_array(x_ft: Any) -> Any:
    """NumPy mirror of :func:`_sail_half_width_ft`."""
    import numpy as np

    x = np.asarray(x_ft, dtype=np.float64)
    out = np.zeros_like(x)
    entrance = (x > _SAIL_X_START_FT) & (x < _SAIL_X1_END_FT)
    if np.any(entrance):
        d = _SAIL_ENTRANCE_RATE * (x[entrance] - _SAIL_X_START_FT)
        c = 1.0 - (d - 1.0) ** 4 * (4.0 * d + 1.0)
        b = (d * d * (d - 1.0) ** 3) / 3.0
        a = 2.0 * d * (d - 1.0) ** 4
        out[entrance] = _SAIL_HALF_WIDTH_FT * np.sqrt(np.maximum(2.094759 * a + 0.2071781 * b + c, 0.0))
    middle = (x >= _SAIL_X1_END_FT) & (x <= _SAIL_X2_END_FT)
    out[middle] = _SAIL_HALF_WIDTH_FT
    exit_ = (x > _SAIL_X2_END_FT) & (x <= _SAIL_X_END_FT)
    if np.any(exit_):
        e = (_SAIL_X_END_FT - x[exit_]) / _SAIL_EXIT_LENGTH_FT
        f = e - 1.0
        out[exit_] = _SAIL_HALF_WIDTH_FT * np.maximum(
            2.238361 * e * f**4 + 3.106529 * e * e * f**3 + 1.0 - f**4 * (4.0 * e + 1.0),
            0.0,
        )
    return out


def _sail_half_width_slope_array(x_ft: Any) -> Any:
    """NumPy mirror of :func:`_sail_half_width_slope_ft`."""
    import numpy as np

    x = np.asarray(x_ft, dtype=np.float64)
    out = np.zeros_like(x)
    entrance = (x > _SAIL_X_START_FT) & (x < _SAIL_X1_END_FT)
    if np.any(entrance):
        d = _SAIL_ENTRANCE_RATE * (x[entrance] - _SAIL_X_START_FT)
        c = 1.0 - (d - 1.0) ** 4 * (4.0 * d + 1.0)
        b = (d * d * (d - 1.0) ** 3) / 3.0
        a = 2.0 * d * (d - 1.0) ** 4
        dc = -(4.0 * (d - 1.0) ** 3 * (4.0 * d + 1.0) + 4.0 * (d - 1.0) ** 4)
        db = (2.0 * d * (d - 1.0) ** 3 + 3.0 * d * d * (d - 1.0) ** 2) / 3.0
        da = 2.0 * (d - 1.0) ** 4 + 8.0 * d * (d - 1.0) ** 3
        value = np.maximum(2.094759 * a + 0.2071781 * b + c, 1e-300)
        slope = 2.094759 * da + 0.2071781 * db + dc
        out[entrance] = _SAIL_HALF_WIDTH_FT * slope * _SAIL_ENTRANCE_RATE / (2.0 * np.sqrt(value))
    exit_ = (x > _SAIL_X2_END_FT) & (x <= _SAIL_X_END_FT)
    if np.any(exit_):
        e = (_SAIL_X_END_FT - x[exit_]) / _SAIL_EXIT_LENGTH_FT
        f = e - 1.0
        slope = (
            2.238361 * (f**4 + 4.0 * e * f**3)
            + 3.106529 * (2.0 * e * f**3 + 3.0 * e * e * f * f)
            - (4.0 * f**3 * (4.0 * e + 1.0) + 4.0 * f**4)
        )
        out[exit_] = _SAIL_HALF_WIDTH_FT * slope * (-1.0 / _SAIL_EXIT_LENGTH_FT)
    return out


# --------------------------------------------------------------------------- #
# The sail implicit function, in feet
# --------------------------------------------------------------------------- #


def _sail_implicit_ft(x_ft: float, y_ft: float, z_ft: float) -> float:
    """AFF-8 sail implicit value at a point in the DARPA feet frame."""
    width = _sail_half_width_ft(x_ft)
    slope = _sail_half_width_slope_ft(x_ft)
    rise = z_ft - _SAIL_BODY_TOP_FT
    cap = rise > 0.0
    rho = math.hypot(y_ft, 2.0 * rise) if cap else abs(y_ft)
    # |grad(rho - a)|: rho is |y| in the straight-sided body and
    # sqrt(y^2 + (2 rise)^2) under the cap, so the radial part of the gradient is
    # a unit y below the body top and (y/rho, 4 rise/rho) above it.
    radial = math.hypot(y_ft / rho, 4.0 * rise / rho) if cap and rho > 0.0 else 1.0
    gap = (rho - width) / math.hypot(radial, slope)
    return max(_SAIL_X_START_FT - x_ft, x_ft - _SAIL_X_END_FT, -z_ft, gap)


def _sail_implicit_array_ft(x_ft: Any, y_ft: Any, z_ft: Any) -> Any:
    """NumPy mirror of :func:`_sail_implicit_ft`."""
    import numpy as np

    x = np.asarray(x_ft, dtype=np.float64)
    y = np.asarray(y_ft, dtype=np.float64)
    z = np.asarray(z_ft, dtype=np.float64)
    width = _sail_half_width_array(x)
    slope = _sail_half_width_slope_array(x)
    rise = z - _SAIL_BODY_TOP_FT
    above = rise > 0.0
    rho = np.where(above, np.hypot(y, 2.0 * rise), np.abs(y))
    safe = np.maximum(rho, 1e-300)
    radial = np.where(
        above & (rho > 0.0),
        np.hypot(y / safe, 4.0 * rise / safe),
        1.0,
    )
    gap = (rho - width) / np.hypot(radial, slope)
    out = np.maximum(_SAIL_X_START_FT - x, x - _SAIL_X_END_FT)
    out = np.maximum(out, -z)
    return np.maximum(out, gap)


# --------------------------------------------------------------------------- #
# The fin implicit function, in feet
# --------------------------------------------------------------------------- #


def _naca_half_thickness_scalar(s: float) -> float:
    value = min(max(s, 0.0), 1.0)
    a, b, c, d, e = _NACA_COEFFS
    return max(a * math.sqrt(value) - b * value - c * value**2 + d * value**3 - e * value**4, 0.0)


def _naca_half_thickness_slope_scalar(s: float) -> float:
    value = min(max(s, 0.0), 1.0)
    a, b, c, d, e = _NACA_COEFFS
    if value <= 0.0:  # the leading edge: the section has an infinite slope there
        return math.inf
    return a * 0.5 / math.sqrt(value) - b - 2.0 * c * value + 3.0 * d * value**2 - 4.0 * e * value**3


def _naca_half_thickness_array(s: Any) -> Any:
    """NumPy mirror of :func:`_naca_half_thickness_scalar`."""
    import numpy as np

    value = np.clip(np.asarray(s, dtype=np.float64), 0.0, 1.0)
    a, b, c, d, e = _NACA_COEFFS
    return np.maximum(a * np.sqrt(value) - b * value - c * value**2 + d * value**3 - e * value**4, 0.0)


def _naca_half_thickness_slope_array(s: Any) -> Any:
    """NumPy mirror of :func:`_naca_half_thickness_slope_scalar`."""
    import numpy as np

    value = np.clip(np.asarray(s, dtype=np.float64), 0.0, 1.0)
    a, b, c, d, e = _NACA_COEFFS
    return (
        a * 0.5 / np.sqrt(np.maximum(value, 1e-300))
        - b
        - 2.0 * c * value
        + 3.0 * d * value**2
        - 4.0 * e * value**3
    )


def _fin_chord_coefficients(scale: float) -> tuple[float, float]:
    """Chord law ``(C_scale, r_outer)`` of a fin grown by ``scale``, ft.

    Scaling a fin about its trailing edge (chord) and root radius (span) keeps
    the chord linear in the span radius -- the scaled fin is the same trapezoid
    with a shifted intercept, which is why ``scale`` costs nothing here.
    """
    intercept = scale * _FIN_SWEEP_C + _FIN_SWEEP_K * _FIN_R_INNER_FT * (scale - 1.0)
    r_outer = _FIN_R_INNER_FT + scale * (_FIN_R_OUTER_FT - _FIN_R_INNER_FT)
    return intercept, r_outer


def _fin_planform_distance_ft(x_ft: float, r_ft: float, scale: float) -> float:
    """Exact 2D signed distance to the fin planform in the ``(x, r)`` plane, ft."""
    intercept, r_outer = _fin_chord_coefficients(scale)
    # Convex quadrilateral, counter-clockwise in (x, r): root leading edge ->
    # root trailing edge -> tip trailing edge -> tip leading edge.
    corners = (
        (_FIN_TRAILING_FT - (_FIN_SWEEP_K * _FIN_R_INNER_FT + intercept), _FIN_R_INNER_FT),
        (_FIN_TRAILING_FT, _FIN_R_INNER_FT),
        (_FIN_TRAILING_FT, r_outer),
        (_FIN_TRAILING_FT - (_FIN_SWEEP_K * r_outer + intercept), r_outer),
    )
    best = math.inf
    inside = True
    for index, (ax, ar) in enumerate(corners):
        bx, br = corners[(index + 1) % len(corners)]
        vx, vr = bx - ax, br - ar
        px, pr = x_ft - ax, r_ft - ar
        length_squared = vx * vx + vr * vr
        t = max(0.0, min(1.0, (px * vx + pr * vr) / length_squared))
        best = min(best, (px - t * vx) ** 2 + (pr - t * vr) ** 2)
        if vx * pr - vr * px < 0.0:
            inside = False
    distance = math.sqrt(best)
    return -distance if inside else distance


def _fin_planform_distance_array_ft(x_ft: Any, r_ft: Any, scale: float) -> Any:
    """NumPy mirror of :func:`_fin_planform_distance_ft`."""
    import numpy as np

    x = np.asarray(x_ft, dtype=np.float64)
    r = np.asarray(r_ft, dtype=np.float64)
    intercept, r_outer = _fin_chord_coefficients(scale)
    corners = (
        (_FIN_TRAILING_FT - (_FIN_SWEEP_K * _FIN_R_INNER_FT + intercept), _FIN_R_INNER_FT),
        (_FIN_TRAILING_FT, _FIN_R_INNER_FT),
        (_FIN_TRAILING_FT, r_outer),
        (_FIN_TRAILING_FT - (_FIN_SWEEP_K * r_outer + intercept), r_outer),
    )
    best = np.full(np.broadcast(x, r).shape, np.inf)
    signs = []
    for index, (ax, ar) in enumerate(corners):
        bx, br = corners[(index + 1) % len(corners)]
        vx, vr = bx - ax, br - ar
        px, pr = x - ax, r - ar
        t = np.clip((px * vx + pr * vr) / (vx * vx + vr * vr), 0.0, 1.0)
        best = np.minimum(best, (px - t * vx) ** 2 + (pr - t * vr) ** 2)
        signs.append(vx * pr - vr * px)
    stacked = np.stack([np.broadcast_to(sign, best.shape) for sign in signs])
    inside = np.all(stacked >= 0.0, axis=0) | np.all(stacked <= 0.0, axis=0)
    return np.where(inside, -np.sqrt(best), np.sqrt(best))


def _fin_plate_implicit_ft(x_ft: float, r_ft: float, t_ft: float, scale: float) -> float:
    """One fin plate: planform in ``(x, r)``, NACA thickness along ``t``, ft."""
    intercept, _ = _fin_chord_coefficients(scale)
    planform = _fin_planform_distance_ft(x_ft, r_ft, scale)
    chord = _FIN_SWEEP_K * r_ft + intercept
    if chord <= 0.0:  # pragma: no cover - the chord never closes for these constants
        return planform
    s = min(max((x_ft - _FIN_TRAILING_FT) / chord + 1.0, 0.0), 1.0)
    half = scale * _naca_half_thickness_scalar(s)
    slope = scale * _naca_half_thickness_slope_scalar(s)
    ds_dx = 1.0 / chord
    ds_dr = -(x_ft - _FIN_TRAILING_FT) * _FIN_SWEEP_K / (chord * chord)
    if not math.isfinite(slope):
        return planform
    denom = math.sqrt(1.0 + (slope * ds_dx) ** 2 + (slope * ds_dr) ** 2)
    return max(planform, (abs(t_ft) - half) / denom)


def _fin_plate_implicit_array_ft(x_ft: Any, r_ft: Any, t_ft: Any, scale: float) -> Any:
    """NumPy mirror of :func:`_fin_plate_implicit_ft`."""
    import numpy as np

    x = np.asarray(x_ft, dtype=np.float64)
    r = np.asarray(r_ft, dtype=np.float64)
    t = np.asarray(t_ft, dtype=np.float64)
    intercept, _ = _fin_chord_coefficients(scale)
    planform = _fin_planform_distance_array_ft(x, r, scale)
    chord = _FIN_SWEEP_K * r + intercept
    s = np.clip((x - _FIN_TRAILING_FT) / chord + 1.0, 0.0, 1.0)
    half = scale * _naca_half_thickness_array(s)
    slope = scale * _naca_half_thickness_slope_array(s)
    ds_dx = 1.0 / chord
    ds_dr = -(x - _FIN_TRAILING_FT) * _FIN_SWEEP_K / (chord * chord)
    denom = np.sqrt(1.0 + (slope * ds_dx) ** 2 + (slope * ds_dr) ** 2)
    return np.maximum(planform, (np.abs(t) - half) / denom)


def _fins_implicit_ft(x_ft: float, y_ft: float, z_ft: float, scale: float) -> float:
    """AFF-8 cruciform fins implicit value at a point in the DARPA feet frame."""
    ay, az = abs(y_ft), abs(z_ft)
    return min(
        _fin_plate_implicit_ft(x_ft, ay, az, scale),
        _fin_plate_implicit_ft(x_ft, az, ay, scale),
    )


def _fins_implicit_array_ft(x_ft: Any, y_ft: Any, z_ft: Any, scale: float) -> Any:
    """NumPy mirror of :func:`_fins_implicit_ft`."""
    import numpy as np

    x = np.asarray(x_ft, dtype=np.float64)
    ay = np.abs(np.asarray(y_ft, dtype=np.float64))
    az = np.abs(np.asarray(z_ft, dtype=np.float64))
    return np.minimum(
        _fin_plate_implicit_array_ft(x, ay, az, scale),
        _fin_plate_implicit_array_ft(x, az, ay, scale),
    )


# --------------------------------------------------------------------------- #
# The fields
# --------------------------------------------------------------------------- #


def _positive_scale(name: str, value: float) -> float:
    number = float(value)
    if not math.isfinite(number) or number <= 0.0:
        raise FieldError(f"{name} must be a finite positive number, got {value!r}")
    return number


def _positive_length(value: float) -> float:
    number = float(value)
    if not math.isfinite(number) or number <= 0.0:
        raise FieldError(f"length must be a finite positive number, got {value!r}")
    return number


@dataclass(frozen=True)
class SuboffSailField(GeometryField):
    """The AFF-8 sail as a field: three polynomial segments plus an elliptical cap.

    ``scale`` grows the sail about its own anchors (axial centre, centreplane and
    deck plane), exactly like TensorLBM's ``sail_scale``, so the sail keeps its
    foot on the hull.
    """

    length: float = SUBOFF_LENGTH_M
    scale: float = 1.0
    units: str = DEFAULT_UNITS
    coordinate_frame: str = DEFAULT_FRAME

    kind: ClassVar[str] = "suboff_sail"
    quantity: ClassVar[str] = IMPLICIT_FUNCTION

    def __post_init__(self) -> None:
        object.__setattr__(self, "length", _positive_length(self.length))
        object.__setattr__(self, "scale", _positive_scale("scale", self.scale))

    @property
    def feet_per_length_unit(self) -> float:
        """The DARPA feet per field unit: how the whole model scales with ``length``."""
        return _LENGTH_FT / self.length

    @property
    def height(self) -> float:
        """Height of the cap apex above the hull axis, in the field's units.

        This is the sail's topmost z: the body top (1.507813 ft) plus half the
        maximum half width (the cap's semi-axis), grown by ``scale`` about the
        deck plane.
        """
        feet = self.feet_per_length_unit
        top_ft = _SAIL_BODY_TOP_FT + 0.5 * _SAIL_HALF_WIDTH_FT
        return (_DECK_FT + self.scale * (top_ft - _DECK_FT)) / feet

    def signed_distance(self, points: Iterable[Sequence[float]]) -> list[float]:
        batch = _points(points)
        if not batch:
            return []
        evaluate = self.numpy_callable()
        if evaluate is not None:
            import numpy as np

            return [float(value) for value in evaluate(np.asarray(batch, dtype=np.float64))]
        return [self._value_scalar(point) for point in batch]

    def _value_scalar(self, point: Point) -> float:
        feet = self.feet_per_length_unit
        inverse = 1.0 / self.scale
        x_ft = point[0] * feet
        y_ft = point[1] * feet
        z_ft = point[2] * feet
        value = _sail_implicit_ft(
            _SAIL_X_CENTER_FT + (x_ft - _SAIL_X_CENTER_FT) * inverse,
            y_ft * inverse,
            _DECK_FT + (z_ft - _DECK_FT) * inverse,
        )
        return self.scale * value / feet

    def contains(self, points: Iterable[Sequence[float]]) -> list[bool]:
        """Exact containment by the published predicate (independent of the distance)."""
        feet = self.feet_per_length_unit
        inverse = 1.0 / self.scale
        out = []
        for point in _points(points):
            x_ft = _SAIL_X_CENTER_FT + (point[0] * feet - _SAIL_X_CENTER_FT) * inverse
            y_ft = point[1] * feet * inverse
            z_ft = _DECK_FT + (point[2] * feet - _DECK_FT) * inverse
            if not z_ft > 0.0 or not _SAIL_X_START_FT < x_ft < _SAIL_X_END_FT:
                out.append(False)
                continue
            width = _sail_half_width_ft(x_ft)
            if z_ft > _SAIL_BODY_TOP_FT:
                rho = math.hypot(y_ft, 2.0 * (z_ft - _SAIL_BODY_TOP_FT))
            else:
                rho = abs(y_ft)
            out.append(bool(rho < width))
        return out

    def gradient(self, points: Iterable[Sequence[float]]) -> list[Point]:
        """Central differences of the implicit function, at 1e-6 of the model length."""
        batch = _points(points)
        if not batch:
            return []
        return self.numerical_gradient(batch, step=1e-6 * self.length)

    def bounds(self) -> BoundingBox:
        feet = self.feet_per_length_unit
        scale = self.scale
        # The cap's apex is the sail's topmost point: body top plus half the
        # maximum half width, both grown by `scale` about the deck plane.
        top_ft = _SAIL_BODY_TOP_FT + 0.5 * _SAIL_HALF_WIDTH_FT
        return BoundingBox(
            (
                (_SAIL_X_CENTER_FT + scale * (_SAIL_X_START_FT - _SAIL_X_CENTER_FT)) / feet,
                -scale * _SAIL_HALF_WIDTH_FT / feet,
                (_DECK_FT + scale * (0.0 - _DECK_FT)) / feet,
            ),
            (
                (_SAIL_X_CENTER_FT + scale * (_SAIL_X_END_FT - _SAIL_X_CENTER_FT)) / feet,
                scale * _SAIL_HALF_WIDTH_FT / feet,
                (_DECK_FT + scale * (top_ft - _DECK_FT)) / feet,
            ),
        )

    def state(self) -> dict[str, float]:
        return {"length": self.length, "scale": self.scale}

    def with_state(self, **values: float) -> SuboffSailField:
        unknown = set(values) - {"length", "scale"}
        if unknown:
            raise FieldError(f"{self.kind} has no control(s) {sorted(unknown)}")
        return SuboffSailField(
            length=float(values.get("length", self.length)),
            scale=float(values.get("scale", self.scale)),
            units=self.units,
            coordinate_frame=self.coordinate_frame,
        )

    def parameters(self) -> dict[str, Any]:
        feet = self.feet_per_length_unit
        return {
            "length": self.length,
            "scale": self.scale,
            "axial_extent_m": [
                (_SAIL_X_CENTER_FT + self.scale * (_SAIL_X_START_FT - _SAIL_X_CENTER_FT)) / feet,
                (_SAIL_X_CENTER_FT + self.scale * (_SAIL_X_END_FT - _SAIL_X_CENTER_FT)) / feet,
            ],
            "half_width_m": _SAIL_HALF_WIDTH_FT / feet,
            "body_top_m": _SAIL_BODY_TOP_FT / feet,
            "segments": "entrance / parallel middle / exit polynomial, elliptical cap",
            "profile_source": (
                "Groves, Huang & Chang (1989), DTRC/SHD-1298-01, via "
                "TensorLBM src/tensorlbm/suboff_cad.py::suboff_sail_contains_points"
            ),
        }

    def provenance(self) -> dict[str, Any]:
        return {
            "producer": "TensorCAD",
            "representation": "analytic",
            "kind": self.kind,
            "geometry_source": "DARPA SUBOFF AFF-8 sail, 3-segment polynomial with elliptical cap",
            "transcription": "TensorLBM suboff_cad._SAIL_* / suboff_sail_contains_points, read 2026-10-24",
            "units": self.units,
            "coordinate_frame": self.coordinate_frame,
        }

    def numpy_callable(self) -> Callable[[Any], Any] | None:
        try:
            import numpy as np
        except ImportError:  # pragma: no cover - numpy is optional for the core
            return None
        feet = self.feet_per_length_unit
        inverse = 1.0 / self.scale
        scale = self.scale
        center = _SAIL_X_CENTER_FT
        deck = _DECK_FT

        def evaluate(points: Any) -> Any:
            arrays = np.asarray(points, dtype=np.float64).reshape(-1, 3)
            x_ft = arrays[:, 0] * feet
            y_ft = arrays[:, 1] * feet
            z_ft = arrays[:, 2] * feet
            return (
                scale
                * _sail_implicit_array_ft(
                    center + (x_ft - center) * inverse,
                    y_ft * inverse,
                    deck + (z_ft - deck) * inverse,
                )
                / feet
            )

        return evaluate


@dataclass(frozen=True)
class SuboffFinsField(GeometryField):
    """The AFF-8 cruciform stern fins: a swept NACA section on four blades.

    ``scale`` grows each fin about its trailing edge and root radius, exactly like
    TensorLBM's ``fin_scale``, so the fins ride the stern instead of detaching.
    """

    length: float = SUBOFF_LENGTH_M
    scale: float = 1.0
    units: str = DEFAULT_UNITS
    coordinate_frame: str = DEFAULT_FRAME

    kind: ClassVar[str] = "suboff_fins"
    quantity: ClassVar[str] = IMPLICIT_FUNCTION

    def __post_init__(self) -> None:
        object.__setattr__(self, "length", _positive_length(self.length))
        object.__setattr__(self, "scale", _positive_scale("scale", self.scale))

    @property
    def feet_per_length_unit(self) -> float:
        return _LENGTH_FT / self.length

    def signed_distance(self, points: Iterable[Sequence[float]]) -> list[float]:
        batch = _points(points)
        if not batch:
            return []
        evaluate = self.numpy_callable()
        if evaluate is not None:
            import numpy as np

            return [float(value) for value in evaluate(np.asarray(batch, dtype=np.float64))]
        return [self._value_scalar(point) for point in batch]

    def _value_scalar(self, point: Point) -> float:
        feet = self.feet_per_length_unit
        return (
            _fins_implicit_ft(
                point[0] * feet,
                point[1] * feet,
                point[2] * feet,
                self.scale,
            )
            / feet
        )

    def contains(self, points: Iterable[Sequence[float]]) -> list[bool]:
        """Exact containment by the published predicate (independent of the distance)."""
        feet = self.feet_per_length_unit
        intercept, r_outer = _fin_chord_coefficients(self.scale)
        out = []
        for point in _points(points):
            x_ft = point[0] * feet
            y_ft = abs(point[1] * feet)
            z_ft = abs(point[2] * feet)
            inside = False
            for span, thickness in ((y_ft, z_ft), (z_ft, y_ft)):
                if not _FIN_R_INNER_FT < span < r_outer:
                    continue
                chord = _FIN_SWEEP_K * span + intercept
                s = (x_ft - _FIN_TRAILING_FT) / chord + 1.0
                if not 0.0 < s < 1.0:
                    continue
                if thickness < self.scale * _naca_half_thickness_scalar(s):
                    inside = True
                    break
            out.append(inside)
        return out

    def gradient(self, points: Iterable[Sequence[float]]) -> list[Point]:
        """Central differences of the implicit function, at 1e-6 of the model length."""
        batch = _points(points)
        if not batch:
            return []
        return self.numerical_gradient(batch, step=1e-6 * self.length)

    def bounds(self) -> BoundingBox:
        feet = self.feet_per_length_unit
        intercept, r_outer = _fin_chord_coefficients(self.scale)
        # The longest chord is at the root (the section is swept back-to-front), so
        # the most forward point of the planform is the root leading edge.
        leading = _FIN_TRAILING_FT - (_FIN_SWEEP_K * _FIN_R_INNER_FT + intercept)
        return BoundingBox(
            (leading / feet, -r_outer / feet, -r_outer / feet),
            (_FIN_TRAILING_FT / feet, r_outer / feet, r_outer / feet),
        )

    def state(self) -> dict[str, float]:
        return {"length": self.length, "scale": self.scale}

    def with_state(self, **values: float) -> SuboffFinsField:
        unknown = set(values) - {"length", "scale"}
        if unknown:
            raise FieldError(f"{self.kind} has no control(s) {sorted(unknown)}")
        return SuboffFinsField(
            length=float(values.get("length", self.length)),
            scale=float(values.get("scale", self.scale)),
            units=self.units,
            coordinate_frame=self.coordinate_frame,
        )

    def parameters(self) -> dict[str, Any]:
        feet = self.feet_per_length_unit
        intercept, r_outer = _fin_chord_coefficients(self.scale)
        return {
            "length": self.length,
            "scale": self.scale,
            "root_radius_m": _FIN_R_INNER_FT / feet,
            "tip_radius_m": r_outer / feet,
            "span_m": (r_outer - _FIN_R_INNER_FT) / feet,
            "chord_root_m": (_FIN_SWEEP_K * _FIN_R_INNER_FT + intercept) / feet,
            "chord_tip_m": (_FIN_SWEEP_K * r_outer + intercept) / feet,
            "sweep_slope": _FIN_SWEEP_K,
            "profile": "NACA 4-digit half thickness (0.2969, 0.126, 0.3516, 0.2852, 0.1045)",
            "profile_source": (
                "Groves, Huang & Chang (1989), DTRC/SHD-1298-01, via "
                "TensorLBM src/tensorlbm/suboff_cad.py::suboff_fins_contain_points"
            ),
        }

    def provenance(self) -> dict[str, Any]:
        return {
            "producer": "TensorCAD",
            "representation": "analytic",
            "kind": self.kind,
            "geometry_source": "DARPA SUBOFF AFF-8 cruciform stern fins, swept NACA section",
            "transcription": "TensorLBM suboff_cad._FIN_* / suboff_fins_contain_points, read 2026-10-24",
            "units": self.units,
            "coordinate_frame": self.coordinate_frame,
        }

    def numpy_callable(self) -> Callable[[Any], Any] | None:
        try:
            import numpy as np
        except ImportError:  # pragma: no cover
            return None
        feet = self.feet_per_length_unit
        scale = self.scale

        def evaluate(points: Any) -> Any:
            arrays = np.asarray(points, dtype=np.float64).reshape(-1, 3)
            return (
                _fins_implicit_array_ft(
                    arrays[:, 0] * feet,
                    arrays[:, 1] * feet,
                    arrays[:, 2] * feet,
                    scale,
                )
                / feet
            )

        return evaluate


@dataclass(frozen=True)
class SuboffAppendagesField(GeometryField):
    """The AFF-8 appendage set: sail and cruciform fins, unioned."""

    length: float = SUBOFF_LENGTH_M
    sail_scale: float = 1.0
    fin_scale: float = 1.0
    units: str = DEFAULT_UNITS
    coordinate_frame: str = DEFAULT_FRAME

    kind: ClassVar[str] = "suboff_appendages"
    quantity: ClassVar[str] = IMPLICIT_FUNCTION

    def __post_init__(self) -> None:
        object.__setattr__(self, "length", _positive_length(self.length))
        object.__setattr__(self, "sail_scale", _positive_scale("sail_scale", self.sail_scale))
        object.__setattr__(self, "fin_scale", _positive_scale("fin_scale", self.fin_scale))

    def sail(self) -> SuboffSailField:
        """The sail on its own, at this field's length and scale."""
        return SuboffSailField(self.length, self.sail_scale, self.units, self.coordinate_frame)

    def fins(self) -> SuboffFinsField:
        """The fins on their own, at this field's length and scale."""
        return SuboffFinsField(self.length, self.fin_scale, self.units, self.coordinate_frame)

    def pieces(self) -> tuple[GeometryField, ...]:
        return (self.sail(), self.fins())

    def signed_distance(self, points: Iterable[Sequence[float]]) -> list[float]:
        batch = _points(points)
        if not batch:
            return []
        sail = self.sail().signed_distance(batch)
        fins = self.fins().signed_distance(batch)
        return [min(value, other) for value, other in zip(sail, fins, strict=True)]

    def contains(self, points: Iterable[Sequence[float]]) -> list[bool]:
        batch = _points(points)
        sail = self.sail().contains(batch)
        fins = self.fins().contains(batch)
        return [a or b for a, b in zip(sail, fins, strict=True)]

    def gradient(self, points: Iterable[Sequence[float]]) -> list[Point]:
        batch = _points(points)
        if not batch:
            return []
        sail = self.sail().signed_distance(batch)
        fins = self.fins().signed_distance(batch)
        sail_gradient = self.sail().gradient(batch)
        fins_gradient = self.fins().gradient(batch)
        return [
            sail_gradient[index] if sail[index] <= fins[index] else fins_gradient[index]
            for index in range(len(batch))
        ]

    def bounds(self) -> BoundingBox:
        sail = self.sail().bounds()
        fins = self.fins().bounds()
        assert sail is not None and fins is not None
        return BoundingBox(
            (
                min(sail.minimum[0], fins.minimum[0]),
                min(sail.minimum[1], fins.minimum[1]),
                min(sail.minimum[2], fins.minimum[2]),
            ),
            (
                max(sail.maximum[0], fins.maximum[0]),
                max(sail.maximum[1], fins.maximum[1]),
                max(sail.maximum[2], fins.maximum[2]),
            ),
        )

    def state(self) -> dict[str, float]:
        return {
            "length": self.length,
            "sail_scale": self.sail_scale,
            "fin_scale": self.fin_scale,
        }

    def with_state(self, **values: float) -> SuboffAppendagesField:
        unknown = set(values) - {"length", "sail_scale", "fin_scale"}
        if unknown:
            raise FieldError(f"{self.kind} has no control(s) {sorted(unknown)}")
        return SuboffAppendagesField(
            length=float(values.get("length", self.length)),
            sail_scale=float(values.get("sail_scale", self.sail_scale)),
            fin_scale=float(values.get("fin_scale", self.fin_scale)),
            units=self.units,
            coordinate_frame=self.coordinate_frame,
        )

    def parameters(self) -> dict[str, Any]:
        return {
            "length": self.length,
            "sail_scale": self.sail_scale,
            "fin_scale": self.fin_scale,
            "pieces": [self.sail().parameters(), self.fins().parameters()],
            "combination": "min of the sail and fin implicit functions (exact union sign)",
        }

    def provenance(self) -> dict[str, Any]:
        return {
            "producer": "TensorCAD",
            "representation": "analytic",
            "kind": self.kind,
            "geometry_source": "DARPA SUBOFF AFF-8 sail and cruciform stern fins",
            "transcription": "TensorLBM suboff_cad.suboff_appendages_contain_points, read 2026-10-24",
            "units": self.units,
            "coordinate_frame": self.coordinate_frame,
        }

    def numpy_callable(self) -> Callable[[Any], Any] | None:
        try:
            import numpy as np
        except ImportError:  # pragma: no cover
            return None
        sail = self.sail().numpy_callable()
        fins = self.fins().numpy_callable()
        if sail is None or fins is None:  # pragma: no cover
            return None

        def evaluate(points: Any) -> Any:
            arrays = np.asarray(points, dtype=np.float64).reshape(-1, 3)
            return np.minimum(sail(arrays), fins(arrays))

        return evaluate


@dataclass(frozen=True)
class SuboffFullField(GeometryField):
    """The complete DARPA SUBOFF: bare hull (AFF-1) with sail and fins (AFF-8).

    The three solids are unioned by ``min`` of their fields, which is the exact
    signed distance outside the solid and a bounded underestimate in the concave
    creases where the appendages meet the hull -- the standard caveat for a
    min-union, and the reason this field declares ``implicit_function`` rather
    than claiming an exact distance everywhere.

    The hull radius is the DARPA one for this length (``length * R / L``): the
    appendages are published as absolute coordinates in the DARPA frame, so a
    hull with a different radius would leave the sail and fins floating. Use
    :class:`tensorcad.suboff.SuboffHullField` for a rescaled radius, or
    ``length`` / ``scale`` here to resize the whole model.
    """

    length: float = SUBOFF_LENGTH_M
    sail_scale: float = 1.0
    fin_scale: float = 1.0
    units: str = DEFAULT_UNITS
    coordinate_frame: str = DEFAULT_FRAME

    kind: ClassVar[str] = "suboff_full"
    quantity: ClassVar[str] = IMPLICIT_FUNCTION

    def __post_init__(self) -> None:
        object.__setattr__(self, "length", _positive_length(self.length))
        object.__setattr__(self, "sail_scale", _positive_scale("sail_scale", self.sail_scale))
        object.__setattr__(self, "fin_scale", _positive_scale("fin_scale", self.fin_scale))

    def hull(self) -> SuboffHullField:
        """The bare hull at this field's length (DARPA radius)."""
        return SuboffHullField(self.length, None, self.units, self.coordinate_frame)

    def appendages(self) -> SuboffAppendagesField:
        """The sail and fins at this field's length and scales."""
        return SuboffAppendagesField(
            self.length, self.sail_scale, self.fin_scale, self.units, self.coordinate_frame
        )

    def pieces(self) -> tuple[GeometryField, ...]:
        return (self.hull(), self.appendages())

    def signed_distance(self, points: Iterable[Sequence[float]]) -> list[float]:
        batch = _points(points)
        if not batch:
            return []
        hull = self.hull().signed_distance(batch)
        appendages = self.appendages().signed_distance(batch)
        return [min(value, other) for value, other in zip(hull, appendages, strict=True)]

    def contains(self, points: Iterable[Sequence[float]]) -> list[bool]:
        batch = _points(points)
        hull = self.hull().contains(batch)
        appendages = self.appendages().contains(batch)
        return [a or b for a, b in zip(hull, appendages, strict=True)]

    def gradient(self, points: Iterable[Sequence[float]]) -> list[Point]:
        batch = _points(points)
        if not batch:
            return []
        hull = self.hull().signed_distance(batch)
        appendages = self.appendages().signed_distance(batch)
        hull_gradient = self.hull().gradient(batch)
        appendages_gradient = self.appendages().gradient(batch)
        return [
            hull_gradient[index] if hull[index] <= appendages[index] else appendages_gradient[index]
            for index in range(len(batch))
        ]

    def bounds(self) -> BoundingBox:
        hull = self.hull().bounds()
        appendages = self.appendages().bounds()
        assert hull is not None and appendages is not None
        return BoundingBox(
            (
                min(hull.minimum[0], appendages.minimum[0]),
                min(hull.minimum[1], appendages.minimum[1]),
                min(hull.minimum[2], appendages.minimum[2]),
            ),
            (
                max(hull.maximum[0], appendages.maximum[0]),
                max(hull.maximum[1], appendages.maximum[1]),
                max(hull.maximum[2], appendages.maximum[2]),
            ),
        )

    def state(self) -> dict[str, float]:
        return {
            "length": self.length,
            "sail_scale": self.sail_scale,
            "fin_scale": self.fin_scale,
        }

    def with_state(self, **values: float) -> SuboffFullField:
        unknown = set(values) - {"length", "sail_scale", "fin_scale"}
        if unknown:
            raise FieldError(f"{self.kind} has no control(s) {sorted(unknown)}")
        return SuboffFullField(
            length=float(values.get("length", self.length)),
            sail_scale=float(values.get("sail_scale", self.sail_scale)),
            fin_scale=float(values.get("fin_scale", self.fin_scale)),
            units=self.units,
            coordinate_frame=self.coordinate_frame,
        )

    def parameters(self) -> dict[str, Any]:
        return {
            "length": self.length,
            "radius": self.hull().max_radius,
            "sail_scale": self.sail_scale,
            "fin_scale": self.fin_scale,
            "pieces": [self.hull().parameters(), self.appendages().parameters()],
            "combination": "min of the hull and appendage fields (exact union sign)",
        }

    def provenance(self) -> dict[str, Any]:
        return {
            "producer": "TensorCAD",
            "representation": "analytic",
            "kind": self.kind,
            "geometry_source": "DARPA SUBOFF AFF-1 bare hull + AFF-8 sail and stern fins",
            "transcription": (
                "TensorLBM suboff_cad.suboff_radius_profile / suboff_sail_contains_points / "
                "suboff_fins_contain_points, read 2026-10-24"
            ),
            "units": self.units,
            "coordinate_frame": self.coordinate_frame,
        }

    def numpy_callable(self) -> Callable[[Any], Any] | None:
        try:
            import numpy as np
        except ImportError:  # pragma: no cover
            return None
        hull = self.hull().numpy_callable()
        appendages = self.appendages().numpy_callable()
        if hull is None or appendages is None:  # pragma: no cover
            return None

        def evaluate(points: Any) -> Any:
            arrays = np.asarray(points, dtype=np.float64).reshape(-1, 3)
            return np.minimum(hull(arrays), appendages(arrays))

        return evaluate


# --------------------------------------------------------------------------- #
# Reference helpers
# --------------------------------------------------------------------------- #


def suboff_sail_half_width(x: float, length: float = SUBOFF_LENGTH_M) -> float:
    """Sail half width at axial position ``x``, in the field's units (metres).

    ``x`` runs from the bow tip. This is the published ``a(x)`` of the AFF-8 sail
    (0 outside the three segments), in the DARPA frame scaled to ``length``.
    """
    feet = _LENGTH_FT / _positive_length(length)
    return _sail_half_width_ft(float(x) * feet) / feet


def suboff_fin_chord(radius: float, length: float = SUBOFF_LENGTH_M) -> float:
    """Fin chord at span ``radius`` from the hull axis, in the field's units.

    The chord law is ``c(r) = -0.466308 r + 0.88859`` in feet: the fins are swept
    so the leading edge moves forward as the span grows, and the tip chord is
    about 41 % shorter than the root chord.
    """
    feet = _LENGTH_FT / _positive_length(length)
    return (_FIN_SWEEP_K * float(radius) * feet + _FIN_SWEEP_C) / feet
