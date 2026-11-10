"""Parametric hull surface generation: NURBS cross-sections + longitudinal loft.

This module provides a complete pipeline from mathematical hull definition to
solver-ready surface mesh:

1. **Cross-section definition** -- from offset tables (station-half-breadth-height)
   or parametric functions (e.g. Wigley, SUBOFF radius profiles).
2. **NURBS curve interpolation** -- each cross-section is fitted with a NURBS curve;
   longitudinal curves (waterlines, buttock lines) are interpolated through the
   section control points.
3. **Surface construction** -- a tensor-product NURBS surface is lofted through the
   section curves, producing a smooth, watertight hull surface.
4. **Geometric transformations** -- scale, translate, rotate, mirror (port<->starboard)
   to support twin-screw / single-screw, bulbous bow / cruiser stern variants.
5. **Surface mesh generation** -- automatic triangulation of the parametric domain
   with density control (uniform or curvature-adaptive).
6. **Hydrostatic properties** -- displacement volume, centre of buoyancy, waterplane
   area, transverse metacentric radius, computed by divergence-theorem integration
   over the triangulated surface.

Coordinate frame
----------------
``x`` longitudinal (bow -> stern positive), ``y`` transverse (starboard positive),
``z`` vertical (up positive).  The waterline is at ``z = 0`` by default.

Dependencies
------------
Only the Python standard library.  The module re-uses :mod:`tensorcad.nurbs` for
B-spline evaluation and :mod:`tensorcad.geometry` / :mod:`tensorcad.mesh` for the
mesh contract and topology checks.
"""
from __future__ import annotations

import math
from collections.abc import Sequence
from dataclasses import dataclass, field
from typing import Any

from .geometry import Mesh, MeshBuilder, Vec3
from .nurbs import NURBSCurve, NURBSSurface, NURBSError, find_span, _basis_values
from . import mesh as meshlib

__all__ = [
    "CrossSection",
    "OffsetTable",
    "ParametricHull",
    "ParametricHullError",
    "SurfaceMesh",
    "HydrostaticResult",
    "Transform",
    "build_parametric_hull",
    "cross_section_from_offsets",
    "cross_section_from_function",
    "wigley_cross_sections",
    "suboff_cross_sections",
    "loft_nurbs_surface",
    "tessellate_surface",
    "apply_transform",
    "compute_hydrostatics",
    "mesh_quality_report",
    "mesh_density_for",
    "interpolate_control_points",
    "wigley_family_sections",
    "suboff_station_fractions",
    "wigley_analytic_properties",
    "depth_exponent_for_block_coefficient",
]


class ParametricHullError(ValueError):
    """Raised when hull definition or transformation is geometrically invalid."""


# ---------------------------------------------------------------------------
# Cross-section definition
# ---------------------------------------------------------------------------
@dataclass(frozen=True)
class CrossSection:
    """A single transverse section at a given longitudinal station.

    Attributes
    ----------
    x : float
        Longitudinal position of the section.
    points : tuple[tuple[float, float], ...]
        Ordered (y, z) pairs defining the half-breadth profile.  ``y >= 0``
        (starboard half); ``z`` runs from keel (lowest) to deck (highest).
    """
    x: float
    points: tuple[tuple[float, float], ...]

    @property
    def max_half_breadth(self) -> float:
        return max(p[0] for p in self.points) if self.points else 0.0

    @property
    def max_height(self) -> float:
        return max(p[1] for p in self.points) if self.points else 0.0

    @property
    def min_z(self) -> float:
        return min(p[1] for p in self.points) if self.points else 0.0

    def area(self) -> float:
        """Section area by the trapezoidal rule (half-section, starboard)."""
        if len(self.points) < 2:
            return 0.0
        sorted_pts = sorted(self.points, key=lambda p: p[1])
        area = 0.0
        for i in range(len(sorted_pts) - 1):
            y0, z0 = sorted_pts[i]
            y1, z1 = sorted_pts[i + 1]
            area += 0.5 * (y0 + y1) * (z1 - z0)
        return area


@dataclass(frozen=True)
class OffsetTable:
    """Station-based offset table.

    Each entry is a tuple ``(x_frac, [(y0,z0), (y1,z1), ...])`` where ``x_frac``
    is the fractional station position in [0, 1] (0 = AP, 1 = FP).
    """
    stations: tuple[tuple[float, tuple[tuple[float, float], ...]], ...]

    def __post_init__(self) -> None:
        for frac, pts in self.stations:
            if not (0.0 <= frac <= 1.0):
                raise ParametricHullError(f"station fraction {frac} not in [0, 1]")
            if len(pts) < 2:
                raise ParametricHullError(
                    f"station at x={frac} needs >= 2 points, got {len(pts)}"
                )

    @property
    def num_stations(self) -> int:
        return len(self.stations)

    def to_cross_sections(self, length: float) -> tuple[CrossSection, ...]:
        """Convert to CrossSection objects scaled to the given hull length."""
        sections = []
        for frac, pts in self.stations:
            x = frac * length
            sections.append(CrossSection(x=x, points=tuple(pts)))
        return tuple(sections)


def cross_section_from_offsets(
    x: float,
    half_breadths: Sequence[float],
    heights: Sequence[float],
) -> CrossSection:
    """Build a cross-section from paired half-breadth and height arrays."""
    if len(half_breadths) != len(heights):
        raise ParametricHullError("half_breadths and heights must have same length")
    if len(half_breadths) < 2:
        raise ParametricHullError("need at least 2 points per section")
    points = tuple(
        (float(y), float(z)) for y, z in zip(half_breadths, heights)
    )
    return CrossSection(x=float(x), points=points)


def cross_section_from_function(
    x: float,
    profile_func: Any,
    n_points: int = 20,
    z_min: float = 0.0,
    z_max: float = 1.0,
) -> CrossSection:
    """Build a cross-section by sampling a parametric profile function.

    Parameters
    ----------
    x : float
        Longitudinal position.
    profile_func : callable
        ``profile_func(z) -> y`` returning half-breadth at height z.
    n_points : int
        Number of sample points along z.
    z_min, z_max : float
        Vertical extent of the section.
    """
    if n_points < 2:
        raise ParametricHullError("n_points must be >= 2")
    points = []
    for i in range(n_points):
        z = z_min + (z_max - z_min) * i / (n_points - 1)
        y = profile_func(z)
        if y < -1e-15:
            raise ParametricHullError(f"negative half-breadth y={y} at z={z}")
        points.append((max(y, 0.0), float(z)))
    return CrossSection(x=float(x), points=tuple(points))


# ---------------------------------------------------------------------------
# Standard hull families
# ---------------------------------------------------------------------------
def wigley_family_sections(
    length: float = 2.0,
    beam: float = 0.2,
    draft: float = 0.125,
    *,
    fore_exponent: float = 2.0,
    depth_exponent: float = 2.0,
    block_coefficient: float | None = None,
    n_stations: int = 11,
    n_points: int = 10,
) -> tuple[CrossSection, ...]:
    """Cross-sections of the generalised Wigley family.

    The half-breadth is

    .. math::  y(x, z) = \\frac{B}{2}\\left[1 - \\left(\\frac{2x - L}{L}\\right)^p\\right]
               \\left[1 - \\left(\\frac{|z|}{T}\\right)^q\\right]

    with ``x`` measured from the AP (``0``) to the FP (``L``) and ``z`` from the keel
    (``-T``) to the waterline (``0``).  ``p = q = 2`` is the classical Wigley III
    parabolic hull.  ``p`` (``fore_exponent``) drives the plan-form fullness and ``q``
    (``depth_exponent``) the section fullness, so the pair maps one-to-one onto the
    block coefficient ``Cb = p q / ((p+1)(q+1))``; passing ``block_coefficient`` solves
    ``q`` for the requested ``Cb`` at the given ``p`` (see
    :func:`depth_exponent_for_block_coefficient`).

    All properties are available in closed form -- see
    :func:`wigley_analytic_properties`.
    """
    if length <= 0 or beam <= 0 or draft <= 0:
        raise ParametricHullError("length, beam and draft must be positive")
    if n_stations < 2 or n_points < 2:
        raise ParametricHullError("need at least 2 stations and 2 points per section")
    if fore_exponent <= 0:
        raise ParametricHullError("fore_exponent must be positive")
    q = depth_exponent
    if block_coefficient is not None:
        q = depth_exponent_for_block_coefficient(block_coefficient, fore_exponent)
    if q <= 0:
        raise ParametricHullError("depth_exponent must be positive")

    p = fore_exponent
    sections = []
    for i in range(n_stations):
        x = length * i / (n_stations - 1)

        def profile(z: float, _x: float = x) -> float:
            ratio = abs(2.0 * _x - length) / length
            plan = 1.0 - ratio ** p if ratio < 1.0 else 0.0
            depth = 1.0 - (abs(z) / draft) ** q
            return 0.5 * beam * max(plan, 0.0) * max(depth, 0.0)

        sections.append(
            cross_section_from_function(
                x=x,
                profile_func=profile,
                n_points=n_points,
                z_min=-draft,
                z_max=0.0,
            )
        )
    return tuple(sections)


def depth_exponent_for_block_coefficient(
    block_coefficient: float, fore_exponent: float = 2.0
) -> float:
    """Depth exponent ``q`` giving ``Cb = p q / ((p+1)(q+1))`` for a fixed ``p``.

    ``q = Cb (p + 1) / (p - Cb (p + 1))``, defined for ``p > Cb / (1 - Cb)``: below that
    plan-form exponent the requested fullness cannot be reached by deepening the
    sections alone.
    """
    p = fore_exponent
    if not 0.0 < block_coefficient < 1.0:
        raise ParametricHullError("block_coefficient must lie strictly inside (0, 1)")
    denominator = p - block_coefficient * (p + 1.0)
    if denominator <= 1e-12:
        raise ParametricHullError(
            f"fore_exponent={p} cannot reach Cb={block_coefficient}; "
            f"need p > Cb/(1-Cb) = {block_coefficient / (1.0 - block_coefficient):.4f}"
        )
    return block_coefficient * (p + 1.0) / denominator


def wigley_analytic_properties(
    length: float = 2.0,
    beam: float = 0.2,
    draft: float = 0.125,
    *,
    fore_exponent: float = 2.0,
    depth_exponent: float = 2.0,
) -> dict[str, float]:
    """Closed-form hydrostatics of the generalised Wigley hull.

    With ``f(x) = 1 - |2x/L - 1|^p`` and ``g(z) = 1 - (|z|/T)^q``::

        V   = L B T p q / ((p+1)(q+1))          Cb = p q / ((p+1)(q+1))
        Aw  = L B p / (p+1)                     Cw = p / (p+1)
        I_T = L B^3 / 12 * (1 - 3/(p+1) + 3/(2p+1) - 1/(3p+1))
        I_L = L^3 B * (1/12 - 1/(4(p+3)))
        VCB = -T (q+1)/q (1/2 - 1/(q+2))         LCB = LCF = L/2

    For ``p = q = 2`` these reduce to the textbook Wigley III values
    ``V = 4/9 L B T``, ``Aw = 2/3 L B``, ``I_T = 0.0381 L B^3`` and
    ``VCB = -0.375 T``; ``tests/test_parametric_hull_dag.py`` pins them against the
    published values.
    """
    p, q = fore_exponent, depth_exponent
    volume = length * beam * draft * p * q / ((p + 1.0) * (q + 1.0))
    waterplane = length * beam * p / (p + 1.0)
    inertia_t = (
        length
        * beam ** 3
        / 12.0
        * (1.0 - 3.0 / (p + 1.0) + 3.0 / (2.0 * p + 1.0) - 1.0 / (3.0 * p + 1.0))
    )
    inertia_l = beam * length ** 3 * (1.0 / 12.0 - 1.0 / (4.0 * (p + 3.0)))
    vcb = -draft * (q + 1.0) / q * (0.5 - 1.0 / (q + 2.0))
    return {
        "length": length,
        "beam": beam,
        "draft": draft,
        "fore_exponent": p,
        "depth_exponent": q,
        "displaced_volume": volume,
        "waterplane_area": waterplane,
        "block_coefficient": p * q / ((p + 1.0) * (q + 1.0)),
        "waterplane_coefficient": p / (p + 1.0),
        "midship_coefficient": q / (q + 1.0),
        "transverse_waterplane_inertia": inertia_t,
        "longitudinal_waterplane_inertia": inertia_l,
        "transverse_metacentric_radius": inertia_t / volume,
        "lcb": length / 2.0,
        "lcf": length / 2.0,
        "vcb": vcb,
    }


def wigley_cross_sections(
    length: float = 2.0,
    beam: float = 0.2,
    draft: float = 0.125,
    n_stations: int = 11,
    n_points: int = 10,
) -> tuple[CrossSection, ...]:
    """Cross-sections of the Wigley III parabolic mathematical hull.

    The hull is defined by::

        y(x, z) = (B/2) * (1 - (2(x - L/2)/L)^2) * (1 - (|z|/T)^2)

    with ``x`` from the AP (0) to the FP (L) and ``z`` from the keel (-T) to the
    waterline (0).  It is the ``p = q = 2`` member of :func:`wigley_family_sections`.

    Analytic properties (see :func:`wigley_analytic_properties`)::

        Volume        V  = (4/9) * L * B * T      =>  Cb = 4/9  ~ 0.4444
        Waterplane    Aw = (2/3) * L * B          =>  Cw = 2/3  ~ 0.6667
        Inertia       I_T = 0.0381 * L * B^3,  I_L = L^3 B / 30
        Waterplane area was previously documented here as ``(4/9) L B``, which is
        the *volume* coefficient, not the waterplane one; corrected 2026-11-10.
    """
    return wigley_family_sections(
        length=length,
        beam=beam,
        draft=draft,
        fore_exponent=2.0,
        depth_exponent=2.0,
        n_stations=n_stations,
        n_points=n_points,
    )


def suboff_station_fractions(n_stations: int = 41) -> tuple[float, ...]:
    """Station fractions that actually resolve the DARPA SUBOFF bare hull.

    Uniform stations cannot: the elliptical stern cap only occupies the last 2.2% of
    the length (``x/L`` from 0.978 to 1.0), so a 21-station uniform set puts the whole
    cap inside one interval and the lofted radius there is wrong by tens of percent
    (measured 33.9% at ``x/L = 0.975`` before this helper existed).  This set is
    cosine-clustered towards both ends, pins the three published segment junctions
    (bow end, parallel-midbody end, stern-taper end) and adds four stations inside the
    stern cap.
    """
    if n_stations < 5:
        raise ParametricHullError("need at least 5 stations for the SUBOFF hull")
    from .suboff import _BOW_END_FT, _LENGTH_FT, _MID_END_FT, _STERN_END_FT

    fractions = [
        0.5 * (1.0 - math.cos(math.pi * i / (n_stations - 1))) for i in range(n_stations)
    ]
    fractions.extend(
        [_BOW_END_FT / _LENGTH_FT, _MID_END_FT / _LENGTH_FT, _STERN_END_FT / _LENGTH_FT]
    )
    cap_start = _STERN_END_FT / _LENGTH_FT
    fractions.extend(cap_start + (1.0 - cap_start) * k / 5.0 for k in range(1, 5))
    return tuple(sorted({round(value, 12) for value in fractions if 0.0 <= value <= 1.0}))


def suboff_cross_sections(
    length_m: float | None = None,
    n_stations: int = 21,
    n_points: int = 20,
    *,
    station_fractions: Sequence[float] | None = None,
) -> tuple[CrossSection, ...]:
    """Generate cross-sections for the DARPA SUBOFF bare hull.

    The SUBOFF is axisymmetric, so each section is a circle of radius r(x) taken from
    the four-segment offset table of Groves, Huang & Chang (1989) -- see
    :mod:`tensorcad.suboff`, whose transcription is pinned to the published values in
    ``tests/test_suboff.py``.

    ``station_fractions`` overrides the (uniform by default) station positions; pass
    :func:`suboff_station_fractions` for a set that resolves the stern cap.
    """
    from .suboff import (
        SUBOFF_LENGTH_M,
        SUBOFF_RADIUS_M,
        suboff_radius_fraction,
    )

    L = length_m if length_m is not None else SUBOFF_LENGTH_M
    R = SUBOFF_RADIUS_M * (L / SUBOFF_LENGTH_M)
    if station_fractions is not None:
        fractions = [float(f) for f in station_fractions]
    else:
        fractions = [i / (n_stations - 1) for i in range(n_stations)]
    sections = []
    for xi in fractions:
        x = xi * L
        r = R * suboff_radius_fraction(xi)

        # Circular section: sample the right half (y >= 0)
        points = []
        for j in range(n_points):
            theta = -math.pi / 2 + math.pi * j / (n_points - 1)
            y = r * abs(math.cos(theta))
            z = r * math.sin(theta)
            points.append((y, z))
        sections.append(CrossSection(x=x, points=tuple(points)))
    return tuple(sections)


# ---------------------------------------------------------------------------
# NURBS surface lofting
# ---------------------------------------------------------------------------
def _clamped_knots(degree: int, n_ctrl: int) -> tuple[float, ...]:
    """Uniform clamped knot vector on [0, 1]."""
    if n_ctrl < degree + 1:
        raise ParametricHullError(
            f"{n_ctrl} control points insufficient for degree {degree}"
        )
    n_inner = n_ctrl - degree - 1
    knots = [0.0] * (degree + 1)
    for i in range(1, n_inner + 1):
        knots.append(i / (n_inner + 1))
    knots.extend([1.0] * (degree + 1))
    return tuple(knots)


def _resample_section(section: CrossSection, n_points: int) -> CrossSection:
    """Resample a cross-section to a uniform number of points by linear interpolation."""
    pts = sorted(section.points, key=lambda p: p[1])
    if len(pts) == n_points:
        return section
    chords = [0.0]
    for i in range(1, len(pts)):
        dy = pts[i][0] - pts[i - 1][0]
        dz = pts[i][1] - pts[i - 1][1]
        chords.append(chords[-1] + math.sqrt(dy * dy + dz * dz))
    total = chords[-1]
    if total < 1e-15:
        return CrossSection(x=section.x, points=tuple(pts[:n_points]))

    new_pts = []
    for i in range(n_points):
        target = total * i / (n_points - 1)
        seg = 0
        for k in range(1, len(chords)):
            if chords[k] >= target:
                seg = k - 1
                break
        else:
            seg = len(chords) - 2
        seg_len = chords[seg + 1] - chords[seg]
        if seg_len < 1e-15:
            t = 0.0
        else:
            t = (target - chords[seg]) / seg_len
        y = pts[seg][0] + t * (pts[seg + 1][0] - pts[seg][0])
        z = pts[seg][1] + t * (pts[seg + 1][1] - pts[seg][1])
        new_pts.append((y, z))
    return CrossSection(x=section.x, points=tuple(new_pts))


# ---------------------------------------------------------------------------
# B-spline interpolation (Piegl & Tiller, "The NURBS Book", eqs. 9.1-9.8)
# ---------------------------------------------------------------------------
def _point_distance(a: Sequence[float], b: Sequence[float]) -> float:
    return math.sqrt(sum((a[i] - b[i]) ** 2 for i in range(len(a))))


def _chord_length_parameters(points: Sequence[Sequence[float]]) -> tuple[float, ...]:
    """Normalised chord-length parametrisation of a data point sequence."""
    n = len(points)
    if n < 2:
        raise ParametricHullError("need at least 2 data points to parametrise")
    cumulative = [0.0]
    for i in range(1, n):
        cumulative.append(cumulative[-1] + _point_distance(points[i], points[i - 1]))
    total = cumulative[-1]
    if total <= 1e-30:
        # Coincident points: fall back to a uniform parametrisation.
        return tuple(i / (n - 1) for i in range(n))
    return tuple(value / total for value in cumulative)


def _averaged_knots(degree: int, params: Sequence[float]) -> tuple[float, ...]:
    """Averaging knot vector ``U_{j+p} = (u_{j-1} + ... + u_{j+p-2}) / p`` (Piegl 9.8)."""
    n_ctrl = len(params)
    if n_ctrl < degree + 1:
        raise ParametricHullError(
            f"{n_ctrl} data points insufficient for degree {degree} interpolation"
        )
    last = n_ctrl - 1
    knots = [0.0] * (degree + 1)
    for j in range(1, last - degree + 1):
        total = 0.0
        for i in range(1, degree + 1):
            total += params[j + i - 1]
        knots.append(total / degree)
    knots.extend([1.0] * (degree + 1))
    return tuple(knots)


def _basis_matrix(
    params: Sequence[float],
    degree: int,
    knots: Sequence[float],
) -> list[list[float]]:
    """Square collocation matrix ``N_i,p(u_k)`` for global interpolation."""
    n = len(params)
    matrix: list[list[float]] = []
    for u in params:
        span = find_span(knots, degree, u)
        values = _basis_values(span, u, degree, knots)
        row = [0.0] * n
        for offset, value in enumerate(values):
            row[span - degree + offset] = value
        matrix.append(row)
    return matrix


def _solve_dense(matrix: Sequence[Sequence[float]], rhs: Sequence[float]) -> list[float]:
    """Gaussian elimination with partial pivoting (standard library only)."""
    n = len(matrix)
    a = [list(matrix[i]) + [rhs[i]] for i in range(n)]
    for col in range(n):
        pivot = max(range(col, n), key=lambda r: abs(a[r][col]))
        if abs(a[pivot][col]) < 1e-13:
            raise ParametricHullError(
                "singular interpolation system: duplicate or degenerate data points"
            )
        if pivot != col:
            a[col], a[pivot] = a[pivot], a[col]
        pivot_value = a[col][col]
        for row in range(col + 1, n):
            factor = a[row][col] / pivot_value
            if factor == 0.0:
                continue
            for k in range(col, n + 1):
                a[row][k] -= factor * a[col][k]
    solution = [0.0] * n
    for row in range(n - 1, -1, -1):
        acc = a[row][n]
        for k in range(row + 1, n):
            acc -= a[row][k] * solution[k]
        solution[row] = acc / a[row][row]
    return solution


def interpolate_control_points(
    points: Sequence[Sequence[float]],
    degree: int = 3,
    *,
    params: Sequence[float] | None = None,
    knots: Sequence[float] | None = None,
) -> tuple[tuple[tuple[float, ...], ...], tuple[float, ...]]:
    """Control points of a B-spline curve that **interpolates** ``points``.

    Solves the global interpolation system of Piegl & Tiller (eqs. 9.1-9.8): with the
    chord-length parametrisation ``u_k``, the averaged knot vector and the collocation
    matrix ``N_i,p(u_k)``, the control points are ``P = N^-1 Q``.  Returns the control
    points and the knot vector so a caller can loft many curves on a shared
    parametrisation (skinning).

    Using the sampled points *directly* as control points (the historical behaviour of
    :func:`loft_nurbs_surface`) only reproduces the data at the clamped ends, which
    biases every interior ordinate towards the control polygon (for a Wigley section
    that is a ~3% displacement error -- see ``docs/parametric_hull`` verification).
    """
    n = len(points)
    if n < 2:
        raise ParametricHullError("need at least 2 data points to interpolate")
    actual_degree = max(1, min(degree, n - 1))
    if params is None:
        params = _chord_length_parameters(points)
    elif len(params) != n:
        raise ParametricHullError("parametrisation length does not match data points")
    if knots is None:
        knots = _averaged_knots(actual_degree, params)
    elif len(knots) != n + actual_degree + 1:
        raise ParametricHullError("knot vector length does not match data points")
    matrix = _basis_matrix(params, actual_degree, knots)
    dim = len(points[0])
    coords = [[float(p[c]) for p in points] for c in range(dim)]
    control = [[0.0] * dim for _ in range(n)]
    for c in range(dim):
        solution = _solve_dense(matrix, coords[c])
        for i in range(n):
            control[i][c] = solution[i]
    return tuple(tuple(row) for row in control), tuple(knots)


def _common_chord_parameters(columns: Sequence[Sequence[Sequence[float]]]) -> tuple[float, ...]:
    """Mean normalised chord-length parametrisation shared by parallel data columns.

    A tensor-product loft needs *one* parametrisation per direction; averaging the
    normalised chord lengths of every column keeps the knots identical for all
    interpolated curves, which is what makes the resulting surface pass through every
    data point instead of only through the boundary rows.
    """
    n = len(columns[0])
    accumulator = [0.0] * n
    counted = 0
    for column in columns:
        if len(column) != n:
            raise ParametricHullError("columns must have equal length")
        if any(_point_distance(column[i], column[0]) > 1e-14 for i in range(1, n)):
            params = _chord_length_parameters(column)
            for i in range(n):
                accumulator[i] += params[i]
            counted += 1
    if counted == 0:
        return tuple(i / (n - 1) for i in range(n))
    return tuple(value / counted for value in accumulator)


def loft_nurbs_surface(
    sections: Sequence[CrossSection],
    degree_u: int = 3,
    degree_v: int = 3,
    n_points_per_section: int | None = None,
    *,
    interpolate: bool = True,
) -> NURBSSurface:
    """Loft a NURBS surface through a sequence of cross-sections.

    With ``interpolate=True`` (default) the surface **interpolates** every resampled
    section point: each section is first interpolated in the transverse direction
    (``v``), then the resulting control nets are interpolated along the hull (``u``).
    The result passes through all section data, so ordinates, volume and waterplane
    area converge to the analytic geometry instead of to a smoothed control polygon.

    Parameters
    ----------
    sections : sequence of CrossSection
        Ordered from AP to FP.
    degree_u, degree_v : int
        Spline degrees, clamped to ``n - 1`` when fewer data are available.
    n_points_per_section : int, optional
        If given, resample each section to this many points before lofting.
    interpolate : bool
        ``False`` keeps the legacy behaviour of using the section points directly as
        control points (a smoothed approximant, not an interpolant).  Kept for
        regression comparison only.
    """
    if len(sections) < 2:
        raise ParametricHullError("need at least 2 sections to loft a surface")

    n_v = n_points_per_section if n_points_per_section is not None else len(sections[0].points)
    if n_v < 2:
        raise ParametricHullError("need at least 2 points per section")
    resampled = [_resample_section(sec, n_v) for sec in sections]
    n_u = len(resampled)

    degree_u_actual = min(degree_u, n_u - 1)
    degree_v_actual = min(degree_v, n_v - 1)

    if not interpolate:
        control_net = []
        for section in resampled:
            row = []
            for y, z in section.points:
                row.append((section.x, y, z))
            control_net.append(tuple(row))
        knots_u = _clamped_knots(degree_u_actual, n_u)
        knots_v = _clamped_knots(degree_v_actual, n_v)
        weights = [[1.0] * n_v for _ in range(n_u)]
        return NURBSSurface(
            degree=(degree_u_actual, degree_v_actual),
            knots=(tuple(knots_u), tuple(knots_v)),
            control_net=tuple(control_net),
            weights=weights,
        )

    # --- v direction: interpolate each section on a shared parametrisation -------
    section_points = [
        [(section.x, y, z) for (y, z) in section.points] for section in resampled
    ]
    # A *shared* uniform parametrisation across the profile: constant-v rows then sit
    # at the same relative position on every section, so the tensor product does not
    # mix sections evaluated at different heights.  Chord-length parameterisation here
    # (per section, then averaged) leaves a fixed ~4e-4 displacement floor that no
    # amount of refinement can remove -- measured, then removed by this choice.
    params_v = tuple(i / (n_v - 1) for i in range(n_v))
    knots_v = _averaged_knots(degree_v_actual, params_v)
    section_control = [
        interpolate_control_points(points, degree_v_actual, params=params_v, knots=knots_v)[0]
        for points in section_points
    ]

    # --- u direction: interpolate the control nets along the hull -----------------
    columns = [[section_control[i][j] for i in range(n_u)] for j in range(n_v)]
    params_u = _common_chord_parameters(columns)
    knots_u = _averaged_knots(degree_u_actual, params_u)
    column_control = [
        interpolate_control_points(column, degree_u_actual, params=params_u, knots=knots_u)[0]
        for column in columns
    ]

    control_net = tuple(
        tuple(column_control[j][i] for j in range(n_v)) for i in range(n_u)
    )
    weights = [[1.0] * n_v for _ in range(n_u)]
    return NURBSSurface(
        degree=(degree_u_actual, degree_v_actual),
        knots=(knots_u, knots_v),
        control_net=control_net,
        weights=weights,
    )


# ---------------------------------------------------------------------------
# Surface tessellation (mesh generation)
# ---------------------------------------------------------------------------
@dataclass(frozen=True)
class SurfaceMesh:
    """A triangulated surface mesh with quality metadata.

    ``hull_face_count`` separates the hull surface from closure faces appended
    afterwards (the waterplane lid of :func:`tessellate_surface`); hydrostatics and
    mesh-quality claims are always evaluated on the hull faces only.
    """
    mesh: Mesh
    n_u: int
    n_v: int
    hull_face_count: int = -1
    waterplane_z: float | None = None

    @property
    def n_vertices(self) -> int:
        return len(self.mesh.vertices)

    @property
    def n_faces(self) -> int:
        return len(self.mesh.faces)

    @property
    def hull_mesh(self) -> Mesh:
        """The hull surface without closure lids."""
        if self.hull_face_count < 0 or self.hull_face_count >= len(self.mesh.faces):
            return self.mesh
        return Mesh(
            vertices=self.mesh.vertices,
            faces=self.mesh.faces[: self.hull_face_count],
            units=self.mesh.units,
            coordinate_frame=self.mesh.coordinate_frame,
        )


def _x_parameter_table(
    surface: NURBSSurface, samples: int = 2001
) -> tuple[list[float], list[float]]:
    """Tabulate ``x(u)`` at ``v = 0.5`` for parametrisation inversion."""
    us = [i / (samples - 1) for i in range(samples)]
    xs = [surface.evaluate(u, 0.5)[0] for u in us]
    return us, xs


def _invert_x(us: Sequence[float], xs: Sequence[float], target: float) -> float:
    """Parameter ``u`` whose surface ``x`` equals ``target`` (monotone table)."""
    if target <= xs[0]:
        return us[0]
    if target >= xs[-1]:
        return us[-1]
    lo, hi = 0, len(xs) - 1
    while hi - lo > 1:
        mid = (lo + hi) // 2
        if xs[mid] <= target:
            lo = mid
        else:
            hi = mid
    span = xs[hi] - xs[lo]
    if span <= 1e-300:
        return us[lo]
    return us[lo] + (target - xs[lo]) / span * (us[hi] - us[lo])


def _monotone_x_table(
    surface: NURBSSurface, samples: int = 2001
) -> tuple[list[float], list[float]] | None:
    """``(u, x)`` table with ``x`` forced non-decreasing, or ``None`` if not usable.

    ``None`` means the surface genuinely folds back in ``x`` (an overhang or a
    re-entrant patch) and physical-abscissa sampling is not meaningful.  Sub-micron
    cubic overshoot at a station jump is repaired instead of rejected: the SUBOFF loft
    dips 4e-5 m below ``x = 0`` at the bow tip, which is a numerical artifact of the
    interpolant, not a shape feature.
    """
    us, xs = _x_parameter_table(surface, samples)
    if xs[-1] < xs[0]:
        us = us[::-1]
        xs = xs[::-1]
    span = abs(xs[-1] - xs[0])
    if span <= 0.0:
        return None
    allowed_backstep = 1e-4 * span
    for i in range(1, len(xs)):
        if xs[i] < xs[i - 1]:
            if xs[i - 1] - xs[i] > allowed_backstep:
                return None
            xs[i] = xs[i - 1]
    return us, xs


def _cosine_x_parameters(surface: NURBSSurface, n_u: int) -> list[float]:
    """``n_u`` parameters clustered towards both ends (cosine station spacing).

    Required for hulls that *close* at the ends with a square-root radius profile (the
    DARPA SUBOFF bow has ``r ~ x**0.476``, i.e. an infinite slope at the tip): the
    physical spacing there has to shrink like ``sqrt(x)`` for the parametric cells to
    stay well shaped.  Uniform spacing leaves stretched cells at such a tip; cosine
    clustering reduces the worst aspect ratio by 3-4x on the same body.
    """
    if n_u < 3:
        raise ParametricHullError("cosine spacing needs at least 3 samples")
    table = _monotone_x_table(surface)
    if table is None:
        return [i / (n_u - 1) for i in range(n_u)]
    us, xs = table
    x_lo, x_hi = xs[0], xs[-1]
    targets = [
        x_lo + (x_hi - x_lo) * 0.5 * (1.0 - math.cos(math.pi * i / (n_u - 1)))
        for i in range(n_u)
    ]
    return [_invert_x(us, xs, target) for target in targets]


def _uniform_x_parameters(surface: NURBSSurface, n_u: int) -> list[float]:
    """``n_u`` parameters whose surface abscissae ``x`` are uniformly spaced.

    A clamped B-spline's parameter is *not* proportional to the physical abscissa: on
    a lofted hull the first knot span can cover several times the mean station
    spacing, which produced 14:1 sliver triangles at the stem in the pre-fix
    tessellator.  Sampling in physical ``x`` removes that distortion at the source.
    """
    if n_u < 2:
        raise ParametricHullError("need at least 2 samples along the hull")
    table = _monotone_x_table(surface)
    if table is None:
        return [i / (n_u - 1) for i in range(n_u)]
    us, xs = table
    x_lo, x_hi = xs[0], xs[-1]
    return [
        _invert_x(us, xs, x_lo + (x_hi - x_lo) * i / (n_u - 1)) for i in range(n_u)
    ]


def tessellate_surface(
    surface: NURBSSurface,
    n_u: int = 40,
    n_v: int = 14,
    *,
    mirror_y: bool = True,
    close_waterplane: bool = False,
    u_spacing: str = "uniform_x",
    weld: bool = True,
    weld_tolerance: float | None = None,
) -> SurfaceMesh:
    """Tessellate a NURBS hull patch into a triangulated surface mesh.

    Parameters
    ----------
    surface : NURBSSurface
        The half-breadth (starboard) hull patch.
    n_u, n_v : int
        Mesh density control: evaluation points along the hull (``u``) and across the
        section profile (``v``).  Triangles are well shaped when the physical spacing
        ``dL/(n_u-1)`` stays within ~2-3x the profile spacing ``dT/(n_v-1)``;
        :func:`mesh_density_for` derives a consistent pair.
    mirror_y : bool
        Mirror about the xz plane to build the full (port + starboard) body.  Vertices
        on the centre plane are shared, so the keel seam is watertight.
    close_waterplane : bool
        Append a planar lid at the top edge (``v = 1``) so the body is a closed
        manifold.  The lid is *not* counted in ``hull_face_count`` and is excluded
        from hull mesh-quality statistics.
    u_spacing : {"uniform_x", "uniform"}
        ``uniform_x`` samples the surface at uniformly spaced physical ``x``
        (default); ``uniform`` samples uniformly in the spline parameter.
    weld : bool
        Merge coincident vertices (centre-line seam, degenerate stem columns).
    """
    if n_u < 2 or n_v < 2:
        raise ParametricHullError("need at least 2x2 evaluation points")
    if u_spacing == "uniform_x":
        u_samples = _uniform_x_parameters(surface, n_u)
    elif u_spacing == "cosine":
        u_samples = _cosine_x_parameters(surface, n_u)
    elif u_spacing == "uniform":
        u_samples = [i / (n_u - 1) for i in range(n_u)]
    else:
        raise ParametricHullError(f"unknown u_spacing {u_spacing!r}")
    v_samples = [j / (n_v - 1) for j in range(n_v)]

    builder = MeshBuilder()
    grid: list[list[Vec3]] = [
        [surface.evaluate(u, v) for v in v_samples] for u in u_samples
    ]
    extent = max(max(abs(c) for pt in row for c in pt) for row in grid) or 1.0
    tol = weld_tolerance if weld_tolerance is not None else 1e-12 * extent

    starboard = [[builder.add_vertex(pt) for pt in row] for row in grid]
    on_seam = [[abs(pt[1]) <= tol for pt in row] for row in grid] if mirror_y else None
    if mirror_y:
        port: list[list[int]] = []
        for row_index, row in enumerate(grid):
            port_row: list[int] = []
            for column_index, pt in enumerate(row):
                if abs(pt[1]) <= tol:
                    # On the centre plane: share the starboard vertex (keel seam).
                    port_row.append(starboard[row_index][column_index])
                else:
                    port_row.append(builder.add_vertex((pt[0], -pt[1], pt[2])))
            port.append(port_row)

    def _quad_triangles(i: int, j: int) -> tuple[tuple[int, int, int], ...]:
        """Split a quad along its *shorter* diagonal.

        Choosing the shorter diagonal improves the worst triangle and, on a
        mirror-symmetric hull, keeps the triangulation itself symmetric: a fixed
        corner-to-corner split reverses its direction under ``x -> L - x``, which
        biases the longitudinal centre of buoyancy by ~1 mm on a 10 m hull (measured)
        even though the vertex set is perfectly symmetric.
        """
        p00, p10 = grid[i][j], grid[i + 1][j]
        p01, p11 = grid[i][j + 1], grid[i + 1][j + 1]
        d02 = math.dist(p00, p11)
        d13 = math.dist(p10, p01)
        if d02 <= d13:
            return ((0, 1, 3), (0, 3, 2))
        return ((0, 1, 2), (1, 3, 2))

    for i in range(n_u - 1):
        for j in range(n_v - 1):
            corners_starboard = (
                starboard[i][j],
                starboard[i + 1][j],
                starboard[i][j + 1],
                starboard[i + 1][j + 1],
            )
            corners_port = None
            if mirror_y:
                corners_port = (
                    port[i][j],
                    port[i + 1][j],
                    port[i][j + 1],
                    port[i + 1][j + 1],
                )
            for a, b, c in _quad_triangles(i, j):
                s_tri = (corners_starboard[a], corners_starboard[b], corners_starboard[c])
                if len(set(s_tri)) < 3:
                    continue
                if mirror_y:
                    assert on_seam is not None
                    if on_seam[i + (1 if a in (1, 3) else 0)][j + (1 if a in (2, 3) else 0)] and \
                       on_seam[i + (1 if b in (1, 3) else 0)][j + (1 if b in (2, 3) else 0)] and \
                       on_seam[i + (1 if c in (1, 3) else 0)][j + (1 if c in (2, 3) else 0)]:
                        # Meridian-plane flap: it would be mirrored onto itself.
                        continue
                    builder.add_triangle(*s_tri)
                    p_tri = (corners_port[a], corners_port[c], corners_port[b])
                    if len(set(p_tri)) == 3:
                        builder.add_triangle(*p_tri)
                else:
                    builder.add_triangle(*s_tri)

    mesh = builder.build()
    if weld:
        mesh = meshlib.weld(mesh, tolerance=tol)
    mesh = meshlib.remove_degenerate_faces(mesh, tolerance=0.0)
    hull_face_count = len(mesh.faces)

    waterplane_z: float | None = None
    if close_waterplane:
        top = [row[-1] for row in grid]
        zs = [pt[2] for pt in top]
        waterplane_z = sum(zs) / len(zs)
        if max(abs(z - waterplane_z) for z in zs) > max(tol, 1e-9 * extent):
            raise ParametricHullError(
                "close_waterplane requires a planar top edge (constant z at v=1)"
            )
        mesh = _append_waterplane_lid(mesh, top, mirror_y)
        # Stitch the lid to the hull top row (exact coordinates) and drop any lid
        # sliver at the stem; face order is preserved so ``hull_face_count`` stays
        # valid.
        mesh = meshlib.weld(mesh, tolerance=tol)
        mesh = meshlib.remove_degenerate_faces(mesh, tolerance=0.0)

    return SurfaceMesh(
        mesh=mesh,
        n_u=n_u,
        n_v=n_v,
        hull_face_count=hull_face_count,
        waterplane_z=waterplane_z,
    )


def _append_waterplane_lid(
    mesh: Mesh,
    top_edge: Sequence[Vec3],
    mirror_y: bool,
) -> Mesh:
    """Close the top edge of a hull patch with a waterplane lid (outward +z).

    The lid reuses the exact top-edge coordinates so a following :func:`mesh.weld`
    stitches it to the hull without slivers at the seam.
    """
    vertices = list(mesh.vertices)
    faces = list(mesh.faces)

    def add(pt: Vec3) -> int:
        vertices.append(pt)
        return len(vertices) - 1

    seams = [add((pt[0], 0.0, pt[2])) for pt in top_edge]
    starboard = [add((pt[0], abs(pt[1]), pt[2])) for pt in top_edge]
    port = (
        [add((pt[0], -abs(pt[1]), pt[2])) for pt in top_edge] if mirror_y else []
    )
    for i in range(len(top_edge) - 1):
        a, b = seams[i], seams[i + 1]
        s0, s1 = starboard[i], starboard[i + 1]
        faces.append((a, s1, s0))
        faces.append((a, b, s1))
        if mirror_y:
            p0, p1 = port[i], port[i + 1]
            faces.append((a, p0, p1))
            faces.append((a, p1, b))
    return Mesh(
        vertices=tuple(vertices),
        faces=tuple(faces),
        units=mesh.units,
        coordinate_frame=mesh.coordinate_frame,
    )


def mesh_density_for(
    length: float,
    profile_extent: float,
    *,
    target_aspect: float = 2.5,
    n_v: int | None = None,
) -> tuple[int, int]:
    """Consistent ``(n_u, n_v)`` for near-isotropic triangles on a hull patch.

    A structured grid on a hull with a sharp stem has two competing limits: the stem
    triangles degenerate as ``dz / |dx|``, while the midbody triangles stretch as
    ``dx / dz``.  Both stay acceptable (scaled Jacobian > 0.35, aspect ratio < 3)
    when the physical spacings satisfy ``dx = target_aspect * dz``; solving that for
    the requested ``n_v`` (default: 12 profile rows) gives the pair below.  ``length``
    is the hull length and ``profile_extent`` the profile height (draft for a canoe
    body).
    """
    if length <= 0 or profile_extent <= 0:
        raise ParametricHullError("length and profile_extent must be positive")
    rows = n_v if n_v is not None else 12
    if rows < 2:
        raise ParametricHullError("need at least 2 profile rows")
    dz = profile_extent / (rows - 1)
    columns = int(round(length / (target_aspect * dz))) + 1
    return max(columns, 3), rows


# ---------------------------------------------------------------------------
# Geometric transformations
# ---------------------------------------------------------------------------
@dataclass(frozen=True)
class Transform:
    """Rigid-body + scale transformation for hull surfaces.

    All transformations are applied in order: mirror -> scale -> rotate -> translate.
    """
    scale: tuple[float, float, float] = (1.0, 1.0, 1.0)
    rotation_deg: tuple[float, float, float] = (0.0, 0.0, 0.0)
    translation: tuple[float, float, float] = (0.0, 0.0, 0.0)
    mirror_x: bool = False
    mirror_y: bool = False
    mirror_z: bool = False

    def apply_to_point(self, pt: Vec3) -> Vec3:
        x, y, z = pt
        if self.mirror_x:
            x = -x
        if self.mirror_y:
            y = -y
        if self.mirror_z:
            z = -z
        x *= self.scale[0]
        y *= self.scale[1]
        z *= self.scale[2]
        rx = math.radians(self.rotation_deg[0])
        ry = math.radians(self.rotation_deg[1])
        rz = math.radians(self.rotation_deg[2])
        y1 = y * math.cos(rx) - z * math.sin(rx)
        z1 = y * math.sin(rx) + z * math.cos(rx)
        x2 = x * math.cos(ry) + z1 * math.sin(ry)
        z2 = -x * math.sin(ry) + z1 * math.cos(ry)
        x3 = x2 * math.cos(rz) - y1 * math.sin(rz)
        y3 = x2 * math.sin(rz) + y1 * math.cos(rz)
        return (x3 + self.translation[0], y3 + self.translation[1], z2 + self.translation[2])


def apply_transform(mesh: Mesh, transform: Transform) -> Mesh:
    """Apply a geometric transformation to all vertices of a mesh."""
    new_vertices = tuple(transform.apply_to_point(v) for v in mesh.vertices)
    return Mesh(
        vertices=new_vertices,
        faces=mesh.faces,
        units=mesh.units,
        coordinate_frame=mesh.coordinate_frame,
    )


# ---------------------------------------------------------------------------
# Hydrostatic property computation
# ---------------------------------------------------------------------------
@dataclass(frozen=True)
class HydrostaticResult:
    """Hydrostatic properties computed from a triangulated hull surface."""
    displaced_volume: float
    centre_of_buoyancy: tuple[float, float, float]
    waterplane_area: float
    wetted_surface_area: float
    block_coefficient: float
    transverse_metacentric_radius: float
    longitudinal_waterplane_inertia: float
    transverse_waterplane_inertia: float
    lcb: float
    vcb: float
    lcf: float = 0.0
    displaced_mass: float = 0.0

    def to_dict(self) -> dict[str, float]:
        return {
            "displaced_volume_m3": self.displaced_volume,
            "centre_of_buoyancy": self.centre_of_buoyancy,
            "waterplane_area_m2": self.waterplane_area,
            "wetted_surface_area_m2": self.wetted_surface_area,
            "block_coefficient": self.block_coefficient,
            "transverse_metacentric_radius_m": self.transverse_metacentric_radius,
            "longitudinal_waterplane_inertia_m4": self.longitudinal_waterplane_inertia,
            "transverse_waterplane_inertia_m4": self.transverse_waterplane_inertia,
            "lcb_m": self.lcb,
            "vcb_m": self.vcb,
            "lcf_m": self.lcf,
            "displaced_mass_kg": self.displaced_mass,
        }


def _clip_below_plane(
    vertices: Sequence[Vec3], waterline_z: float
) -> list[Vec3]:
    """Sutherland-Hodgman clip of a convex polygon to ``z <= waterline_z``."""
    out: list[Vec3] = []
    n = len(vertices)
    for i in range(n):
        a = vertices[i]
        b = vertices[(i + 1) % n]
        za = a[2] - waterline_z
        zb = b[2] - waterline_z
        a_in = za <= 0.0
        b_in = zb <= 0.0
        if a_in:
            out.append(a)
        if a_in != b_in:
            t = za / (za - zb)
            out.append(tuple(a[k] + t * (b[k] - a[k]) for k in range(3)))  # type: ignore[misc]
    return out


def _triangle_integrals(vertices: Sequence[Vec3]) -> dict[str, float]:
    """Exact planar integrals of a triangle (all integrands at most quadratic).

    For a triangle with area ``A``::

        integral f dA          = A / 3 * sum(f_i)                  (f linear)
        integral f * g dA      = A / 12 * (sum(f_i) sum(g_i) + sum(f_i g_i))

    which is the standard trapezoid/quadrature rule for linear-``x``-linear products
    and is exact, so the divergence-theorem volume below carries no quadrature error.
    """
    v0, v1, v2 = vertices
    ax, ay, az = v1[0] - v0[0], v1[1] - v0[1], v1[2] - v0[2]
    bx, by, bz = v2[0] - v0[0], v2[1] - v0[1], v2[2] - v0[2]
    cx = ay * bz - az * by
    cy = az * bx - ax * bz
    cz = ax * by - ay * bx
    two_area = math.sqrt(cx * cx + cy * cy + cz * cz)
    area = 0.5 * two_area
    if two_area <= 0.0:
        return {
            "area": 0.0, "nx": 0.0, "ny": 0.0, "nz": 0.0,
            "x": 0.0, "y": 0.0, "z": 0.0,
            "xz": 0.0, "yz": 0.0, "xx": 0.0, "yy": 0.0, "zz": 0.0,
        }
    sx = v0[0] + v1[0] + v2[0]
    sy = v0[1] + v1[1] + v2[1]
    sz = v0[2] + v1[2] + v2[2]
    sxx = v0[0] ** 2 + v1[0] ** 2 + v2[0] ** 2
    syy = v0[1] ** 2 + v1[1] ** 2 + v2[1] ** 2
    szz = v0[2] ** 2 + v1[2] ** 2 + v2[2] ** 2
    sxz = v0[0] * v0[2] + v1[0] * v1[2] + v2[0] * v2[2]
    syz = v0[1] * v0[2] + v1[1] * v1[2] + v2[1] * v2[2]
    return {
        "area": area,
        "nx": cx / two_area,
        "ny": cy / two_area,
        "nz": cz / two_area,
        "x": area * sx / 3.0,
        "y": area * sy / 3.0,
        "z": area * sz / 3.0,
        "xz": area * (sx * sz + sxz) / 12.0,
        "yz": area * (sy * sz + syz) / 12.0,
        "xx": area * (sx * sx + sxx) / 12.0,
        "yy": area * (sy * sy + syy) / 12.0,
        "zz": area * (sz * sz + szz) / 12.0,
    }


def _wetted_integrals(
    mesh: Mesh, waterline_z: float, n_hull_faces: int | None = None
) -> dict[str, float]:
    """Surface integrals of the wetted hull, clipped exactly at the waterline.

    Each face is clipped to ``z <= z_wl`` and contributes its exact area, centroid,
    second moments and normal.  Clipping (instead of dropping faces by centroid) makes
    the volume and waterplane moments exact for any waterline position, not only when
    the waterline happens to fall on a mesh row.
    """
    faces = mesh.faces if n_hull_faces is None else mesh.faces[:n_hull_faces]
    acc = {
        "area": 0.0, "nz_a": 0.0, "z": 0.0, "xz": 0.0, "yz": 0.0,
        "zz": 0.0, "xx": 0.0, "yy": 0.0, "x": 0.0, "y": 0.0,
    }
    for face in faces:
        polygon = _clip_below_plane([mesh.vertices[i] for i in face], waterline_z)
        if len(polygon) < 3:
            continue
        anchor = polygon[0]
        for k in range(1, len(polygon) - 1):
            tri = (anchor, polygon[k], polygon[k + 1])
            terms = _triangle_integrals(tri)
            if terms["area"] == 0.0:
                continue
            nz = terms["nz"]
            acc["area"] += terms["area"]
            acc["nz_a"] += nz * terms["area"]
            acc["z"] += nz * terms["z"]
            acc["xz"] += nz * terms["xz"]
            acc["yz"] += nz * terms["yz"]
            acc["zz"] += nz * terms["zz"]
            acc["xx"] += nz * terms["xx"]
            acc["yy"] += nz * terms["yy"]
            acc["x"] += nz * terms["x"]
            acc["y"] += nz * terms["y"]
    return acc


def compute_hydrostatics(
    mesh: Mesh,
    waterline_z: float = 0.0,
    water_density: float = 1025.0,
    length: float | None = None,
    beam: float | None = None,
    draft: float | None = None,
    n_hull_faces: int | None = None,
) -> HydrostaticResult:
    """Hydrostatic properties of the hull volume below the waterline.

    The wetted surface is clipped at ``waterline_z`` and integrated with the divergence
    theorem.  The unknown lid closing the body (the waterplane) is never meshed: for
    any integrand independent of ``z`` the closed-surface integral vanishes, hence

    * ``Aw   = -integral_hull n_z dA``
    * ``I_T  = -integral_hull y^2 n_z dA``   (transverse waterplane inertia)
    * ``I_L  = -integral_hull x^2 n_z dA``   (longitudinal waterplane inertia)
    * ``V    =  integral_hull z n_z dA + z_wl * Aw``

    This is exact for the triangulated geometry and replaces the previous contour
    tracing of the waterline (which lost the faces that straddle the plane, giving a
    -2.1% waterplane error on Wigley III).

    ``n_hull_faces`` restricts the integration to the first ``n`` faces, which is how
    a closure lid appended by :func:`tessellate_surface` is excluded.
    """
    if not mesh.vertices:
        return HydrostaticResult(
            displaced_volume=0.0,
            centre_of_buoyancy=(0.0, 0.0, 0.0),
            waterplane_area=0.0,
            wetted_surface_area=0.0,
            block_coefficient=0.0,
            transverse_metacentric_radius=0.0,
            longitudinal_waterplane_inertia=0.0,
            transverse_waterplane_inertia=0.0,
            lcb=0.0,
            vcb=0.0,
        )
    xs = [v[0] for v in mesh.vertices]
    ys = [v[1] for v in mesh.vertices]
    zs = [v[2] for v in mesh.vertices]
    if length is None:
        length = max(xs) - min(xs)
    if beam is None:
        beam = max(ys) - min(ys)
    if draft is None:
        draft = waterline_z - min(zs)

    acc = _wetted_integrals(mesh, waterline_z, n_hull_faces)

    waterplane_area = -acc["nz_a"]
    volume = acc["z"] + waterline_z * waterplane_area
    moment_x = acc["xz"] + waterline_z * (-acc["x"])
    moment_y = acc["yz"] + waterline_z * (-acc["y"])
    moment_z = 0.5 * acc["zz"] + 0.5 * waterline_z * waterline_z * waterplane_area
    transverse_inertia = -acc["yy"]
    longitudinal_inertia = -acc["xx"]
    lid_x = -acc["x"]

    # Face winding may be inward; a negative signed volume reveals it.  Every quantity
    # above is linear in the normal, so one sign flip fixes them all.
    sign = -1.0 if volume < 0.0 else 1.0
    volume *= sign
    waterplane_area *= sign
    moment_x *= sign
    moment_y *= sign
    moment_z *= sign
    transverse_inertia *= sign
    longitudinal_inertia *= sign
    lid_x *= sign

    if volume > 1e-15:
        centre_of_buoyancy = (
            moment_x / volume,
            moment_y / volume,
            moment_z / volume,
        )
    else:
        centre_of_buoyancy = (0.0, 0.0, 0.0)

    # The surface integral above carries x about the AP origin; naval architecture
    # quotes the longitudinal waterplane inertia about the transverse axis through the
    # centre of flotation, so transfer with the parallel-axis theorem.
    lcf = lid_x / waterplane_area if abs(waterplane_area) > 1e-15 else 0.0
    longitudinal_inertia -= lcf * lcf * waterplane_area

    lbt = length * beam * draft
    block_coefficient = volume / lbt if lbt > 1e-15 else 0.0
    bm_transverse = transverse_inertia / volume if volume > 1e-15 else 0.0

    return HydrostaticResult(
        displaced_volume=volume,
        centre_of_buoyancy=centre_of_buoyancy,
        waterplane_area=waterplane_area,
        wetted_surface_area=acc["area"],
        block_coefficient=block_coefficient,
        transverse_metacentric_radius=bm_transverse,
        longitudinal_waterplane_inertia=longitudinal_inertia,
        transverse_waterplane_inertia=transverse_inertia,
        lcb=centre_of_buoyancy[0],
        vcb=centre_of_buoyancy[2],
        lcf=lcf,
        displaced_mass=volume * water_density,
    )


def _compute_waterplane_area(mesh: Mesh, waterline_z: float) -> float:
    """Waterplane area (exactly, from the wetted-surface projected area)."""
    return abs(compute_hydrostatics(mesh, waterline_z=waterline_z).waterplane_area)


def _compute_waterplane_inertia(
    mesh: Mesh, waterline_z: float
) -> tuple[float, float]:
    """``(transverse, longitudinal)`` waterplane inertias about the centroid-free axes."""
    result = compute_hydrostatics(mesh, waterline_z=waterline_z)
    return (
        result.transverse_waterplane_inertia,
        result.longitudinal_waterplane_inertia,
    )


def mesh_quality_report(
    mesh: Mesh,
    *,
    n_faces: int | None = None,
) -> dict[str, Any]:
    """Triangle quality statistics for a surface mesh.

    ``min_scaled_jacobian`` is ``2A / (e_mid * e_max)`` per triangle (0 for a
    degenerate sliver, 0.87 for an equilateral triangle) and ``max_aspect_ratio`` is
    ``e_max / e_min``.  ``n_faces`` restricts the statistics to the first ``n`` faces,
    so a closure lid can be excluded from hull-quality claims (its shape is dictated by
    the body plan, not by the solver).
    """
    faces = mesh.faces if n_faces is None else mesh.faces[:n_faces]
    sub = Mesh(
        vertices=mesh.vertices,
        faces=tuple(faces),
        units=mesh.units,
        coordinate_frame=mesh.coordinate_frame,
    )
    is_watertight = meshlib.is_watertight(sub)
    is_manifold = meshlib.is_edge_manifold(sub)

    min_jacobian = 1.0
    max_aspect_ratio = 0.0
    aspect_sum = 0.0
    min_area = float("inf")
    n_degenerate = 0
    n_evaluated = 0

    for face in sub.faces:
        v0 = sub.vertices[face[0]]
        v1 = sub.vertices[face[1]]
        v2 = sub.vertices[face[2]]

        e0 = math.dist(v0, v1)
        e1 = math.dist(v1, v2)
        e2 = math.dist(v2, v0)
        edges = sorted((e0, e1, e2))
        if edges[0] < 1e-15:
            n_degenerate += 1
            min_jacobian = 0.0
            continue

        ux, uy, uz = v1[0] - v0[0], v1[1] - v0[1], v1[2] - v0[2]
        vx, vy, vz = v2[0] - v0[0], v2[1] - v0[1], v2[2] - v0[2]
        nx = uy * vz - uz * vy
        ny = uz * vx - ux * vz
        nz = ux * vy - uy * vx
        two_area = math.sqrt(nx * nx + ny * ny + nz * nz)
        if two_area <= 0.0:
            n_degenerate += 1
            min_jacobian = 0.0
            continue

        aspect = edges[2] / edges[0]
        max_aspect_ratio = max(max_aspect_ratio, aspect)
        aspect_sum += aspect
        min_area = min(min_area, 0.5 * two_area)
        denominator = edges[1] * edges[2]
        if denominator > 1e-15:
            min_jacobian = min(min_jacobian, two_area / denominator)
        n_evaluated += 1

    return {
        "watertight": is_watertight,
        "edge_manifold": is_manifold,
        "n_vertices": len(mesh.vertices),
        "n_faces": len(faces),
        "n_faces_evaluated": n_evaluated,
        "min_scaled_jacobian": min_jacobian if n_evaluated else 0.0,
        "max_aspect_ratio": max_aspect_ratio,
        "mean_aspect_ratio": aspect_sum / n_evaluated if n_evaluated else 0.0,
        "min_triangle_area": min_area if n_evaluated else 0.0,
        "n_degenerate_faces": n_degenerate,
    }


# ---------------------------------------------------------------------------
# High-level builder
# ---------------------------------------------------------------------------
@dataclass
class ParametricHull:
    """Complete parametric hull: sections -> surface -> mesh -> hydrostatics."""
    sections: tuple[CrossSection, ...]
    surface: NURBSSurface
    mesh_data: SurfaceMesh
    hydrostatics: HydrostaticResult
    quality: dict[str, Any]
    parameters: dict[str, float]

    def to_dict(self) -> dict[str, Any]:
        return {
            "parameters": self.parameters,
            "n_sections": len(self.sections),
            "surface_degree": list(self.surface.degree),
            "mesh": {
                "n_vertices": self.mesh_data.n_vertices,
                "n_faces": self.mesh_data.n_faces,
            },
            "hydrostatics": self.hydrostatics.to_dict(),
            "quality": self.quality,
        }


def build_parametric_hull(
    sections: Sequence[CrossSection],
    *,
    length: float | None = None,
    beam: float | None = None,
    draft: float | None = None,
    degree_u: int = 3,
    degree_v: int = 3,
    n_points_per_section: int | None = None,
    mesh_n_u: int = 40,
    mesh_n_v: int = 14,
    waterline_z: float = 0.0,
    mirror_y: bool = True,
    close_waterplane: bool = False,
    interpolate: bool = True,
) -> ParametricHull:
    """Build a complete parametric hull from cross-sections.

    This is the main entry point: loft a NURBS surface through the sections,
    tessellate it, compute hydrostatics on the wetted hull, and assess mesh quality on
    the hull triangles (closure lids excluded).

    Returns
    -------
    ParametricHull
        Sections, interpolating NURBS surface, surface mesh, hydrostatics and quality
        report for this hull.
    """
    xs = [s.x for s in sections]
    if length is None:
        length = max(xs) - min(xs)
    if beam is None:
        beam = 2.0 * max(s.max_half_breadth for s in sections)
    if draft is None:
        min_z = min(s.min_z for s in sections)
        draft = abs(min_z) if min_z < 0 else max(s.max_height for s in sections)

    surface = loft_nurbs_surface(
        sections,
        degree_u=degree_u,
        degree_v=degree_v,
        n_points_per_section=n_points_per_section,
        interpolate=interpolate,
    )
    surface_mesh = tessellate_surface(
        surface,
        n_u=mesh_n_u,
        n_v=mesh_n_v,
        mirror_y=mirror_y,
        close_waterplane=close_waterplane,
    )
    hydro = compute_hydrostatics(
        surface_mesh.mesh,
        waterline_z=waterline_z,
        length=length,
        beam=beam,
        draft=draft,
        n_hull_faces=surface_mesh.hull_face_count,
    )
    quality = mesh_quality_report(
        surface_mesh.mesh, n_faces=surface_mesh.hull_face_count
    )
    quality["closure_faces"] = surface_mesh.n_faces - surface_mesh.hull_face_count

    params = {
        "length": length,
        "beam": beam,
        "draft": draft,
        "n_sections": len(sections),
        "degree_u": degree_u,
        "degree_v": degree_v,
        "mesh_n_u": mesh_n_u,
        "mesh_n_v": mesh_n_v,
        "interpolating_loft": float(1 if interpolate else 0),
        "close_waterplane": float(1 if close_waterplane else 0),
    }

    return ParametricHull(
        sections=tuple(sections),
        surface=surface,
        mesh_data=surface_mesh,
        hydrostatics=hydro,
        quality=quality,
        parameters=params,
    )
