"""CAD -> high-order curved-element handoff for TensorDG (``curved_geometry.v1``).

Where this sits
---------------
:mod:`tensorcad.adapters.fem` hands a surface to a *linear* finite-element
preprocessor, :mod:`tensorcad.volume_mesh` builds a straight-sided volume mesh and
:mod:`tensorcad.solver_interface` writes the ``handoff/1`` bundle with the z-up
frame contract.  None of those carry *curved* element geometry: once a solver
integrates over an isoparametric cell whose wall is a polynomial curve, the wall
has to be described by control points, weights and knots, not by a triangle
table.  That is the one solver-facing route the CAD side never produced, and it
is the third leg of ``GAP-CAD-B3-SOLVER-INTERFACE`` ("watertight -> FVM volume /
FEM / **DG high-order curved** / LBM voxel").

This module produces exactly that description, from a mesh the package already
generated, and verifies it *before* a solver ever sees it:

``section_contour``
    Slice a watertight :class:`~tensorcad.geometry.Mesh` with the plane
    ``z = const`` (the z-up contract's horizontal plane) and return the closed
    polyline the surface cuts.  Manifoldness is checked: every section point must
    join exactly two segments, otherwise the cut is refused.
``fit_circle`` / ``rational_circle_contract``
    Fit the contour's centre and radius (algebraic fit + Gauss-Newton refinement)
    and emit the rational-quadratic NURBS circle as a
    ``tensordg.curved_geometry.v1`` contract.  The fit is *measured*, not
    assumed: the residual of the section points against the fitted circle is
    reported so the caller sees the faceting error it inherits.
``bspline_closed_contract``
    Least-squares clamped closed B-spline for a section that is not a circle
    (an ellipse, a hull station): the honest general path, with its own
    deviation report.
``validate_contract`` / ``curve_report``
    The consumer's rules are re-checked here with an *independent* de Boor
    evaluator (nothing in this module imports TensorDG), and the classic
    analytic invariants of the emitted curve are measured: closure gap,
    radial deviation, enclosed area against ``pi R^2`` (Green's theorem) and
    arc length against ``2 pi R``, each at two sampling densities so the
    reader can see that the residual is sampling-limited and not geometry-
    limited.
``write_handoff``
    ``contract.json`` + ``handoff.json`` (units, frame, source mesh digest,
    measured metrics) in one directory, so an artefact on disk says what it is.

What is *not* claimed
---------------------
The general (non-circular) contract is produced and verified here, but the
consumer side that exists today (``tensordg.curved_geometry.mesh_from_contract``)
builds a polar annulus and therefore only accepts a **circular** inner wall.
:func:`consumability` reports that boundary explicitly instead of pretending the
ellipse is consumable; see ``docs/dg_curved_handoff.md``.
"""
from __future__ import annotations

import hashlib
import json
import math
from collections.abc import Sequence
from pathlib import Path
from typing import Any

from .geometry import Mesh

__all__ = [
    "DG_CONTRACT_SCHEMA",
    "HANDOFF_SCHEMA",
    "ContractError",
    "bspline_closed_contract",
    "curve_report",
    "fit_circle",
    "nurbs_curve_points",
    "rational_circle_contract",
    "section_contour",
    "validate_contract",
    "write_handoff",
]

#: The consumer's schema string, quoted rather than imported (see module docstring).
DG_CONTRACT_SCHEMA = "tensordg.curved_geometry.v1"

#: Bumped when the *CAD-side* handoff record changes shape.
HANDOFF_SCHEMA = "tensorcad.dg_curved/1"

#: Exact rational-quadratic circle: weights alternate ``1, 1/sqrt(2)``.
_SQRT_HALF = 1.0 / math.sqrt(2.0)
_UNIT_DIRECTIONS = (
    (1.0, 0.0), (1.0, 1.0), (0.0, 1.0), (-1.0, 1.0), (-1.0, 0.0),
    (-1.0, -1.0), (0.0, -1.0), (1.0, -1.0), (1.0, 0.0),
)


class ContractError(ValueError):
    """A contract or a sliced section violates the handoff rules."""


# ---------------------------------------------------------------------------
# Slicing a watertight surface
# ---------------------------------------------------------------------------
def _quantise(point: tuple[float, float], scale: float) -> tuple[int, int]:
    return (int(round(point[0] / scale)), int(round(point[1] / scale)))


def section_contour(
    mesh: Mesh, z: float, *, tol: float = 1e-9 * 1.0
) -> list[tuple[float, float]]:
    """Closed polyline where the plane ``z = const`` cuts a watertight mesh.

    Only triangles that straddle the plane contribute; the plane may not pass
    through a vertex (that is an ambiguous cut and is refused, not guessed), and
    every resulting section point must be shared by exactly two segments --
    otherwise the surface is not closed across the plane and no contour exists.
    The returned ring is a list of distinct points in traversal order; the first
    point is *not* repeated at the end.
    """
    vertices = mesh.vertices
    zs = [vertex[2] for vertex in vertices]
    z_min, z_max = min(zs), max(zs)
    if not (z_min < z < z_max):
        raise ContractError(
            f"section plane z={z} is outside the mesh extent [{z_min}, {z_max}]"
        )

    segments: list[tuple[tuple[float, float], tuple[float, float]]] = []
    for face in mesh.faces:
        pts: list[tuple[float, float]] = []
        for a, b in ((face[0], face[1]), (face[1], face[2]), (face[2], face[0])):
            za, zb = zs[a] - z, zs[b] - z
            if abs(za) <= tol or abs(zb) <= tol:
                raise ContractError(
                    "section plane passes through a mesh vertex; shift z so the "
                    "cut is unambiguous"
                )
            if (za > 0.0) == (zb > 0.0):
                continue
            t = za / (za - zb)
            va, vb = vertices[a], vertices[b]
            pts.append((va[0] + t * (vb[0] - va[0]), va[1] + t * (vb[1] - va[1])))
        if len(pts) == 2:
            segments.append((pts[0], pts[1]))
        elif len(pts) != 0:  # pragma: no cover - defensive
            raise ContractError(f"triangle crossed the plane {len(pts)} times")
    if not segments:
        raise ContractError("no triangle crosses the section plane")

    # Chain the segments through quantised endpoints.
    scale = 1e-12 + 1e-9 * max(
        max(abs(p[0]), abs(p[1])) for seg in segments for p in seg
    )
    keys = [_quantise(p, scale) for seg in segments for p in seg]
    lookup: dict[tuple[int, int], tuple[float, float]] = {}
    incident: dict[tuple[int, int], list[tuple[int, int]]] = {}
    for index, seg in enumerate(segments):
        ka, kb = keys[2 * index], keys[2 * index + 1]
        if ka == kb:
            raise ContractError("degenerate section segment (zero length)")
        lookup[ka] = seg[0]
        lookup[kb] = seg[1]
        incident.setdefault(ka, []).append(kb)
        incident.setdefault(kb, []).append(ka)
    for key, neighbours in incident.items():
        if len(neighbours) != 2:
            raise ContractError(
                f"section point {lookup[key]} joins {len(neighbours)} segments; "
                "the mesh is not closed across the plane"
            )

    start = min(incident)
    ring_keys = [start]
    previous, current = None, start
    while True:
        options = incident[current]
        nxt = options[0] if options[0] != previous else options[1]
        if nxt == start:
            break
        ring_keys.append(nxt)
        previous, current = current, nxt
        if len(ring_keys) > len(incident):  # pragma: no cover - defensive
            raise ContractError("section ring did not close")
    if len(ring_keys) != len(incident):
        raise ContractError(
            f"section produced {len(ring_keys)} of {len(incident)} points: "
            "the cut is more than one closed loop"
        )
    return [lookup[key] for key in ring_keys]


# ---------------------------------------------------------------------------
# Circle fit
# ---------------------------------------------------------------------------
def _solve(matrix: list[list[float]], rhs: list[float]) -> list[float]:
    """Gaussian elimination with partial pivoting (stdlib only, small systems)."""
    n = len(rhs)
    a = [list(row) + [rhs[i]] for i, row in enumerate(matrix)]
    for col in range(n):
        pivot = max(range(col, n), key=lambda r: abs(a[r][col]))
        if abs(a[pivot][col]) < 1e-300:
            raise ContractError("singular normal-equation system in the fit")
        a[col], a[pivot] = a[pivot], a[col]
        for row in range(col + 1, n):
            factor = a[row][col] / a[col][col]
            if factor:
                for k in range(col, n + 1):
                    a[row][k] -= factor * a[col][k]
    out = [0.0] * n
    for row in range(n - 1, -1, -1):
        acc = a[row][n] - sum(a[row][k] * out[k] for k in range(row + 1, n))
        out[row] = acc / a[row][row]
    return out


def fit_circle(points: Sequence[Sequence[float]]) -> dict[str, float]:
    """Algebraic (Kasa) circle fit refined by Gauss-Newton on the true distance.

    Returns ``centre_x``, ``centre_y``, ``radius``, ``max_residual`` and
    ``rms_residual`` -- the residuals are *point-to-circle* distances, so a
    caller can compare them with the faceting error of the mesh the section came
    from and see that the fit is not the limiting factor.
    """
    pts = [(float(p[0]), float(p[1])) for p in points]
    if len(pts) < 3:
        raise ContractError("a circle fit needs at least three points")
    n = len(pts)
    mx = sum(p[0] for p in pts) / n
    my = sum(p[1] for p in pts) / n
    xs = [p[0] - mx for p in pts]
    ys = [p[1] - my for p in pts]
    # Kasa: minimise sum (x^2 + y^2 - 2 cx x - 2 cy y + k)^2
    sxx = sum(x * x for x in xs)
    syy = sum(y * y for y in ys)
    sxy = sum(x * y for x, y in zip(xs, ys))
    sx = sum(xs)
    sy = sum(ys)
    sxz = sum(x * (x * x + y * y) for x, y in zip(xs, ys))
    syz = sum(y * (x * x + y * y) for x, y in zip(xs, ys))
    sz = sum(x * x + y * y for x, y in zip(xs, ys))
    a = [[2.0 * sxx, 2.0 * sxy, -sx],
         [2.0 * sxy, 2.0 * syy, -sy],
         [2.0 * sx, 2.0 * sy, -n]]
    cx, cy, k = _solve(a, [sxz, syz, sz])
    centre = [cx + mx, cy + my]
    radius = math.sqrt(max(cx * cx + cy * cy - k, 0.0))

    # Gauss-Newton on the geometric residual |p - c| - r.
    for _ in range(24):
        grad = [[0.0, 0.0, 0.0] for _ in range(3)]
        rhs = [0.0, 0.0, 0.0]
        for p in pts:
            dx, dy = p[0] - centre[0], p[1] - centre[1]
            dist = math.hypot(dx, dy)
            if dist < 1e-300:
                continue
            jac = (-dx / dist, -dy / dist, -1.0)
            res = dist - radius
            for i in range(3):
                rhs[i] -= jac[i] * res
                for j in range(3):
                    grad[i][j] += jac[i] * jac[j]
        step = _solve(grad, rhs)
        centre[0] += step[0]
        centre[1] += step[1]
        radius += step[2]
        if max(abs(s) for s in step) < 1e-15 * max(radius, 1.0):
            break

    residuals = [
        math.hypot(p[0] - centre[0], p[1] - centre[1]) - radius for p in pts
    ]
    return {
        "centre_x": centre[0],
        "centre_y": centre[1],
        "radius": radius,
        "max_residual": max(abs(r) for r in residuals),
        "rms_residual": math.sqrt(sum(r * r for r in residuals) / n),
        "n_points": float(n),
    }


# ---------------------------------------------------------------------------
# NURBS: construction (stdlib) and independent de Boor evaluation
# ---------------------------------------------------------------------------
def rational_circle_contract(
    centre: Sequence[float],
    radius: float,
    *,
    name: str = "cad_section",
    source: dict[str, Any] | None = None,
) -> dict[str, Any]:
    """Exact rational-quadratic NURBS circle as a solver-facing contract.

    Nine control points with weights ``1, 1/sqrt2`` reproduce the circle
    *exactly* (not to a tolerance) -- the underlying conic is a NURBS, so this is
    a representation, not a fit.  Units and frame are declared in the contract
    because a solver that receives geometry without them cannot know either.
    """
    cx, cy = float(centre[0]), float(centre[1])
    r = float(radius)
    if not (math.isfinite(cx) and math.isfinite(cy) and math.isfinite(r)) or r <= 0.0:
        raise ContractError("centre must be finite and radius strictly positive")
    contract: dict[str, Any] = {
        "schema": DG_CONTRACT_SCHEMA,
        "kind": "nurbs_curve",
        "name": name,
        "units": "m",
        "frame": "z-up",
        "centre": [cx, cy],
        "degree": 2,
        "knots": [0.0, 0.0, 0.0, 0.25, 0.25, 0.5, 0.5, 0.75, 0.75, 1.0, 1.0, 1.0],
        "control_points": [[cx + r * ux, cy + r * uy] for ux, uy in _UNIT_DIRECTIONS],
        "weights": [1.0, _SQRT_HALF, 1.0, _SQRT_HALF, 1.0, _SQRT_HALF, 1.0,
                    _SQRT_HALF, 1.0],
        "closed": True,
    }
    if source is not None:
        contract["source"] = dict(source)
    return contract


def _find_span(n: int, degree: int, knots: Sequence[float], t: float) -> int:
    if t >= knots[n]:
        return n - 1
    if t <= knots[degree]:
        return degree
    lo, hi = degree, n
    mid = (lo + hi) // 2
    while t < knots[mid] or t >= knots[mid + 1]:
        if t < knots[mid]:
            hi = mid
        else:
            lo = mid
        mid = (lo + hi) // 2
    return mid


def _basis_funs(span: int, t: float, degree: int, knots: Sequence[float]) -> list[float]:
    values = [0.0] * (degree + 1)
    left = [0.0] * (degree + 1)
    right = [0.0] * (degree + 1)
    values[0] = 1.0
    for j in range(1, degree + 1):
        left[j] = t - knots[span + 1 - j]
        right[j] = knots[span + j] - t
        saved = 0.0
        for r in range(j):
            denom = right[r + 1] + left[j - r]
            temp = values[r] / denom if denom != 0.0 else 0.0
            values[r] = saved + right[r + 1] * temp
            saved = left[j - r] * temp
        values[j] = saved
    return values


def nurbs_curve_points(
    contract: dict[str, Any], t: Sequence[float]
) -> list[tuple[float, float]]:
    """Independent rational de Boor evaluation of a contract's curve.

    Deliberately written here, in the CAD package, over the raw contract fields:
    the numbers a solver derives have to be reproducible without the solver.
    """
    degree = int(contract["degree"])
    knots = [float(k) for k in contract["knots"]]
    cps = [(float(p[0]), float(p[1])) for p in contract["control_points"]]
    wts = [float(w) for w in contract["weights"]]
    n = len(cps)
    out: list[tuple[float, float]] = []
    for tt in t:
        span = _find_span(n, degree, knots, float(tt))
        basis = _basis_funs(span, float(tt), degree, knots)
        num_x = num_y = den = 0.0
        for k, value in enumerate(basis):
            index = span - degree + k
            weight = wts[index] * value
            num_x += weight * cps[index][0]
            num_y += weight * cps[index][1]
            den += weight
        if den == 0.0:  # pragma: no cover - defensive
            raise ContractError(f"degenerate rational basis at t={tt}")
        out.append((num_x / den, num_y / den))
    return out


def _closed_knots(degree: int, n_control: int) -> list[float]:
    interior = n_control - degree - 1
    knots = [0.0] * (degree + 1)
    for j in range(1, interior + 1):
        knots.append(j / (interior + 1))
    knots.extend([1.0] * (degree + 1))
    return knots


def bspline_closed_contract(
    points: Sequence[Sequence[float]],
    *,
    degree: int = 3,
    n_control: int = 12,
    name: str = "cad_section_lsq",
    source: dict[str, Any] | None = None,
) -> dict[str, Any]:
    """Least-squares clamped closed B-spline through a non-circular section.

    The seam sits at ``points[0]``: the first and last control points are pinned
    to it (the consumer requires a clamped knot vector with identical first/last
    control points for a closed curve), and the interior control points are
    solved from the section samples in the least-squares sense.  The result is a
    genuine approximation, so the deviation it achieves is measured in
    :func:`curve_report` rather than assumed small.
    """
    pts = [(float(p[0]), float(p[1])) for p in points]
    if len(pts) < n_control:
        raise ContractError(
            f"need at least n_control={n_control} section points, got {len(pts)}")
    if degree < 1 or n_control <= degree:
        raise ContractError("need degree >= 1 and n_control > degree")
    knots = _closed_knots(degree, n_control)

    # Chord-length parameterisation of the ring, seam at points[0].
    lengths = [0.0]
    for a, b in zip(pts, pts[1:]):
        lengths.append(lengths[-1] + math.hypot(b[0] - a[0], b[1] - a[1]))
    total = lengths[-1]
    params = [length / total for length in lengths]

    n_unknown = n_control - 2
    basis_rows: list[list[float]] = []
    targets: list[tuple[float, float]] = []
    seam_basis: list[float] = []
    n = n_control
    for p, t in zip(pts, params):
        span = _find_span(n, degree, knots, t)
        basis = _basis_funs(span, t, degree, knots)
        row = [0.0] * n
        for k, value in enumerate(basis):
            row[span - degree + k] = value
        basis_rows.append(row)
        seam_basis.append(row[0] + row[-1])
        targets.append(p)

    # Normal equations for the interior control points only.
    gram = [[0.0] * n_unknown for _ in range(n_unknown)]
    rhs_x = [0.0] * n_unknown
    rhs_y = [0.0] * n_unknown
    for row, seam, (px, py) in zip(basis_rows, seam_basis, targets):
        reduced = row[1:-1]
        for i in range(n_unknown):
            if reduced[i] == 0.0:
                continue
            rhs_x[i] += reduced[i] * (px - seam * px)
            rhs_y[i] += reduced[i] * (py - seam * py)
            for j in range(n_unknown):
                gram[i][j] += reduced[i] * reduced[j]
    sol_x = _solve(gram, rhs_x)
    sol_y = _solve(gram, rhs_y)
    control = [list(pts[0])] + [[sol_x[i], sol_y[i]] for i in range(n_unknown)]
    control.append(list(pts[0]))

    contract: dict[str, Any] = {
        "schema": DG_CONTRACT_SCHEMA,
        "kind": "nurbs_curve",
        "name": name,
        "units": "m",
        "frame": "z-up",
        "degree": degree,
        "knots": knots,
        "control_points": control,
        "weights": [1.0] * n_control,
        "closed": True,
    }
    if source is not None:
        contract["source"] = dict(source)
    return contract


# ---------------------------------------------------------------------------
# Validation (the consumer's rules, plus the CAD-side declarations)
# ---------------------------------------------------------------------------
def validate_contract(contract: dict[str, Any]) -> dict[str, Any]:
    """Check a contract against the consumer's rules *and* the CAD-side ones.

    Beyond the structural rules the solver applies (schema, kind, degree, knot
    count, clamped and non-decreasing knots, positive weights, closed-curve
    endpoints), a CAD contract must declare its units and frame: geometry that
    does not say whether it is metres or millimetres, and which way is up, is
    refused here rather than mis-consumed downstream.
    """
    if not isinstance(contract, dict):
        raise ContractError("contract must be a mapping")
    if contract.get("schema") != DG_CONTRACT_SCHEMA:
        raise ContractError(
            f"unsupported schema {contract.get('schema')!r}; expected "
            f"{DG_CONTRACT_SCHEMA!r}")
    if contract.get("kind") != "nurbs_curve":
        raise ContractError(f"unsupported geometry kind {contract.get('kind')!r}")
    if contract.get("units") != "m":
        raise ContractError(
            f"contract must declare metres (units='m'), got {contract.get('units')!r}")
    frame = contract.get("frame")
    if frame != "z-up":
        raise ContractError(
            f"contract must declare the z-up frame, got frame={frame!r}")
    try:
        degree = int(contract["degree"])
    except KeyError as exc:
        raise ContractError(f"missing field {exc.args[0]!r}") from None
    except (TypeError, ValueError) as exc:
        raise ContractError(f"degree is not an integer: {exc}") from None
    if degree < 1:
        raise ContractError("degree must be >= 1")
    for key in ("knots", "control_points", "weights"):
        if key not in contract:
            raise ContractError(f"missing field {key!r}")
    knots = [float(k) for k in contract["knots"]]
    cps = [[float(p[0]), float(p[1])] for p in contract["control_points"]]
    weights = [float(w) for w in contract["weights"]]
    for value in knots + [c for p in cps for c in p] + weights:
        if not math.isfinite(value):
            raise ContractError("non-finite knot, control point or weight")
    if len(knots) != len(cps) + degree + 1:
        raise ContractError(
            f"knot count {len(knots)} != n_control_points + degree + 1 = "
            f"{len(cps) + degree + 1}")
    if any(b < a for a, b in zip(knots, knots[1:])):
        raise ContractError("knot vector must be non-decreasing")
    if knots[0] >= knots[-1]:
        raise ContractError("knot vector must span a non-empty parameter interval")
    if any(abs(k - knots[0]) > 0.0 for k in knots[: degree + 1]) or any(
            abs(k - knots[-1]) > 0.0 for k in knots[-(degree + 1):]):
        raise ContractError("only clamped (open) knot vectors are supported")
    if len(weights) != len(cps):
        raise ContractError("weights must have one entry per control point")
    if any(w <= 0.0 for w in weights):
        raise ContractError("all weights must be strictly positive")
    if contract.get("closed", False):
        if max(abs(cps[0][i] - cps[-1][i]) for i in range(2)) > 1e-12:
            raise ContractError(
                "closed curve requires identical first/last control point")
    if "centre" in contract:
        centre = [float(c) for c in contract["centre"]]
        if len(centre) != 2 or not all(math.isfinite(c) for c in centre):
            raise ContractError("centre must be a finite 2-vector")
    return contract


# ---------------------------------------------------------------------------
# Measured geometry report
# ---------------------------------------------------------------------------
def _rings(contract: dict[str, Any], samples: int) -> list[tuple[float, float]]:
    return nurbs_curve_points(contract, [i / samples for i in range(samples + 1)])


def _green_area(ring: Sequence[tuple[float, float]]) -> float:
    return 0.5 * sum(
        a[0] * b[1] - b[0] * a[1] for a, b in zip(ring, ring[1:])
    )


def _perimeter(ring: Sequence[tuple[float, float]]) -> float:
    return sum(math.hypot(b[0] - a[0], b[1] - a[1]) for a, b in zip(ring, ring[1:]))


def curve_report(
    contract: dict[str, Any],
    *,
    samples: int = 4096,
    analytic_centre: Sequence[float] | None = None,
    analytic_radius: float | None = None,
    section_points: Sequence[Sequence[float]] | None = None,
) -> dict[str, Any]:
    """Measure the emitted curve: closure, area, perimeter, deviation.

    ``samples`` sets the evaluation density; the report includes the same
    quantities at ``samples // 4`` so the residual can be attributed: an O(h^2)
    drop by a factor ~16 between the two densities means the number is limited by
    the polygonal sampling, not by the geometry.  When ``section_points`` is
    given, the distance from every sliced point to the emitted curve is measured
    too -- that is the faceting error the CAD section carries in.
    """
    validate_contract(contract)
    ring = _rings(contract, samples)
    coarse = _rings(contract, max(samples // 4, 8))
    centre = analytic_centre
    if centre is None and "centre" in contract:
        centre = contract["centre"]
    radius = analytic_radius
    report: dict[str, Any] = {
        "schema": HANDOFF_SCHEMA,
        "contract_schema": contract["schema"],
        "degree": int(contract["degree"]),
        "n_control": len(contract["control_points"]),
        "samples": samples,
        "closure_gap": math.hypot(
            ring[0][0] - ring[-1][0], ring[0][1] - ring[-1][1]),
        "area_m2": _green_area(ring),
        "perimeter_m": _perimeter(ring),
        "area_m2_coarse": _green_area(coarse),
        "perimeter_m_coarse": _perimeter(coarse),
    }
    if centre is not None:
        cx, cy = float(centre[0]), float(centre[1])
        radii = [math.hypot(p[0] - cx, p[1] - cy) for p in ring]
        report["radius_mean_m"] = sum(radii) / len(radii)
        report["radius_spread_m"] = max(radii) - min(radii)
        # Does the emitted curve reproduce the circle it *declares*?  That is the
        # representation invariant and it is round-off, because a rational
        # quadratic NURBS circle is the circle.  The deviation from the analytic
        # radius below is a *different* number: it carries the faceting error of
        # the CAD section that the contract was fitted to.
        report["declared_radius_m"] = report["radius_mean_m"]
        report["self_consistency_max_m"] = max(radii) - min(radii)
        if radius is not None:
            report["radial_deviation_max_m"] = max(abs(r - radius) for r in radii)
            report["area_reference_m2"] = math.pi * radius * radius
            report["perimeter_reference_m"] = 2.0 * math.pi * radius
            report["area_rel_error"] = (
                report["area_m2"] / report["area_reference_m2"] - 1.0)
            report["perimeter_rel_error"] = (
                report["perimeter_m"] / report["perimeter_reference_m"] - 1.0)
            report["area_rel_error_coarse"] = (
                report["area_m2_coarse"] / report["area_reference_m2"] - 1.0)
    if section_points is not None:
        dense = _rings(contract, 4096)
        worst = 0.0
        for p in section_points:
            worst = max(worst, min(
                math.hypot(p[0] - q[0], p[1] - q[1]) for q in dense))
        report["section_point_deviation_max_m"] = worst
        # Compare with the section-to-*ideal-circle* residual when known.
        if centre is not None and radius is not None:
            cx, cy = float(centre[0]), float(centre[1])
            report["section_point_circle_residual_max_m"] = max(
                abs(math.hypot(p[0] - cx, p[1] - cy) - radius) for p in section_points)
    payload = json.dumps(contract, sort_keys=True, separators=(",", ":"))
    report["contract_digest"] = hashlib.sha256(payload.encode()).hexdigest()
    return report


# ---------------------------------------------------------------------------
# Writing the handoff
# ---------------------------------------------------------------------------
def write_handoff(outdir: str | Path, contract: dict[str, Any], report: dict[str, Any],
                  *, extra: dict[str, Any] | None = None) -> dict[str, Any]:
    """Write ``contract.json`` + ``handoff.json``; return the handoff record."""
    validate_contract(contract)
    directory = Path(outdir)
    directory.mkdir(parents=True, exist_ok=True)
    (directory / "contract.json").write_text(
        json.dumps(contract, indent=2, sort_keys=True) + "\n")
    record = {
        "schema": HANDOFF_SCHEMA,
        "contract_file": "contract.json",
        "contract_schema": contract["schema"],
        "frame": contract["frame"],
        "units": contract["units"],
        "geometry": report,
    }
    if extra:
        record.update(extra)
    (directory / "handoff.json").write_text(
        json.dumps(record, indent=2, sort_keys=True) + "\n")
    return record