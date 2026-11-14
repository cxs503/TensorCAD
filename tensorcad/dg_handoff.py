"""CAD-side geometry source for TensorDG's curved (isoparametric) elements.

Why this module exists
----------------------
:mod:`tensordg.curved_geometry` consumes a ``tensordg.curved_geometry.v1``
contract (NURBS control points + knots + weights) and turns it into DG element
geometry with ``mesh_from_contract``.  Its own docstring calls that path "the
concrete TensorCAD(N03) -> TensorDG hand-off: CAD supplies control points,
weights and knots, DG samples them at its GLL nodes".  Until now the only
producer of that payload was TensorDG's own ``circle_contract`` -- the solver
was handing geometry to itself, so a CAD error had no way to enter, and a wrong
contract could not be diagnosed as a geometry defect rather than a solver
defect.

This module makes the CAD side authoritative and machine-checked:

``construction``
    :func:`circle_contract` builds the rational-quadratic circle from the
    *tangent-intersection* construction a CAD kernel uses (segment corner =
    intersection of the end tangents, weight ``cos(dtheta/2)``), for any radius,
    centre, segment count and 3D plane.  It is parametric geometry, not a copy
    of the solver's constant table.
``independent evaluation``
    :func:`evaluate_curve` / :func:`evaluate_curve_derivative` implement the
    rational de Boor evaluation (Piegl & Tiller A3.2 / A4.2) in plain double
    precision, written here rather than imported, so a CAD number is never
    produced by the consumer's code path.
``measurement``
    :func:`curve_geometry_report` measures closure gap, maximum radial
    deviation, circumference (Gauss quadrature of ``|C'(t)|`` per knot span) and
    enclosed area (Green's theorem).  The analytic references ``2 pi R`` and
    ``pi R^2`` come from the circle, not from the curve under test.
``frame``
    The DG contract is two-dimensional while the ecosystem frame is z-up
    3-D (:data:`tensorcad.solver_interface.COORDINATE_CONTRACT`).  A
    :class:`PlaneMap` records which two axes the 2-D chart uses, the plane
    offset, and the chart handedness -- ``u x v`` is ``+normal`` for the
    horizontal ``x-y`` chart and ``-normal`` for the centreline ``x-z`` chart.
    Handing over a 2-D table without that record is exactly how a coordinate
    convention gets lost between CAD and a solver.
``verification``
    :func:`verify_bundle` re-reads the written files and recomputes every
    declared number from the stored tables, refusing a bundle whose numbers do
    not survive the re-read.  The 3-D source ring stays in the bundle so the
    3-D -> 2-D projection is re-checked, not remembered.

Limitations (stated, not hidden)
--------------------------------
* The exported wall must be a *circle*: TensorDG's annulus builder refuses
  anything else (radius spread ``> 1e-8 r``), so a generic ship-section curve is
  out of scope until that builder generalises.  :func:`circle_contract` is the
  only constructor; there is no fitting of measured geometry to NURBS.
* The contract is a curve (the wall), not a surface: it exercises the
  curved-wall route of the DG layer, not a curved 3-D volume mesh.
* Everything here is geometry.  Hydrodynamic or structural suitability of the
  exported wall is not claimed and cannot be established at this level.
"""
from __future__ import annotations

import hashlib
import json
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable, Sequence

__all__ = [
    "BUNDLE_SCHEMA",
    "DG_CONTRACT_SCHEMA",
    "PLANES",
    "PlaneMap",
    "circle_contract",
    "check_plane_ring",
    "curve_geometry_report",
    "evaluate_curve",
    "evaluate_curve_derivative",
    "verify_bundle",
    "write_dg_bundle",
]

#: The contract schema TensorDG's ``validate_geometry_contract`` accepts.  Not
#: ours to bump; if TensorDG changes it this module refuses the write instead of
#: emitting a payload nobody validates.
DG_CONTRACT_SCHEMA = "tensordg.curved_geometry.v1"

#: Bumped on any incompatible change to ``handoff.json``.
BUNDLE_SCHEMA = "tensorcad.dg_curved_wall/1"

#: ``(u axis, v axis, normal axis)`` for every 2-D chart we agree to emit.
PLANES: dict[str, tuple[str, str, str]] = {
    "xy": ("x", "y", "z"),   # waterline-parallel (horizontal) chart
    "xz": ("x", "z", "y"),   # centreline (vertical) chart
}

_AXIS_INDEX = {"x": 0, "y": 1, "z": 2}


class HandoffError(ValueError):
    """A handoff that cannot be verified is refused, never silently repaired."""


# --------------------------------------------------------------------------- frame


@dataclass(frozen=True)
class PlaneMap:
    """How a 2-D DG chart sits inside the z-up 3-D contract frame."""

    name: str
    u: str
    v: str
    normal: str
    offset_m: float = 0.0

    def __post_init__(self) -> None:
        if self.name not in PLANES:
            raise HandoffError(f"unknown plane {self.name!r}; expected one of {sorted(PLANES)}")
        expected = PLANES[self.name]
        if (self.u, self.v, self.normal) != expected:
            raise HandoffError(
                f"plane {self.name!r} is fixed to axes {expected}, got "
                f"{(self.u, self.v, self.normal)}")
        if not math.isfinite(self.offset_m):
            raise HandoffError("plane offset must be finite")

    @property
    def chart_handedness(self) -> int:
        """``+1`` if ``u x v`` is ``+normal``, ``-1`` if the chart is mirrored.

        A 2-D solver that integrates counter-clockwise in ``(u, v)`` walks
        clockwise about ``+normal`` when this is ``-1``; that sign decides the
        physical direction of a wall normal and must not be guessed.
        """
        return 1 if self.name == "xy" else -1

    def project(self, point: Sequence[float]) -> tuple[float, float]:
        return (float(point[_AXIS_INDEX[self.u]]), float(point[_AXIS_INDEX[self.v]]))

    def out_of_plane(self, point: Sequence[float]) -> float:
        return float(point[_AXIS_INDEX[self.normal]]) - self.offset_m

    def to_3d(self, uv: Sequence[float]) -> tuple[float, float, float]:
        """Embed a chart point back into the 3-D frame at the plane offset."""
        out = [0.0, 0.0, 0.0]
        out[_AXIS_INDEX[self.u]] = float(uv[0])
        out[_AXIS_INDEX[self.v]] = float(uv[1])
        out[_AXIS_INDEX[self.normal]] = self.offset_m
        return (out[0], out[1], out[2])

    def to_dict(self) -> dict[str, Any]:
        return {
            "name": self.name,
            "u_axis": self.u,
            "v_axis": self.v,
            "normal_axis": self.normal,
            "offset_m": self.offset_m,
            "chart_handedness": self.chart_handedness,
            "note": "2-D DG chart inside the z-up 3-D contract frame",
        }

    @classmethod
    def from_dict(cls, payload: dict[str, Any]) -> "PlaneMap":
        return cls(name=payload["name"], u=payload["u_axis"], v=payload["v_axis"],
                   normal=payload["normal_axis"], offset_m=float(payload["offset_m"]))


# ------------------------------------------------------------- rational evaluation


def _check_contract(contract: dict[str, Any]) -> tuple[int, list[float], list[tuple[float, float]], list[float]]:
    """Reject a malformed contract before any number is produced."""
    if not isinstance(contract, dict):
        raise HandoffError("contract must be a mapping")
    if contract.get("schema") != DG_CONTRACT_SCHEMA:
        raise HandoffError(
            f"unsupported schema {contract.get('schema')!r}; expected {DG_CONTRACT_SCHEMA!r}")
    if contract.get("kind") != "nurbs_curve":
        raise HandoffError(f"unsupported geometry kind {contract.get('kind')!r}")
    try:
        degree = int(contract["degree"])
    except (KeyError, TypeError, ValueError):
        raise HandoffError("contract needs an integer 'degree'") from None
    if degree < 1:
        raise HandoffError("degree must be >= 1")
    knots = [float(k) for k in contract["knots"]]
    cps = [(float(p[0]), float(p[1])) for p in contract["control_points"]]
    weights = [float(w) for w in contract["weights"]]
    if len(knots) != len(cps) + degree + 1:
        raise HandoffError(
            f"knot count {len(knots)} != n_control_points + degree + 1 = {len(cps) + degree + 1}")
    if len(weights) != len(cps):
        raise HandoffError("weights must have one entry per control point")
    if not all(math.isfinite(k) for k in knots) or not all(math.isfinite(c) for p in cps for c in p):
        raise HandoffError("non-finite knot or control point")
    if any(b < a for a, b in zip(knots, knots[1:])):
        raise HandoffError("knot vector must be non-decreasing")
    if not knots[-1] > knots[0]:
        raise HandoffError("knot vector must span a non-empty parameter interval")
    if any(w <= 0.0 for w in weights):
        raise HandoffError("all weights must be strictly positive")
    if not contract.get("closed", False):
        raise HandoffError("a DG curved-wall contract must be a closed curve")
    if max(abs(cps[0][0] - cps[-1][0]), abs(cps[0][1] - cps[-1][1])) > 1e-12:
        raise HandoffError("closed curve requires identical first/last control point")
    return degree, knots, cps, weights


def _homogeneous(cps: Sequence[Sequence[float]], weights: Sequence[float]):
    return [(w * p[0], w * p[1], w) for p, w in zip(cps, weights)]


def _find_span(n_ctrl: int, degree: int, knots: Sequence[float], u: float) -> int:
    """Knot span index (Piegl & Tiller A2.1), clamped to the domain."""
    lo, hi = degree, n_ctrl  # k in [degree, n_ctrl-1] for u in [knots[p], knots[n]]
    span = hi
    for k in range(degree, n_ctrl):
        if knots[k] <= u < knots[k + 1]:
            span = k
            break
    if u >= knots[hi]:
        span = hi - 1
    return max(lo - 1 + 1, min(span, hi - 1))


def _basis(degree: int, knots: Sequence[float], span: int, u: float) -> list[float]:
    """Non-zero B-spline basis values ``N_{span-degree .. span}`` at ``u``."""
    N = [0.0] * (degree + 1)
    left = [0.0] * (degree + 1)
    right = [0.0] * (degree + 1)
    N[0] = 1.0
    for j in range(1, degree + 1):
        left[j] = u - knots[span + 1 - j]
        right[j] = knots[span + j] - u
        saved = 0.0
        for r in range(j):
            den = right[r + 1] + left[j - r]
            temp = N[r] / den if den != 0.0 else 0.0
            N[r] = saved + right[r + 1] * temp
            saved = left[j - r] * temp
        N[j] = saved
    return N


def _eval_homogeneous(degree: int, knots: Sequence[float], H: Sequence[Sequence[float]],
                      u: float) -> tuple[float, float, float]:
    n_ctrl = len(H)
    span = _find_span(n_ctrl, degree, knots, u)
    N = _basis(degree, knots, span, u)
    hx = hy = hw = 0.0
    for i in range(degree + 1):
        c = span - degree + i
        hx += N[i] * H[c][0]
        hy += N[i] * H[c][1]
        hw += N[i] * H[c][2]
    return (hx, hy, hw)


def _derivative_homogeneous(degree: int, knots: Sequence[float],
                            H: Sequence[Sequence[float]]):
    """Homogeneous control points of ``d/dt C_h`` and its knot vector.

    ``P'_i = p (P_{i+1} - P_i) / (k_{i+p+1} - k_{i+1})`` on ``k[1:-1]``
    (Piegl & Tiller A3.2), applied in homogeneous coordinates so the quotient
    rule below stays exact for the rational curve.
    """
    p = degree
    Dh = []
    for i in range(len(H) - 1):
        den = knots[i + p + 1] - knots[i + 1]
        if den <= 0.0:
            Dh.append((0.0, 0.0, 0.0))
            continue
        f = p / den
        Dh.append((f * (H[i + 1][0] - H[i][0]),
                   f * (H[i + 1][1] - H[i][1]),
                   f * (H[i + 1][2] - H[i][2])))
    return Dh, list(knots[1:-1])


def evaluate_curve(contract_or_parts, t_values: Iterable[float]) -> list[tuple[float, float]]:
    """Rational de Boor evaluation, written here (not imported from the consumer).

    Accepts either a contract mapping or ``(degree, knots, control_points,
    weights)``.  Returns Euclidean chart coordinates, i.e. the homogeneous
    components divided by ``w(t)``.
    """
    degree, knots, cps, weights = _parts(contract_or_parts)
    H = _homogeneous(cps, weights)
    out = []
    lo, hi = knots[degree], knots[len(cps)]
    for t in t_values:
        u = min(max(float(t), lo), hi)
        hx, hy, hw = _eval_homogeneous(degree, knots, H, u)
        out.append((hx / hw, hy / hw))
    return out


def evaluate_curve_derivative(contract_or_parts, t_values: Iterable[float]
                              ) -> list[tuple[float, float]]:
    """``dC/dt`` of the rational curve via the homogeneous quotient rule."""
    degree, knots, cps, weights = _parts(contract_or_parts)
    H = _homogeneous(cps, weights)
    Dh, dknots = _derivative_homogeneous(degree, knots, H)
    lo, hi = knots[degree], knots[len(cps)]
    out = []
    for t in t_values:
        u = min(max(float(t), lo), hi)
        hx, hy, hw = _eval_homogeneous(degree, knots, H, u)
        cx, cy = hx / hw, hy / hw
        if degree == 1:
            # derivative is a degree-0 (piecewise constant) curve
            span = _find_span(len(Dh), 0, dknots, u)
            dx, dy, dw = Dh[min(span, len(Dh) - 1)]
        else:
            dx, dy, dw = _eval_homogeneous(degree - 1, dknots, Dh, u)
        out.append(((dx - cx * dw) / hw, (dy - cy * dw) / hw))
    return out


def _parts(contract_or_parts):
    if isinstance(contract_or_parts, dict):
        return _check_contract(contract_or_parts)
    degree, knots, cps, weights = contract_or_parts
    return (int(degree), [float(k) for k in knots],
            [(float(p[0]), float(p[1])) for p in cps], [float(w) for w in weights])


# ------------------------------------------------------------------- quadrature


def _gauss_legendre(n: int) -> tuple[list[float], list[float]]:
    """Nodes/weights on ``[-1, 1]`` by Newton iteration (stdlib only)."""
    nodes, weights = [], []
    for i in range(1, n + 1):
        x = math.cos(math.pi * (i - 0.25) / (n + 0.5))
        dp = 1.0
        for _ in range(100):
            p0, p1 = 1.0, x
            for k in range(2, n + 1):
                p0, p1 = p1, ((2 * k - 1) * x * p1 - (k - 1) * p0) / k
            dp = n * (x * p1 - p0) / (x * x - 1.0)
            dx = -p1 / dp
            x += dx
            if abs(dx) < 1e-16:
                break
        nodes.append(x)
        weights.append(2.0 / ((1.0 - x * x) * dp * dp))
    return nodes, weights


def _spans(knots: Sequence[float]) -> list[tuple[float, float]]:
    """Distinct non-empty knot intervals."""
    out = []
    for a, b in zip(knots, knots[1:]):
        if b > a:
            out.append((a, b))
    return out


# ------------------------------------------------------------------ construction


def circle_contract(radius: float, centre: Sequence[float] = (0.0, 0.0), *,
                    segments: int = 4, plane: str = "xy", offset_m: float = 0.0,
                    name: str = "cad-circle") -> dict[str, Any]:
    """Rational-quadratic NURBS circle from the CAD tangent-intersection construction.

    One quadratic segment per ``1/segments`` of the circle: the segment's middle
    control point is the intersection of the end tangents (distance
    ``radius / cos(dtheta/2)`` from the centre) with weight ``cos(dtheta/2)``,
    which makes the rational quadratic *exactly* the circular arc.  Nothing is
    interpolated and no table is copied: the control polygon is derived here from
    ``radius``, ``centre`` and ``segments``.
    """
    r = float(radius)
    if not (r > 0.0 and math.isfinite(r)):
        raise HandoffError("radius must be positive and finite")
    cx, cy = float(centre[0]), float(centre[1])
    if not (math.isfinite(cx) and math.isfinite(cy)):
        raise HandoffError("centre must be finite")
    seg = int(segments)
    if seg < 1:
        raise HandoffError("segments must be >= 1")
    if plane not in PLANES:
        raise HandoffError(f"unknown plane {plane!r}; expected one of {sorted(PLANES)}")
    dtheta = 2.0 * math.pi / seg
    cos_half = math.cos(dtheta / 2.0)
    if cos_half <= 0.0:
        raise HandoffError("segments too few for a closed circular arc")
    cps: list[list[float]] = []
    weights: list[float] = []
    for i in range(seg):
        a0 = i * dtheta
        cps.append([cx + r * math.cos(a0), cy + r * math.sin(a0)])
        weights.append(1.0)
        am = a0 + 0.5 * dtheta
        rc = r / cos_half
        cps.append([cx + rc * math.cos(am), cy + rc * math.sin(am)])
        weights.append(cos_half)
    cps.append(list(cps[0]))
    weights.append(weights[0])
    knots = [0.0, 0.0, 0.0]
    for i in range(1, seg):
        knots.extend([i / seg, i / seg])
    knots.extend([1.0, 1.0, 1.0])
    return {
        "schema": DG_CONTRACT_SCHEMA,
        "kind": "nurbs_curve",
        "name": name,
        "centre": [cx, cy],
        "degree": 2,
        "knots": knots,
        "control_points": cps,
        "weights": weights,
        "closed": True,
        "radius_m": r,
        "plane": PlaneMap(plane, *PLANES[plane], offset_m).to_dict(),
        "producer": "tensorcad.dg_handoff.circle_contract (tangent intersection)",
    }


# ----------------------------------------------------------------- measurement


def _contract_digest(contract: dict[str, Any]) -> str:
    degree, knots, cps, weights = _check_contract(contract)
    h = hashlib.sha256()
    h.update(f"{degree}|{len(cps)}|".encode())
    h.update((",".join(f"{k:.17g}" for k in knots) + ";").encode())
    for (x, y), w in zip(cps, weights):
        h.update((f"{x:.17g},{y:.17g},{w:.17g};").encode())
    return h.hexdigest()


def curve_geometry_report(contract: dict[str, Any], *, quad_order: int = 20) -> dict[str, Any]:
    """Measure the curve against the circle it is supposed to be.

    Every number is computed from the stored tables by :func:`evaluate_curve`
    and :func:`evaluate_curve_derivative`; the references (``2 pi R``,
    ``pi R^2``) are analytic and independent of the curve under test.
    """
    degree, knots, cps, weights = _check_contract(contract)
    radius = float(contract.get("radius_m", 0.0))
    centre = contract.get("centre")
    if radius <= 0.0:
        # fall back to a measured radius so a foreign contract can still be read
        sample = evaluate_curve(contract, [knots[degree] + (knots[len(cps)] - knots[degree]) * i / 64
                                           for i in range(65)])
        cx, cy = (float(centre[0]), float(centre[1])) if centre else (0.0, 0.0)
        radius = sum(math.hypot(p[0] - cx, p[1] - cy) for p in sample) / len(sample)
    elif centre:
        cx, cy = float(centre[0]), float(centre[1])
    else:
        raise HandoffError("a contract with radius_m must also carry centre")

    nodes, wq = _gauss_legendre(quad_order)
    circumference = 0.0
    area2 = 0.0
    worst_radial = 0.0
    lo, hi = knots[degree], knots[len(cps)]
    for a, b in _spans(knots):
        half = 0.5 * (b - a)
        ts = [0.5 * (a + b) + half * x for x in nodes]
        pts = evaluate_curve(contract, ts)
        der = evaluate_curve_derivative(contract, ts)
        for (x, y), (dx, dy), w in zip(pts, der, wq):
            speed = math.hypot(dx, dy)
            circumference += w * half * speed
            area2 += w * half * (x * dy - y * dx)
            worst_radial = max(worst_radial, abs(math.hypot(x - cx, y - cy) - radius))
    closure = evaluate_curve(contract, [lo, hi])
    closure_gap = math.hypot(closure[1][0] - closure[0][0], closure[1][1] - closure[0][1])
    area = 0.5 * abs(area2)
    analytic_circumference = 2.0 * math.pi * radius
    analytic_area = math.pi * radius * radius
    return {
        "degree": degree,
        "segments": len(_spans(knots)),
        "control_points": len(cps),
        "radius_m": radius,
        "centre_m": [cx, cy],
        "closure_gap_m": closure_gap,
        "max_radial_deviation_m": worst_radial,
        "max_radial_deviation_rel": worst_radial / radius,
        "circumference_m": circumference,
        "circumference_analytic_m": analytic_circumference,
        "circumference_rel_error": abs(circumference - analytic_circumference) / analytic_circumference,
        "enclosed_area_m2": area,
        "enclosed_area_analytic_m2": analytic_area,
        "area_rel_error": abs(area - analytic_area) / analytic_area,
        "quadrature_order": quad_order,
    }


def check_plane_ring(points_3d: Sequence[Sequence[float]], plane: PlaneMap, *,
                     radius: float, centre: Sequence[float],
                     out_of_plane_tol_m: float = 1e-9) -> dict[str, Any]:
    """Check that a 3-D CAD ring projects onto the declared 2-D chart as expected.

    This is the frame gate: a mirrored or rotated chart shows up as a large
    out-of-plane component or as a radial error, not as a silent offset in a
    solver result.
    """
    if len(points_3d) < 3:
        raise HandoffError("need at least three 3-D points to check a plane")
    out_of_plane = max(abs(plane.out_of_plane(p)) for p in points_3d)
    proj = [plane.project(p) for p in points_3d]
    cx, cy = float(centre[0]), float(centre[1])
    radial = [abs(math.hypot(u - cx, v - cy) - radius) for u, v in proj]
    worst = max(radial)
    ok = out_of_plane <= out_of_plane_tol_m and worst <= 1e-8 * radius
    return {
        "plane": plane.name,
        "chart_handedness": plane.chart_handedness,
        "n_points": len(points_3d),
        "out_of_plane_max_m": out_of_plane,
        "wall_radial_error_max_m": worst,
        "wall_radial_error_rel": worst / radius,
        "within_tolerance": bool(ok),
        "tolerances": {"out_of_plane_m": out_of_plane_tol_m, "radial_rel": 1e-8},
    }


# -------------------------------------------------------------------- bundle I/O


def write_dg_bundle(outdir: str | Path, contract: dict[str, Any], *, r_outer_m: float,
                    plane_map: PlaneMap | None = None,
                    source_3d_ring: Sequence[Sequence[float]] | None = None,
                    case: str = "", extra: dict[str, Any] | None = None) -> dict[str, Any]:
    """Write the DG curved-wall bundle: contract table + declared numbers.

    ``curved_wall.json`` is what TensorDG reads; ``handoff.json`` is what a
    reader (or :func:`verify_bundle`) checks.  Both are plain JSON so the
    contract is a file, never a private import.
    """
    out = Path(outdir)
    out.mkdir(parents=True, exist_ok=True)
    degree, knots, cps, weights = _check_contract(contract)
    report = curve_geometry_report(contract)
    radius = report["radius_m"]
    if not float(r_outer_m) > radius:
        raise HandoffError(
            f"r_outer_m={r_outer_m!r} must exceed the wall radius {radius!r}")
    if plane_map is None:
        payload = contract.get("plane")
        plane_map = (PlaneMap.from_dict(payload) if payload
                     else PlaneMap("xy", *PLANES["xy"], 0.0))
    plane_check = None
    ring = None
    if source_3d_ring is not None:
        ring = [[float(c) for c in p] for p in source_3d_ring]
        plane_check = check_plane_ring(ring, plane_map, radius=radius,
                                       centre=report["centre_m"])
        if not plane_check["within_tolerance"]:
            raise HandoffError(
                "the 3-D CAD ring does not project onto the declared chart: "
                f"out-of-plane {plane_check['out_of_plane_max_m']:.3e} m, "
                f"radial error {plane_check['wall_radial_error_rel']:.3e} rel")

    (out / "curved_wall.json").write_text(
        json.dumps(contract, indent=2, sort_keys=True) + "\n")
    manifest = {
        "schema": BUNDLE_SCHEMA,
        "consumer_schema": DG_CONTRACT_SCHEMA,
        "consumer": "tensordg.curved_geometry.mesh_from_contract",
        "case": case,
        "tables": {"curved_wall_json": "curved_wall.json"},
        "plane": plane_map.to_dict(),
        "geometry": {
            "digest": _contract_digest(contract),
            "degree": degree,
            "control_points": len(cps),
            "weights": len(weights),
            "closed": bool(contract.get("closed", False)),
            "radius_m": radius,
            "centre_m": report["centre_m"],
            "r_inner_m": radius,
            "r_outer_m": float(r_outer_m),
        },
        "measured": report,
        "plane_check": plane_check,
        "source_3d_ring": ring,
        "units": "m",
        "extra": extra or {},
    }
    (out / "handoff.json").write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    return manifest


def verify_bundle(outdir: str | Path) -> dict[str, Any]:
    """Re-read the bundle and recompute every declared number.

    Raises :class:`HandoffError` if the contract no longer evaluates to what the
    handoff declared, if the digest changed, or if the recorded 3-D ring no
    longer projects onto the declared chart.
    """
    out = Path(outdir)
    handoff_path = out / "handoff.json"
    wall_path = out / "curved_wall.json"
    if not handoff_path.is_file() or not wall_path.is_file():
        raise HandoffError(f"bundle at {out} is missing handoff.json or curved_wall.json")
    manifest = json.loads(handoff_path.read_text())
    contract = json.loads(wall_path.read_text())
    if manifest.get("schema") != BUNDLE_SCHEMA:
        raise HandoffError(f"unsupported bundle schema {manifest.get('schema')!r}")
    if manifest.get("consumer_schema") != DG_CONTRACT_SCHEMA:
        raise HandoffError(
            f"bundle targets {manifest.get('consumer_schema')!r}, this module writes "
            f"{DG_CONTRACT_SCHEMA!r}")

    digest = _contract_digest(contract)
    if digest != manifest["geometry"]["digest"]:
        raise HandoffError("curved_wall.json digest does not match handoff.json")
    report = curve_geometry_report(contract)
    declared = manifest["measured"]
    for key, value in report.items():
        if key not in declared:
            raise HandoffError(f"declared numbers are missing {key!r}")
        if isinstance(value, str):
            if value != declared[key]:
                raise HandoffError(f"{key} changed on re-read: {declared[key]!r} -> {value!r}")
            continue
        if isinstance(value, (list, tuple)):
            if len(value) != len(declared[key]) or any(
                    abs(float(a) - float(b)) > 1e-12 * max(1.0, abs(float(b)))
                    for a, b in zip(value, declared[key])):
                raise HandoffError(f"{key} changed on re-read")
            continue
        if abs(float(value) - float(declared[key])) > 1e-12 * max(1.0, abs(float(declared[key]))):
            raise HandoffError(
                f"{key} changed on re-read: declared {declared[key]!r}, re-measured {value!r}")

    plane = PlaneMap.from_dict(manifest["plane"])
    plane_check = None
    if manifest.get("source_3d_ring"):
        plane_check = check_plane_ring(manifest["source_3d_ring"], plane,
                                       radius=report["radius_m"],
                                       centre=report["centre_m"])
        if not plane_check["within_tolerance"]:
            raise HandoffError("the stored 3-D ring no longer projects onto the declared chart")
        for key in ("out_of_plane_max_m", "wall_radial_error_rel"):
            if abs(plane_check[key] - manifest["plane_check"][key]) > 1e-12:
                raise HandoffError(f"plane check {key} changed on re-read")

    geom = manifest["geometry"]
    if abs(float(geom["radius_m"]) - report["radius_m"]) > 1e-12 * report["radius_m"]:
        raise HandoffError("declared radius does not match the curve")
    if not axis_close(geom["centre_m"], report["centre_m"]):
        raise HandoffError("declared centre does not match the curve")
    if int(geom["degree"]) != report["degree"]:
        raise HandoffError("declared degree does not match the curve")
    if not float(geom["r_outer_m"]) > report["radius_m"]:
        raise HandoffError("declared r_outer_m does not enclose the wall radius")
    return {
        "verified": True,
        "digest": digest,
        "plane": plane.to_dict(),
        "plane_check": plane_check,
        "measured": report,
        "r_inner_m": report["radius_m"],
        "r_outer_m": float(geom["r_outer_m"]),
    }


def axis_close(a: Sequence[float], b: Sequence[float], tol: float = 1e-12) -> bool:
    return len(a) == len(b) and all(
        abs(float(x) - float(y)) <= tol * max(1.0, abs(float(y))) for x, y in zip(a, b))