"""Fit a NURBS B-rep to the platform's ``DTMB5415_Hull_exact.stl`` gold surface.

The gold regression named by GAP-CAD-T1-IGES-PIPE is the platform's own hull mesh
``DTMB5415_Hull_exact.stl``.  It is a triangle soup (no topology), so it cannot be fed
to the IGES import pipeline directly.  This module turns it into an IGES patch set:

1. read the STL triangles (``read_stl``);
2. cut ``n_stations`` transverse sections and resample each on arclength
   (``section_polyline`` / ``girth_grid``);
3. interpolate one hull patch through the station/girth grid (separable B-spline
   interpolation at Greville abscissae, exact at the grid points);
4. close the body with a ruled deck patch and two ruled end caps whose boundaries are
   the *same* NURBS curves as the hull patch edges, so sewing has zero gaps;
5. write entity-128 records with ``geomquality.io.write_iges``.

The fit is an *interpolation* of the gold surface, not a re-invention of it: every
reported deviation is therefore a property of the IGES/import route under test, not of
a different hull definition.  Provenance of the gold file (sha256, face count) is
recorded by the caller in ``artifacts/gap_t1_iges/provenance.json``.
"""
from __future__ import annotations

import struct
from pathlib import Path
from typing import Sequence

import numpy as np

from tensorcad.geomquality.baselines import greville_abscissae
from tensorcad.geomquality.brep import BRep, Face, sew
from tensorcad.nurbs import NURBSSurface, _differentiate_axis, tensor_point

__all__ = [
    "read_stl", "section_polyline", "girth_grid", "interpolate_surface",
    "loft_gold_hull", "brep_quadrature", "graded_panels", "enclosed_volume_flags",
    "PatchDerivatives",
]


# --------------------------------------------------------------------------- #
# STL
# --------------------------------------------------------------------------- #
def read_stl(path: str | Path) -> np.ndarray:
    """Binary STL -> ``(n_triangles, 3, 3)`` float64 array."""
    data = Path(path).read_bytes()
    if len(data) < 84:
        raise ValueError(f"{path}: too short for a binary STL")
    count = struct.unpack("<I", data[80:84])[0]
    if len(data) != 84 + 50 * count:
        raise ValueError(
            f"{path}: {len(data)} bytes does not match a binary STL with {count} triangles"
        )
    raw = np.frombuffer(data[84:84 + 50 * count], dtype=np.uint8).reshape(count, 50)
    vertices = raw[:, 12:48].copy().view("<f4").reshape(count, 3, 3).astype(np.float64)
    return vertices


# --------------------------------------------------------------------------- #
# sections
# --------------------------------------------------------------------------- #
def section_polyline(triangles: np.ndarray, x0: float) -> list[np.ndarray]:
    """All polylines of the intersection of the triangle soup with the plane ``x = x0``.

    Returns the linked chains as ``(k, 3)`` arrays.  Chains are built by welding segment
    end points on a tolerance of ``1e-9`` times the body extent.
    """
    lo = triangles[:, :, 0] - x0
    inside = lo <= 0.0
    crossing = ~(inside.all(axis=1) | (~inside).all(axis=1))
    segments: list[tuple[np.ndarray, np.ndarray]] = []
    for tri, sign in zip(triangles[crossing], lo[crossing]):
        points = []
        for a in range(3):
            b = (a + 1) % 3
            if (sign[a] > 0.0) != (sign[b] > 0.0):
                t = sign[a] / (sign[a] - sign[b])
                points.append(tri[a] + t * (tri[b] - tri[a]))
        if len(points) == 2:
            p, q = points
            segments.append((p, q))
    if not segments:
        return []
    scale = float(np.abs(triangles).max()) or 1.0
    tol = 1e-9 * scale
    buckets: dict[tuple[int, int, int], list[int]] = {}
    keys = []

    def key(point: np.ndarray) -> tuple[int, int, int]:
        return tuple(int(round(value / tol)) for value in point)

    for index, (p, q) in enumerate(segments):
        for point in (p, q):
            buckets.setdefault(key(point), []).append(index)
        keys.append((key(p), key(q)))
    used = [False] * len(segments)
    chains: list[np.ndarray] = []
    for start in range(len(segments)):
        if used[start]:
            continue
        used[start] = True
        chain = [segments[start][0], segments[start][1]]
        grow = True
        while grow:
            grow = False
            for candidate in buckets.get(key(chain[-1]), []):
                if used[candidate]:
                    continue
                p, q = segments[candidate]
                nxt = q if key(p) == key(chain[-1]) else p
                chain.append(nxt)
                used[candidate] = True
                grow = True
                break
        chains.append(np.asarray(chain))
    return chains


def resample(points: np.ndarray, count: int) -> np.ndarray:
    """Resample a polyline on normalised arclength to ``count`` points."""
    steps = np.linalg.norm(np.diff(points, axis=0), axis=1)
    total = steps.sum()
    if total <= 0.0:
        raise ValueError("degenerate polyline")
    cumulative = np.concatenate([[0.0], np.cumsum(steps)]) / total
    targets = np.linspace(0.0, 1.0, count)
    return np.stack(
        [np.interp(targets, cumulative, points[:, axis]) for axis in range(3)], axis=1
    )


def girth_grid(triangles: np.ndarray, stations: int, girth_samples: int,
               ) -> tuple[np.ndarray, np.ndarray, list[float]]:
    """Station x-grid and a ``(stations, girth_samples, 3)`` point grid on the girth.

    The main girth at each station is the chain with the largest transverse span; it is
    oriented so that it starts at the +y deck edge and ends at the -y deck edge.
    """
    x_min = float(triangles[:, :, 0].min())
    x_max = float(triangles[:, :, 0].max())
    length = x_max - x_min
    # cosine spacing clusters stations at the ends, where curvature is highest
    theta = np.linspace(0.0, np.pi, stations)
    xs = x_min + 0.5 * length * (1.0 - np.cos(theta))
    xs = np.clip(xs, x_min + 1e-6 * length, x_max - 1e-6 * length)
    grid = np.zeros((stations, girth_samples, 3))
    spans = []
    for index, x0 in enumerate(xs):
        chains = section_polyline(triangles, float(x0))
        if not chains:
            raise ValueError(f"no section at x={x0}")
        main = max(chains, key=lambda chain: np.ptp(chain[:, 1]))
        if main[0, 1] < main[-1, 1]:
            main = main[::-1]
        spans.append(float(np.ptp(main[:, 1])))
        grid[index] = resample(main, girth_samples)
    return xs, grid, spans


# --------------------------------------------------------------------------- #
# B-spline interpolation
# --------------------------------------------------------------------------- #
def clamped_knots(degree: int, count: int) -> list[float]:
    interior = count - degree - 1
    return [0.0] * (degree + 1) + [i / (interior + 1) for i in range(1, interior + 1)] + [1.0] * (degree + 1)


def _basis(knots: list[float], degree: int, u: float) -> np.ndarray:
    """All non-zero-degree B-spline values at ``u`` (Cox-de Boor, A2.2)."""
    count = len(knots) - degree - 1
    values = np.zeros(count)
    if u >= knots[-1]:
        span = count - 1
    else:
        span = degree
        for i in range(degree, count):
            if knots[i] <= u < knots[i + 1]:
                span = i
                break
    n = np.zeros(degree + 1)
    left = np.zeros(degree + 1)
    right = np.zeros(degree + 1)
    n[0] = 1.0
    for j in range(1, degree + 1):
        left[j] = u - knots[span + 1 - j]
        right[j] = knots[span + j] - u
        saved = 0.0
        for r in range(j):
            denominator = right[r + 1] + left[j - r]
            temp = n[r] / denominator if denominator != 0.0 else 0.0
            n[r] = saved + right[r + 1] * temp
            saved = left[j - r] * temp
        n[j] = saved
    values[span - degree:span + 1] = n
    return values


def interpolate_surface(grid: np.ndarray, degree: tuple[int, int] = (3, 3)) -> NURBSSurface:
    """B-spline surface interpolating ``grid`` (stations x girth) at Greville abscissae.

    Separable: solve ``M_u C M_v^T = P`` with ``M`` the collocation matrices at the
    Greville abscissae (which makes both matrices banded and well conditioned).
    """
    stations, girth, _ = grid.shape
    pu, pv = degree
    knots_u = clamped_knots(pu, stations)
    knots_v = clamped_knots(pv, girth)
    gu = np.asarray(greville_abscissae(knots_u, pu))
    gv = np.asarray(greville_abscissae(knots_v, pv))
    mu = np.stack([_basis(knots_u, pu, float(u)) for u in gu])
    mv = np.stack([_basis(knots_v, pv, float(v)) for v in gv])
    # C = Mu^{-1} P Mv^{-T}
    solved_u = np.linalg.solve(mu, grid.reshape(stations, -1)).reshape(grid.shape)
    net = np.linalg.solve(mv, solved_u.transpose(1, 0, 2).reshape(girth, -1))
    net = net.reshape(girth, stations, 3).transpose(1, 0, 2)
    surface = NURBSSurface(
        degree=(pu, pv),
        control_net=[[tuple(float(value) for value in point) for point in row] for row in net],
        weights=[[1.0] * girth for _ in range(stations)],
        knots=[knots_u, knots_v],
    )
    return surface


def _ruled(edge_a: list[tuple[float, float, float]], edge_b: list[tuple[float, float, float]],
           knots_v: list[float], degree_v: int) -> NURBSSurface:
    """Degree-1-in-u ruled patch between two control polygons with the same v-knots."""
    # control_net is indexed [u][v]: u has two control points (degree 1), v carries the
    # full control polygon of the two edge curves.
    net = [list(edge_a), list(edge_b)]
    return NURBSSurface(
        degree=(1, degree_v),
        control_net=net,
        weights=[[1.0] * len(edge_a) for _ in net],
        knots=[[0.0, 0.0, 1.0, 1.0], list(knots_v)],
    )


def loft_gold_hull(triangles: np.ndarray, *, stations: int = 41, girth_samples: int = 41,
                   degree: int = 3, orient: bool = True) -> tuple[BRep, dict]:
    """Full hull B-rep: girth patch + ruled deck + ruled bow/stern caps.

    ``orient`` derives the ``same_sense`` flags from the enclosed volume
    (:func:`enclosed_volume_flags`)
    instead of guessing them.  The first version of this module hard-coded ``True`` for
    every face, and the resulting patch set was *not* consistently oriented: its signed
    volume came out at +2398.9 m3 where the closed body holds ~13922 m3, which is what
    made the earlier gold-regression numbers meaningless (documented in
    ``docs/gap_cad_t1_iges_pipe.md``).  The flags are recorded in the returned info dict
    so the evidence file carries the orientation that was actually used.
    """
    xs, grid, spans = girth_grid(triangles, stations, girth_samples)
    hull = interpolate_surface(grid, degree=(degree, degree))
    knots_u, knots_v = hull.knots
    net = hull.control_net
    # Declared deck: ruled between the +y and -y girth end points along u, using the same
    # u-knots as the hull so the u=0 / u=1 edges are exactly the hull corner curves.
    edge_plus = [tuple(row[0]) for row in net]
    edge_minus = [tuple(row[-1]) for row in net]
    deck = NURBSSurface(
        degree=(degree, 1),
        control_net=[[a, b] for a, b in zip(edge_plus, edge_minus)],
        weights=[[1.0, 1.0] for _ in edge_plus],
        knots=[list(knots_u), [0.0, 0.0, 1.0, 1.0]],
    )
    # End caps: ruled between the hull's end edge (a NURBS curve in v) and the straight
    # deck edge at that station, parametrised on the same v-knots.
    def cap(row: list[tuple[float, float, float]]) -> NURBSSurface:
        first, last = np.asarray(row[0]), np.asarray(row[-1])
        line = [tuple(first + (last - first) * (j / (len(row) - 1))) for j in range(len(row))]
        return _ruled(list(row), line, knots_v, degree)

    names = ["hull", "deck", "cap_aft", "cap_fore"]
    surfaces = [hull, deck, cap(net[0]), cap(net[-1])]
    orientation: dict = {}
    if orient:
        provisional = sew(surfaces, names=names, source="DTMB5415_Hull_exact.stl loft")
        flags, orientation = enclosed_volume_flags(provisional)
    else:
        flags = [True] * len(surfaces)
    brep = sew(
        surfaces,
        names=names,
        same_sense=flags,
        source="DTMB5415_Hull_exact.stl loft",
    )
    info = {
        "stations": stations,
        "girth_samples": girth_samples,
        "degree": degree,
        "x_range": [float(xs[0]), float(xs[-1])],
        "girth_span_min": float(min(spans)),
        "faces": len(brep.faces),
        "hull_control_points": [len(net), len(net[0])],
        "same_sense": flags,
        "orientation": orientation,
    }
    return brep, info


# --------------------------------------------------------------------------- #
# exact patch derivatives, with the derivative control nets cached
# --------------------------------------------------------------------------- #
class PatchDerivatives:
    """``S``, ``S_u`` and ``S_v`` of a ``NURBSSurface`` with its nets built once.

    ``NURBSSurface.derivatives`` rebuilds the homogeneous (derivative) control nets on
    every call, which costs ~6 ms per point -- enough to make a converged reference
    integral unaffordable.  This class performs the *same* exact computation (A3.6
    rational quotient recursion on the same derivative nets, built with
    ``tensorcad.nurbs._differentiate_axis``) but builds the nets in ``__init__``.
    ``tests/test_iges_gold_regression.py`` asserts equality with the library to 1e-12,
    so this is a cache, not a second implementation.
    """

    def __init__(self, surface: NURBSSurface) -> None:
        self.surface = surface
        degree_u, degree_v = surface.degree
        homogeneous = tuple(
            tuple(
                (
                    surface.weights[i][j] * surface.control_net[i][j][0],
                    surface.weights[i][j] * surface.control_net[i][j][1],
                    surface.weights[i][j] * surface.control_net[i][j][2],
                    surface.weights[i][j],
                )
                for j in range(len(surface.control_net[i]))
            )
            for i in range(len(surface.control_net))
        )
        knots_u = tuple(surface.knots[0])
        knots_v = tuple(surface.knots[1])
        self._u0 = (homogeneous, degree_u, knots_u)
        reduced, degree, knots = _differentiate_axis(homogeneous, 0, degree_u, knots_u)
        self._u1 = (reduced, degree, knots)
        reduced, degree, knots = _differentiate_axis(homogeneous, 1, degree_v, knots_v)
        self._v1 = (reduced, degree, knots)
        self._v0 = (homogeneous, degree_v, knots_v)

    def __call__(self, u: float, v: float) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        (net, degree_u, knots_u), (net_u, degree_u1, knots_u1) = self._u0, self._u1
        (_, degree_v, knots_v), (net_v, degree_v1, knots_v1) = self._v0, self._v1
        h00 = tensor_point(net, (degree_u, degree_v), (knots_u, knots_v), u, v)
        h10 = tensor_point(net_u, (degree_u1, degree_v), (knots_u1, knots_v), u, v)
        h01 = tensor_point(net_v, (degree_u, degree_v1), (knots_u, knots_v1), u, v)
        w0 = h00[3]
        if w0 == 0.0:
            raise ValueError(f"degenerate evaluation at (u, v)=({u!r}, {v!r})")
        point = np.asarray(h00[:3], dtype=float) / w0
        du = (np.asarray(h10[:3], dtype=float) - h10[3] * point) / w0
        dv = (np.asarray(h01[:3], dtype=float) - h01[3] * point) / w0
        return point, du, dv


# --------------------------------------------------------------------------- #
# orientation: outward face normals without trusting a caller-supplied flag
# --------------------------------------------------------------------------- #
def _grid_points(surface: NURBSSurface, nu: int, nv: int) -> np.ndarray:
    (u0, u1), (v0, v1) = surface.domain
    return np.array([
        [np.asarray(surface.evaluate(u, v)) for v in np.linspace(v0, v1, nv)]
        for u in np.linspace(u0, u1, nu)
    ])


def enclosed_volume_flags(brep: BRep, *, order: int = 5, panels: int = 4) -> tuple[list[bool], dict]:
    """``same_sense`` per face from the enclosed volume: keep the sign set that maximises it.

    A consistently oriented closed patch set has one outward orientation and one inward
    one; any single reversed face subtracts its own (double) contribution, so the outward
    set is the one whose ``|V| = 1/3 |sum_faces s_f int S.(S_u x S_v)|`` is largest.  The
    argument is not a theorem for a general surface, so the rule is *validated* instead:
    ``tests/test_iges_gold_regression.py`` checks that it reproduces the hand-set flags of
    the analytic baselines (Series 60 Cb=0.60/0.65, KVLCC2), whose orientation is known
    from their construction, before it is used for the DTMB5415 loft.

    This replaces the earlier hard-coded ``same_sense=[True, True, True, True]``, which was
    not consistently oriented: the loft's signed volume came out as +2398.9 m3 against the
    ~13922 m3 the closed body actually encloses.
    """
    count = len(brep.faces)
    best: tuple[float, list[bool]] | None = None
    table: list[dict] = []
    for mask in range(2 ** count):
        flags = [bool(mask >> index & 1) for index in range(count)]
        volume = brep_quadrature(brep, order=order, panels=panels, flags=flags)["volume"]
        table.append({"flags": flags, "volume": volume})
        if best is None or abs(volume) > abs(best[0]) + 1e-12:
            best = (volume, flags)
    assert best is not None
    return best[1], {
        "rule": "maximise |V| = |1/3 sum_faces s_f int S.(S_u x S_v)| over the 2**n sign sets",
        "volume": best[0],
        "combinations": len(table),
        "runner_up": sorted((abs(row["volume"]) for row in table), reverse=True)[1],
    }


# --------------------------------------------------------------------------- #
# independent reference integrals (graded composite Gauss on the patches)
# --------------------------------------------------------------------------- #
_GAUSS = {
    1: ([0.0], [2.0]),
    3: ([-0.7745966692414834, 0.0, 0.7745966692414834],
        [0.5555555555555556, 0.8888888888888888, 0.5555555555555556]),
    5: ([-0.9061798459386640, -0.5384693101056831, 0.0, 0.5384693101056831, 0.9061798459386640],
        [0.2369268850561891, 0.4786286704993665, 0.5688888888888889, 0.4786286704993665,
         0.2369268850561891]),
    8: ([-0.9602898564975363, -0.7966664774136267, -0.5255324099163290, -0.1834346424956498,
         0.1834346424956498, 0.5255324099163290, 0.7966664774136267, 0.9602898564975363],
        [0.1012285362903763, 0.2223810344533745, 0.3137066458778873, 0.3626837833783620,
         0.3626837833783620, 0.3137066458778873, 0.2223810344533745, 0.1012285362903763]),
}


def graded_panels(low: float, high: float, count: int, ratio: float = 0.5) -> list[tuple[float, float]]:
    """Break ``[low, high]`` into ``count`` panels clustered geometrically at both ends.

    The analytic hull baselines carry an algebraic singularity at the keel
    (``y ~ zeta**(1/m)`` with ``m`` up to 24), so a uniform panel grid -- or a single
    panel -- under-integrates by whole per-cent: that was the defect in the first
    version of this reference.  Geometric clustering towards *both* parameter ends makes
    the composite rule converge to the exact value regardless of which end is singular.
    """
    if count <= 1:
        return [(low, high)]
    half = count // 2
    # ``ratio < 1`` must make the panels *small* at the ends and coarse in the middle:
    # the singular end of a hull patch is an endpoint, so that is where the resolution has
    # to go (``weights`` is therefore reversed relative to the panel order).
    weights = [ratio ** (half - 1 - i) for i in range(half)]
    total = sum(weights)
    fractions = [0.0]
    accumulated = 0.0
    for weight in weights:
        accumulated += 0.5 * weight / total
        fractions.append(accumulated)
    mirrored = [1.0 - fraction for fraction in reversed(fractions)]
    if count % 2:
        # odd count: use the even subdivision of ``count - 1`` and split the widest
        # (central) panel, so the panel count is exact for every count.
        even = graded_panels(low, high, count - 1, ratio)
        widest = max(range(len(even)), key=lambda i: even[i][1] - even[i][0])
        a, b = even[widest]
        middle = 0.5 * (a + b)
        return even[:widest] + [(a, middle), (middle, b)] + even[widest + 1:]
    breakpoints = fractions + mirrored[1:]
    breakpoints = sorted({round(fraction, 15) for fraction in breakpoints})
    breakpoints[0], breakpoints[-1] = 0.0, 1.0
    return [(low + (high - low) * a, low + (high - low) * b)
            for a, b in zip(breakpoints[:-1], breakpoints[1:])]


def brep_quadrature(brep: BRep, *, order: int = 5, panels: int = 4,
                    ratio: float = 0.5, flags: Sequence[bool] | None = None) -> dict:
    """Area and signed volume of a full-domain patch set by graded composite Gauss.

    ``V = 1/3 * sum_faces s_f * int int S . (S_u x S_v) du dv`` and
    ``A = sum_faces int int |S_u x S_v| du dv``, with ``s_f = +1/-1`` from
    ``face.same_sense`` (or from ``flags``).  Deliberately independent of the
    tessellation/trim-inversion path under test: it integrates the *same* NURBS patches
    that were written to the file.

    Reference quality is measured, not assumed: ``panels`` and ``order`` are the
    convergence ladder rungs reported by the benchmark driver, and the harness checks the
    result against the analytic closed forms of the baselines (Series 60, KVLCC2) where
    those exist.  A reference that does not reproduce an analytic volume to well under
    1 % is not allowed to judge the pipeline.
    """
    nodes, weights = _GAUSS[order]
    area = 0.0
    volume = 0.0
    per_face: dict[str, dict] = {}
    for index, face in enumerate(brep.faces):
        surface = face.surface
        derivative = PatchDerivatives(surface)
        (u0, u1), (v0, v1) = surface.domain
        sign = 1.0 if (flags[index] if flags is not None else face.same_sense) else -1.0
        face_area = 0.0
        face_volume = 0.0
        for ua, ub in graded_panels(u0, u1, panels, ratio):
            half_u, mid_u = 0.5 * (ub - ua), 0.5 * (ua + ub)
            for va, vb in graded_panels(v0, v1, panels, ratio):
                half_v, mid_v = 0.5 * (vb - va), 0.5 * (va + vb)
                for i, node_u in enumerate(nodes):
                    u = mid_u + half_u * node_u
                    for j, node_v in enumerate(nodes):
                        v = mid_v + half_v * node_v
                        point, du, dv = derivative(u, v)
                        cross = np.cross(du, dv)
                        weight = half_u * half_v * weights[i] * weights[j]
                        face_area += float(np.linalg.norm(cross)) * weight
                        face_volume += sign * float(np.dot(point, cross)) / 3.0 * weight
        area += face_area
        volume += face_volume
        per_face[face.name] = {"area": face_area, "volume": face_volume,
                               "same_sense": bool(flags[index] if flags is not None
                                                  else face.same_sense)}
    return {
        "volume": volume,
        "area": area,
        "order": order,
        "panels": panels,
        "grading_ratio": ratio,
        "per_face": per_face,
        "method": f"{order}x{order} Gauss per panel, {panels}x{panels} panels graded "
                  f"(ratio {ratio}) toward both parameter ends",
    }
