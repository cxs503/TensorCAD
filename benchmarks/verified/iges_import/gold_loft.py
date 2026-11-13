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
from dataclasses import dataclass
from pathlib import Path

import numpy as np

from tensorcad.geomquality.baselines import greville_abscissae
from tensorcad.geomquality.brep import BRep, Face, sew
from tensorcad.nurbs import NURBSSurface

__all__ = [
    "read_stl", "section_polyline", "girth_grid", "interpolate_surface",
    "loft_gold_hull", "brep_quadrature",
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
                   degree: int = 3) -> tuple[BRep, dict]:
    """Full hull B-rep: girth patch + ruled deck + ruled bow/stern caps."""
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

    brep = sew(
        [hull, deck, cap(net[0]), cap(net[-1])],
        names=["hull", "deck", "cap_aft", "cap_fore"],
        same_sense=[True, True, True, True],
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
    }
    return brep, info


# --------------------------------------------------------------------------- #
# independent reference integrals (Gauss quadrature on the patches)
# --------------------------------------------------------------------------- #
_GAUSS = {
    1: ([0.0], [2.0]),
    3: ([-0.7745966692414834, 0.0, 0.7745966692414834],
        [0.5555555555555556, 0.8888888888888888, 0.5555555555555556]),
    5: ([-0.9061798459386640, -0.5384693101056831, 0.0, 0.5384693101056831, 0.9061798459386640],
        [0.2369268850561891, 0.4786286704993665, 0.5688888888888889, 0.4786286704993665,
         0.2369268850561891]),
}


def brep_quadrature(brep: BRep, *, order: int = 5) -> dict[str, float]:
    """Area and signed volume of a full-domain patch set by Gauss-Legendre quadrature.

    Deliberately independent of the tessellation path: ``V = 1/3 ∮ S·(S_u × S_v) du dv``
    and ``A = ∫∫ |S_u × S_v| du dv`` with a 5x5 rule per patch.
    """
    nodes, weights = _GAUSS[order]
    area = 0.0
    volume = 0.0
    for face in brep.faces:
        surface = face.surface
        (u0, u1), (v0, v1) = surface.domain
        ju = 0.5 * (u1 - u0)
        jv = 0.5 * (v1 - v0)
        for i, nu in enumerate(nodes):
            u = 0.5 * (u0 + u1) + ju * nu
            for j, nv in enumerate(nodes):
                v = 0.5 * (v0 + v1) + jv * nv
                derivatives = surface.derivatives(u, v, order=1)
                point = np.asarray(surface.evaluate(u, v))
                du = np.asarray(derivatives[(1, 0)])
                dv = np.asarray(derivatives[(0, 1)])
                cross = np.cross(du, dv)
                weight = ju * jv * weights[i] * weights[j]
                area += float(np.linalg.norm(cross)) * weight
                signed = float(np.dot(point, cross)) / 3.0
                sign = 1.0 if face.same_sense else -1.0
                volume += sign * signed * weight
    return {"area": area, "volume": volume}