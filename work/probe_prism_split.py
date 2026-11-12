"""Probe: can a prism-based split of the wall-layer hexes beat the Kuhn split's tet quality?

Read-only exploration; nothing here is imported by the package.
"""
from __future__ import annotations

from collections import Counter
from itertools import combinations

from tensorcad.hull import HullParameters
from tensorcad.volume_mesh import (
    _corner_index_map,
    _tet_signed_volume,
    build_volume_mesh,
    tet_scaled_jacobian,
    volume_mesh_from_hull,
)

L, B, T = 10.0, 2.0, 0.8


def tetra_prism(P, Q):
    """Three tets tiling the prism (P0,P1,P2) low -> (Q0,Q1,Q2) high."""
    p0, p1, p2 = P
    q0, q1, q2 = Q
    return [(p0, p1, p2, q0), (p1, p2, q0, q1), (p2, q0, q1, q2)]


def prism_tets(cell, along, diag=0):
    """Split a hex into two prisms extruded along index axis ``along``, then into tets."""
    m = _corner_index_map(cell)
    axes = [a for a in (0, 1, 2) if a != along]
    p, q = axes

    def node(a, b, t):
        off = [0, 0, 0]
        off[p] = a
        off[q] = b
        off[along] = t
        return m[(off[0], off[1], off[2])]

    lo = [node(0, 0, 0), node(1, 0, 0), node(1, 1, 0), node(0, 1, 0)]
    hi = [node(0, 0, 1), node(1, 0, 1), node(1, 1, 1), node(0, 1, 1)]
    if diag == 0:  # diagonal 0-2 of the base quad
        tri_lo = ((0, 1, 2), (0, 2, 3))
    else:  # diagonal 1-3
        tri_lo = ((0, 1, 3), (1, 2, 3))
    out = []
    for tri in tri_lo:
        out += tetra_prism([lo[i] for i in tri], [hi[i] for i in tri])
    return out


def sign_fix(verts, tet):
    return tet if _tet_signed_volume(verts, tet) > 0 else (tet[0], tet[2], tet[1], tet[3])


def stats(values, thr=0.1):
    vals = sorted(values)
    n = len(vals)
    return {
        "count": n,
        "min": round(vals[0], 4) if n else None,
        "p05": round(vals[int(0.05 * (n - 1))], 4) if n else None,
        "median": round(vals[n // 2], 4) if n else None,
        "below_thr_%": round(100.0 * sum(1 for v in vals if v < thr) / n, 2) if n else None,
    }


def watertight(tets):
    faces = Counter()
    for tet in tets:
        for face in combinations(tet, 3):
            faces[tuple(sorted(face))] += 1
    good = sum(1 for v in faces.values() if v == 2)
    return {"faces": len(faces), "paired": good, "unpaired": len(faces) - good,
            "max_mult": max(faces.values())}


def measure(mesh, tets, label, wall_cells, order):
    pos = [t for t in tets if abs(_tet_signed_volume(mesh.vertices, t)) > 1e-14]
    zeros = len(tets) - len(pos)
    vals = [tet_scaled_jacobian(mesh.vertices, t) for t in pos]
    wall = [v for v, t in zip(vals, pos) if all(c in wall_cells for c in t)]
    core = [v for v, t in zip(vals, pos) if all(c not in wall_cells for c in t)]
    vol = sum(abs(_tet_signed_volume(mesh.vertices, t)) for t in pos)
    print(f"-- {label}: tets={len(tets)} zero={zeros} vol={vol:.6f}")
    print("   all ", stats(vals))
    if wall:
        print("   wall", stats(wall))
    if core:
        print("   core", stats(core))
    print("   watertight", watertight(pos))
    return pos


if __name__ == "__main__":
    params = HullParameters(length=L, beam=B, draft=T, nx=41, nz=13)
    mesh = volume_mesh_from_hull(params, first_layer=0.005, growth=1.2, layers=4)
    half = 4
    bc = 2 * half
    # cell index order is i, j, k with k fastest
    wall_cells = set()
    for idx in range(len(mesh.hexes)):
        k = idx % bc
        if k < half - 1 or k >= bc - half + 1:
            wall_cells.add(idx)
    print("hexes", len(mesh.hexes), "volume", round(mesh.volume(), 6), "tets", len(mesh.tets))
    measure(mesh, list(mesh.tets), "Kuhn (shipped)", wall_cells, None)
    for along in (1, 0, 2):
        for diag in (0, 1):
            tets = []
            for cell in mesh.hexes:
                tets += [sign_fix(mesh.vertices, t) for t in prism_tets(cell, along, diag)]
            measure(mesh, tets, f"prism along={along} diag={diag}", wall_cells, None)