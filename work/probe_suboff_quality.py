"""Probe: where does the SUBOFF parametric mesh lose its quality?"""
from __future__ import annotations

import math
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from tensorcad.parametric_hull import (  # noqa: E402
    build_parametric_hull,
    mesh_quality_report,
    suboff_cross_sections,
    suboff_station_fractions,
)

stations = suboff_station_fractions(41)
sections = suboff_cross_sections(n_points=32, station_fractions=stations)
hull = build_parametric_hull(sections, mesh_n_u=160, mesh_n_v=16, close_waterplane=False)
q = hull.quality
print("n_u,n_v =", hull.mesh_data.n_u, hull.mesh_data.n_v)
for k in ("min_scaled_jacobian", "max_aspect_ratio", "mean_aspect_ratio",
          "min_triangle_area", "n_degenerate_faces", "watertight", "edge_manifold"):
    print(f"  {k}: {q[k]}")

mesh = hull.mesh_data.mesh
faces = mesh.faces[: hull.mesh_data.hull_face_count]
rows = []
for idx, f in enumerate(faces):
    v0, v1, v2 = (mesh.vertices[i] for i in f)
    e = sorted((math.dist(v0, v1), math.dist(v1, v2), math.dist(v2, v0)))
    if e[0] <= 1e-15:
        continue
    ux, uy, uz = (v1[i] - v0[i] for i in range(3))
    vx, vy, vz = (v2[i] - v0[i] for i in range(3))
    nx, ny, nz = (uy * vz - uz * vy, uz * vx - ux * vz, ux * vy - uy * vx)
    twoA = math.sqrt(nx * nx + ny * ny + nz * nz)
    jac = twoA / (e[1] * e[2]) if e[1] * e[2] > 0 else 0.0
    xc = (v0[0] + v1[0] + v2[0]) / 3
    rows.append((jac, e[2] / e[0], xc, idx))

rows.sort()
print("\nworst 12 by scaled Jacobian (jac, aspect, x_center, face):")
for r in rows[:12]:
    print("   %.4f %8.2f x=%.4f f=%d" % r)
print("\nx distribution of the 200 worst faces:")
worst200 = rows[:200]
print("  min x %.4f max x %.4f" % (min(r[2] for r in worst200), max(r[2] for r in worst200)))
buckets = {"x<0.1": 0, "0.1-0.5": 0, "mid": 0, ">L-0.5": 0, ">L-0.1": 0}
L = 4.3561
for jac, asp, xc, idx in worst200:
    if xc < 0.1:
        buckets["x<0.1"] += 1
    elif xc < 0.5:
        buckets["0.1-0.5"] += 1
    elif xc > L - 0.1:
        buckets[">L-0.1"] += 1
    elif xc > L - 0.5:
        buckets[">L-0.5"] += 1
    else:
        buckets["mid"] += 1
print("  ", buckets)
rows.sort(key=lambda r: -r[1])
print("\nworst 6 by aspect (aspect, jac, x):")
for jac, asp, xc, idx in rows[:6]:
    print("   %.4f %8.2f x=%.4f" % (jac, asp, xc))