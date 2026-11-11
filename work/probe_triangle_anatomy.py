"""Probe 3: anatomy of the worst SUBOFF triangles (edges, direction, location)."""
from __future__ import annotations

import math
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from tensorcad.parametric_hull import (  # noqa: E402
    loft_nurbs_surface,
    suboff_cross_sections,
    suboff_station_fractions,
    tessellate_surface,
)

stations = suboff_station_fractions(41)
sections = suboff_cross_sections(n_points=32, station_fractions=stations)
surface = loft_nurbs_surface(sections)
m = tessellate_surface(surface, n_u=160, n_v=16)
mesh = m.mesh
faces = mesh.faces[: m.hull_face_count]

def describe(fidx):
    f = faces[fidx]
    pts = [mesh.vertices[i] for i in f]
    e = []
    for a, b in ((0, 1), (1, 2), (2, 0)):
        e.append((math.dist(pts[a], pts[b]), a, b))
    e.sort()
    print(" face %d: %s" % (fidx, [(round(p[0], 4), round(p[1], 4), round(p[2], 4)) for p in pts]))
    for ln, a, b in e:
        pa, pb = pts[a], pts[b]
        dx = abs(pa[0] - pb[0])
        circ = math.hypot(pa[1] - pb[1], pa[2] - pb[2])
        print("   edge len=%.5f  |dx|=%.5f  circumferential=%.5f  %s" %
              (ln, dx, circ, "MERIDIONAL" if dx > circ else "CIRCUMFERENTIAL"))
    return e

# worst faces found earlier
for fidx in (9464, 9460):
    describe(fidx)

# bow region worst
best = None
for idx, f in enumerate(faces):
    pts = [mesh.vertices[i] for i in f]
    e = sorted(math.dist(pts[a], pts[b]) for a, b in ((0, 1), (1, 2), (2, 0)))
    if e[0] <= 1e-15:
        continue
    ux, uy, uz = (pts[1][i] - pts[0][i] for i in range(3))
    vx, vy, vz = (pts[2][i] - pts[0][i] for i in range(3))
    nx, ny, nz = (uy * vz - uz * vy, uz * vx - ux * vz, ux * vy - uy * vx)
    twoA = math.sqrt(nx * nx + ny * ny + nz * nz)
    jac = twoA / (e[1] * e[2])
    xc = sum(p[0] for p in pts) / 3
    if xc < 0.5 and (best is None or jac < best[0]):
        best = (jac, idx, xc)
print("\nworst bow face: jac=%.4f idx=%d x=%.4f" % best)
describe(best[1])

# grid ring radii near the ends (to see the actual cap ring radii)
grid = [[surface.evaluate(i / 159, j / 15) for j in range(16)] for i in range(160)]
print("\nnear-stern rings (x, max half-breadth, min edge in ring):")
for i in range(150, 160):
    row = grid[i]
    r = max(p[1] for p in row)
    arcs = [math.dist(row[j], row[j + 1]) for j in range(15)]
    print("  i=%3d x=%.4f r=%.4f arc_min=%.5f arc_max=%.5f dx=%.5f" %
          (i, row[0][0], r, min(arcs), max(arcs), abs(grid[i + 1][0][0] - row[0][0]) if i < 159 else 0))
print("\nnear-bow rings:")
for i in range(0, 6):
    row = grid[i]
    r = max(p[1] for p in row)
    arcs = [math.dist(row[j], row[j + 1]) for j in range(15)]
    print("  i=%3d x=%.4f r=%.4f arc_min=%.5f arc_max=%.5f dx=%.5f" %
          (i, row[0][0], r, min(arcs), max(arcs), abs(grid[i + 1][0][0] - row[0][0])))