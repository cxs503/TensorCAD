"""Structured/adversarial differential checks: TriangleGrid vs brute force."""
import math, random, sys
sys.path.insert(0, '/srv/tensor/agent/work/m22896')
from tensorcad.geometry import Mesh
from tensorcad import spatial_index as si

bad = 0


def report(tag, msg):
    global bad
    bad += 1
    print("MISMATCH", tag, msg)


def check_pair(mesh, grid, origin, d, tag):
    got = grid.cast_ray(origin, d)
    ref = si.brute_force_nearest(mesh, origin, d)
    g = None if got is None else got.distance
    if (g is None) != (ref is None) or (g is not None and abs(g - ref) > 1e-9 * max(1.0, abs(ref))):
        report(tag, f"cast o={origin} d={d} grid={g} ref={ref} res={grid.resolution} over={grid.n_oversized}")
    all_g = grid.cast_all(origin, d)
    all_r = si.brute_force_all(mesh, origin, d)
    if len(all_g) != len(sorted(round(t, 9) for t in all_r)):
        report(tag, f"cast_all count o={origin} d={d} grid={len(all_g)} ref={len(all_r)}")
    hs = grid.nearest_surface(origin)
    rd = si.brute_force_surface_distance(mesh, origin)
    if (hs is None) != (rd == math.inf) or (hs is not None and abs(hs.distance - rd) > 1e-9 * max(1.0, rd)):
        report(tag, f"nearest o={origin} grid={None if hs is None else hs.distance} ref={rd}")


def box_mesh(lo=(0.0, 0.0, 0.0), hi=(1.0, 1.0, 1.0)):
    x0, y0, z0 = lo
    x1, y1, z1 = hi
    v = [(x0, y0, z0), (x1, y0, z0), (x1, y1, z0), (x0, y1, z0),
         (x0, y0, z1), (x1, y0, z1), (x1, y1, z1), (x0, y1, z1)]
    q = [(0, 2, 1), (0, 3, 2), (4, 5, 6), (4, 6, 7),
         (0, 1, 5), (0, 5, 4), (1, 2, 6), (1, 6, 5),
         (2, 3, 7), (2, 7, 6), (3, 0, 4), (3, 4, 7)]
    return Mesh(vertices=v, faces=q)


def sheet_mesh(n=6, extents=(1.0, 1.0, 1e-9)):
    """Two triangles spanning a very thin slab: extreme aspect ratio."""
    sx, sy, sz = extents
    v = [(0, 0, 0), (sx, 0, 0), (sx, sy, 0), (0, sy, 0),
         (0, 0, sz), (sx, 0, sz), (sx, sy, sz), (0, sy, sz)]
    q = [(0, 1, 2), (0, 2, 3), (4, 6, 5), (4, 7, 6)]
    return Mesh(vertices=v, faces=q)


def degenerate_mesh():
    v = [(0, 0, 0), (1, 0, 0), (2, 0, 0),      # collinear
         (0, 1, 0), (0, 1, 0), (0, 1, 0),      # all-equal
         (0, 0, 1), (1, 0, 1), (0, 1, 1)]
    q = [(0, 1, 2), (3, 4, 5), (6, 7, 8)]
    return Mesh(vertices=v, faces=q)


# 1. exact box hits: rays along axes, through vertices/edges, from inside
mesh = box_mesh()
for res in (None, (1, 1, 1), (2, 2, 2), (7, 3, 5), (64, 64, 64)):
    grid = si.TriangleGrid(mesh) if res is None else si.TriangleGrid(mesh, resolution=res)
    for origin in [(0.5, 0.5, 0.5), (-1, 0.5, 0.5), (0.5, -1, 0.5), (-1, -1, -1),
                   (0.0, 0.0, 0.0), (1.0, 1.0, 1.0), (0.0, 0.5, 0.5), (1e-12, 0.5, 0.5)]:
        for d in [(1, 0, 0), (0, 1, 0), (0, 0, 1), (-1, 0, 0), (0, -1, 0), (0, 0, -1),
                  (1, 1, 0), (1, 0, 1), (1, 1, 1), (1, -1, 0)]:
            check_pair(mesh, grid, origin, d, f"box res={res}")

# 2. contains() parity on the unit box
grid = si.TriangleGrid(mesh)
rng = random.Random(7)
inside_bad = 0
for _ in range(400):
    p = tuple(rng.uniform(-0.5, 1.5) for _ in range(3))
    truth = all(0.0 < c < 1.0 for c in p)
    if grid.contains(p) != truth:
        inside_bad += 1
        if inside_bad <= 5:
            report("contains", f"p={p} grid={grid.contains(p)} truth={truth}")
print("contains mismatches:", inside_bad)

# 3. thin slab / extreme aspect ratios
for ext in [(1.0, 1.0, 1e-9), (1e-6, 1.0, 1.0), (1.0, 1e-7, 1.0), (1e3, 1e-3, 1.0)]:
    m = sheet_mesh(extents=ext)
    g = si.TriangleGrid(m)
    for origin in [(0.5, 0.5, -1.0), (0.5, 0.5, 0.5), (-1, 0.5, 0.5), (0.5, 0.5, 1e-10)]:
        for d in [(0, 0, 1), (1, 0, 0), (0, 1, 0), (1, 1, 1)]:
            check_pair(m, g, origin, d, f"slab ext={ext}")
    print("slab", ext, "res", g.resolution, "over", g.n_oversized)

# 4. fully degenerate triangles
m = degenerate_mesh()
for res in (None, (2, 2, 2)):
    g = si.TriangleGrid(m) if res is None else si.TriangleGrid(m, resolution=res)
    for origin in [(0, 0, -1), (0.5, 0.1, 0.5), (0, 0, 0), (3, 3, 3)]:
        for d in [(0, 0, 1), (1, 0, 0), (0, 0, -1), (-1, -1, -1)]:
            check_pair(m, g, origin, d, f"degenerate res={res}")

# 5. random rays grazing exact edges of a box grid surface
hi_res_box = box_mesh(hi=(10.0, 4.0, 2.0))
g = si.TriangleGrid(hi_res_box, resolution=(10, 4, 2))
for _ in range(300):
    axis = rng.randrange(3)
    o = list(rng.uniform(-2, 12) for _ in range(3))
    d = [0.0, 0.0, 0.0]
    d[axis] = rng.choice([-1.0, 1.0])
    check_pair(hi_res_box, g, tuple(o), tuple(d), "grazing")

# 6. nearest_surface with band -------------------------------------------------
m = box_mesh(hi=(2.0, 1.0, 1.0))
g = si.TriangleGrid(m)
for band in (0.0, 0.05, 0.5, 10.0):
    for _ in range(200):
        p = tuple(rng.uniform(-2, 3) for _ in range(3))
        hs = g.nearest_surface(p, band=band)
        rd = si.brute_force_surface_distance(m, p)
        if band > 0 and rd > band:
            if hs is not None:
                report("band", f"p={p} band={band} ref={rd} grid={hs.distance} (expected None)")
            continue
        if hs is None:
            report("band", f"p={p} band={band} ref={rd} grid=None")
        elif abs(hs.distance - rd) > 1e-9 * max(1.0, rd):
            report("band", f"p={p} band={band} ref={rd} grid={hs.distance}")

print("bad =", bad)
sys.exit(1 if bad else 0)