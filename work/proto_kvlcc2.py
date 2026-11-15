"""Prototype: KVLCC2 official IGES -> trim-aware tessellation -> height-field displacement."""
import math, sys, time
from tensorcad.iges import read_iges
from tensorcad.fleet.trimfield import trim_polygons, point_in_polygon

PATH = '/root/sim-assets/m20163-fedab06121c6/artifacts/baseline_fleet/raw/KVLCC2/kvlcc2.igs'
T = 20.8
PUB = 312622.0
GRID = int(sys.argv[1]) if len(sys.argv) > 1 else 60
RASTER = float(sys.argv[2]) if len(sys.argv) > 2 else 1.0


def trimmed_tris(f, grid):
    records = trim_polygons(f)
    tris = []
    stats = {"decided": 0, "undecided": 0, "quads_kept": 0, "quads_dropped": 0}
    for rec in records:
        basis, loops = rec.get("basis"), rec.get("loops") or []
        if basis is None or not loops:
            stats["undecided"] += 1
            continue
        stats["decided"] += 1
        u0, u1 = basis.domain_u; v0, v1 = basis.domain_v
        pts = [[basis.evaluate(u0 + (u1-u0)*i/grid, v0 + (v1-v0)*j/grid) for i in range(grid+1)]
               for j in range(grid+1)]
        for j in range(grid):
            for i in range(grid):
                uc = u0 + (u1-u0)*(i+0.5)/grid; vc = v0 + (v1-v0)*(j+0.5)/grid
                if sum(1 for lp in loops if point_in_polygon(uc, vc, lp)) % 2 == 0:
                    stats["quads_dropped"] += 1; continue
                stats["quads_kept"] += 1
                a, b = pts[j][i], pts[j][i+1]
                c, d = pts[j+1][i+1], pts[j+1][i]
                tris.append((a, b, c)); tris.append((a, c, d))
    return tris, stats


def zfield(tris, nx, ny, xlo, xhi, yhi):
    dx = (xhi-xlo)/nx; dy = yhi/ny; cell = {}
    for tri in tris:
        xs = [p[0] for p in tri]; ys = [p[1] for p in tri]
        i0 = max(0, int((min(xs)-xlo)/dx)); i1 = min(nx-1, int((max(xs)-xlo)/dx)+1)
        j0 = max(0, int(min(ys)/dy)); j1 = min(ny-1, int(max(ys)/dy)+1)
        (x1, y1, z1), (x2, y2, z2), (x3, y3, z3) = tri
        det = (y2-y3)*(x1-x3) + (x3-x2)*(y1-y3)
        if abs(det) < 1e-14: continue
        for i in range(i0, i1+1):
            cx = xlo + (i+0.5)*dx
            for j in range(j0, j1+1):
                cy = (j+0.5)*dy
                l1 = ((y2-y3)*(cx-x3) + (x3-x2)*(cy-y3))/det
                l2 = ((y3-y1)*(cx-x3) + (x1-x3)*(cy-y3))/det
                l3 = 1.0-l1-l2
                if l1 < -1e-9 or l2 < -1e-9 or l3 < -1e-9: continue
                z = l1*z1+l2*z2+l3*z3
                k = (i, j)
                if k not in cell or z < cell[k]: cell[k] = z
    return cell, dx, dy


f = read_iges(PATH)
t0 = time.time()
tris, st = trimmed_tris(f, GRID)
# flip frame: waterline z=0, hull below -> z' = -z
# frame: waterline at z=0, hull below (measured: patch z-ranges bottom out at -20.80 = draft)
xs = [p[0] for t in tris for p in t]; ys = [p[1] for t in tris for p in t]; zs = [p[2] for t in tris for p in t]
print(f"trim: {st} tris={len(tris)} t={time.time()-t0:.1f}s")
print("trimmed bbox x[%.3f, %.3f] (L=%.3f) y[%.3f, %.3f] z[%.3f, %.3f]" %
      (min(xs), max(xs), max(xs)-min(xs), min(ys), max(ys), min(zs), max(zs)))
cell, dx, dy = zfield(tris, int(round((max(xs)-min(xs))/RASTER)), int(round(max(ys)/RASTER)), min(xs), max(xs), max(ys))
wet = {k: z for k, z in cell.items() if z < 0.0}
V = 2.0*sum(max(0.0-z, 0.0) for z in wet.values())*dx*dy
Awp = 2.0*len(wet)*dx*dy
print(f"raster={RASTER} cols={len(cell)} wet={len(wet)} V={V:,.1f} vs pub {PUB:,.0f}: {100*(V-PUB)/PUB:+.3f}%")
print(f"Awp={Awp:,.1f}  Cb={V/(320.0*58.0*20.8):.4f}")