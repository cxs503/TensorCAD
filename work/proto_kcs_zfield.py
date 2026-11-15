"""Prototype v2: height-field (z-buffer) integration of the official KCS half-hull IGES."""
import math, sys, time
from tensorcad.iges import read_iges

PATH = '/root/sim-assets/m21253-379b12515423/artifacts/baseline_fleet/raw/closure/kcs/KCS_hull_SVA.igs'
T = 10.8


def tess(f, nu, nv):
    tris = []
    for s in f.surfaces():
        g = s.sample(nu, nv)
        for j in range(nv - 1):
            for i in range(nu - 1):
                a, b = g[j*nu+i], g[j*nu+i+1]
                c, d = g[(j+1)*nu+i+1], g[(j+1)*nu+i]
                tris.append((a, b, c)); tris.append((a, c, d))
    return tris


def zfield(tris, nx, ny, xlo, xhi, yhi):
    """Lowest surface height per (x,y) cell centre; returns dict cell->(zmin, tri)."""
    dx = (xhi - xlo) / nx; dy = yhi / ny
    cell = {}
    for tri in tris:
        xs = [p[0] for p in tri]; ys = [p[1] for p in tri]; zs = [p[2] for p in tri]
        i0 = max(0, int((min(xs)-xlo)/dx)); i1 = min(nx-1, int((max(xs)-xlo)/dx)+1)
        j0 = max(0, int(min(ys)/dy)); j1 = min(ny-1, int(max(ys)/dy)+1)
        if i1 < i0 or j1 < j0:
            continue
        x1, y1, z1 = tri[0]; x2, y2, z2 = tri[1]; x3, y3, z3 = tri[2]
        det = (y2-y3)*(x1-x3) + (x3-x2)*(y1-y3)
        if abs(det) < 1e-14:
            continue
        for i in range(i0, i1+1):
            cx = xlo + (i+0.5)*dx
            for j in range(j0, j1+1):
                cy = (j+0.5)*dy
                l1 = ((y2-y3)*(cx-x3) + (x3-x2)*(cy-y3)) / det
                l2 = ((y3-y1)*(cx-x3) + (x1-x3)*(cy-y3)) / det
                l3 = 1.0 - l1 - l2
                if l1 < -1e-9 or l2 < -1e-9 or l3 < -1e-9:
                    continue
                z = l1*z1 + l2*z2 + l3*z3
                k = (i, j)
                if k not in cell or z < cell[k]:
                    cell[k] = z
    return cell, dx, dy


def volume(cell, dx, dy, draft=T):
    return sum(max(draft - z, 0.0) * dx * dy for z in cell.values())


if __name__ == '__main__':
    f = read_iges(PATH)
    for grid in (12, 20, 32):
        t0 = time.time()
        tris = tess(f, grid, grid)
        for raster in (0.5, 0.25):
            nx = int(round(250.0/raster)); ny = int(round(18.0/raster))
            cell, dx, dy = zfield(tris, nx, ny, -8.0, 240.0, 18.0)
            V = 2.0*volume(cell, dx, dy)
            Awp = 2.0*len([1 for z in cell.values() if z < T])*dx*dy
            print(f"grid={grid}^2 tris={len(tris):6d} raster={raster} nx*ny={nx}x{ny} "
                  f"V={V:10.1f} (pub 52030: {100*(V-52030)/52030:+.3f}%) Awp={Awp:8.1f} "
                  f"t={time.time()-t0:.1f}s")