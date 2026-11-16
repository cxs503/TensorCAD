"""Randomised differential check: TriangleGrid vs brute force."""
import math, random, sys
sys.path.insert(0, '/root/tp-agent/work/m22896')
from tensorcad.geometry import Mesh
from tensorcad import spatial_index as si


def rand_mesh(rng, n_tris, scale=1.0):
    verts, faces = [], []
    for _ in range(n_tris):
        cx, cy, cz = (rng.uniform(-1, 1) * scale for _ in range(3))
        r = scale * 0.3
        vs = [(cx + rng.uniform(-r, r), cy + rng.uniform(-r, r), cz + rng.uniform(-r, r)) for _ in range(3)]
        i0 = len(verts)
        verts.extend(vs)
        faces.append((i0, i0 + 1, i0 + 2))
    return Mesh(vertices=verts, faces=faces)


def main(seed=0, n_mesh=40, n_rays=60):
    rng = random.Random(seed)
    bad = 0
    for m in range(n_mesh):
        n_tris = rng.choice([1, 2, 5, 20, 100, 400])
        scale = rng.choice([1e-3, 1.0, 1e3, 1e6])
        mesh = rand_mesh(rng, n_tris, scale)
        try:
            grid = si.TriangleGrid(mesh)
        except Exception as e:
            print(f"BUILD-FAIL mesh{m} n={n_tris} scale={scale}: {type(e).__name__}: {e}")
            bad += 1
            continue
        for r in range(n_rays):
            kind = rng.random()
            if kind < 0.4:
                origin = tuple(rng.uniform(-2, 2) * scale for _ in range(3))
            elif kind < 0.7:
                origin = tuple(grid._lo[i] + rng.random() * grid._ext[i] for i in range(3))
            else:
                origin = tuple(rng.uniform(-0.5, 0.5) * scale for _ in range(3))
            if rng.random() < 0.25:
                axis = rng.randrange(3)
                sign = rng.choice([-1.0, 1.0])
                d = [0.0, 0.0, 0.0]
                d[axis] = sign
            else:
                d = [rng.gauss(0, 1) for _ in range(3)]
                if d == [0.0, 0.0, 0.0]:
                    d = [1.0, 0.0, 0.0]
            try:
                got = grid.cast_ray(origin, tuple(d))
            except Exception as e:
                print(f"CAST-FAIL mesh{m} n={n_tris} scale={scale} origin={origin} d={d}: {type(e).__name__}: {e}")
                bad += 1
                continue
            ref = si.brute_force_nearest(mesh, origin, tuple(d))
            g = None if got is None else got.distance
            if (g is None) != (ref is None) or (g is not None and abs(g - ref) > 1e-9 * max(1.0, abs(ref))):
                print(f"MISMATCH cast mesh{m} n={n_tris} scale={scale} res={grid.resolution} "
                      f"oversized={grid.n_oversized} origin={origin} d={d} grid={g} ref={ref}")
                bad += 1
            all_g = grid.cast_all(origin, tuple(d))
            all_ref = si.brute_force_all(mesh, origin, tuple(d))
            refd = sorted(round(t, 9) for t in all_ref)
            if len(all_g) != len(refd):
                print(f"MISMATCH cast_all COUNT mesh{m} n={n_tris} origin={origin} d={d} grid={len(all_g)} ref={len(refd)}")
                bad += 1
            hs = grid.nearest_surface(origin)
            rd = si.brute_force_surface_distance(mesh, origin)
            if (hs is None) != (rd == math.inf) or (hs is not None and abs(hs.distance - rd) > 1e-9 * max(1.0, rd)):
                print(f"MISMATCH nearest mesh{m} n={n_tris} scale={scale} res={grid.resolution} "
                      f"origin={origin} grid={None if hs is None else hs.distance} ref={rd}")
                bad += 1
    print("bad =", bad)
    return bad


if __name__ == '__main__':
    sys.exit(1 if main(int(sys.argv[1]) if len(sys.argv) > 1 else 0) else 0)