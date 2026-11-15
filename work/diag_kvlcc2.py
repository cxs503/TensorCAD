"""Diagnose KVLCC2 trim-loop vs surface: bounding boxes and best achievable residual."""
import sys

sys.path.insert(0, "/root/tp-agent/work/m22464")
from tensorcad.geomquality.iges_trim import read_iges_trimmed  # noqa: E402
from tensorcad.geomquality.trim_mesh import _grid_seeds, invert_point  # noqa: E402


def bbox(points):
    xs = [p[0] for p in points]
    ys = [p[1] for p in points]
    zs = [p[2] for p in points]
    return [(round(min(v), 3), round(max(v), 3)) for v in (xs, ys, zs)]


path = sys.argv[1] if len(sys.argv) > 1 else "work/assets/kvlcc2.igs"
model = read_iges_trimmed(path)
print("faces", len(model.faces), "units", model.units_name, "metres_per_unit", model.metres_per_unit)
for face in model.faces[:3]:
    (u0, u1), (v0, v1) = face.surface.domain
    print(f"--- face {face.index} surface={face.surface_kind} domain=({u0:.4g},{u1:.4g})x({v0:.4g},{v1:.4g})")
    grid = [face.surface.evaluate(u, v) for u, v in _grid_seeds(u0, u1, v0, v1, 8)]
    print("    surface bbox", bbox(grid))
    for li, loop in enumerate(face.loops):
        print(f"    loop {li}: {len(loop)} points bbox {bbox(loop)}")
        for pi in (0, len(loop) // 2):
            p = loop[pi]
            uv = invert_point(face.surface, p, seeds=_grid_seeds(u0, u1, v0, v1, 8),
                              iterations=40, tolerance=1e-9)
            print(f"      point[{pi}]={tuple(round(c, 3) for c in p)} best uv=({uv[0]:.4g},{uv[1]:.4g}) residual={uv[2]:.4e}")
    if face.uv_loops:
        print("    uv_loops provided:", [len(x) if x else None for x in face.uv_loops])