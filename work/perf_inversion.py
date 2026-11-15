"""Micro-benchmark: cost of 3-D trim-loop -> (u,v) inversion, per loop point.

Usage: python work/perf_inversion.py FILE [MAX_FACES] [REPEAT]

Reads the trimmed model once, then inverts the loops of the first MAX_FACES faces
and reports wall time, loop-point count, per-point cost and the worst residual.
This is the hot spot of `import_cad` on publisher IGES files.
"""
import sys
import time

sys.path.insert(0, "/srv/tensor/agent/work/m22464")
from tensorcad.geomquality.iges_trim import read_iges_trimmed  # noqa: E402
from tensorcad.geomquality.trim_mesh import _loop_to_uv, _model_extent_from_faces  # noqa: E402

path = sys.argv[1]
max_faces = int(sys.argv[2]) if len(sys.argv) > 2 else 2
repeat = int(sys.argv[3]) if len(sys.argv) > 3 else 1

t0 = time.perf_counter()
model = read_iges_trimmed(path)
t_read = time.perf_counter() - t0
extent = _model_extent_from_faces(model.faces)
tolerance = 1e-6 * max(extent, 1e-9)

selected = model.faces[:max_faces]
loops = [(face.index, face.surface, loop)
         for face in selected for loop in face.loops]

points = 0
worst = 0.0
t0 = time.perf_counter()
for _ in range(repeat):
    for _idx, surface, loop in loops:
        params, residual = _loop_to_uv(surface, loop, tolerance)
        points += len(params)
        worst = max(worst, residual)
t_inv = time.perf_counter() - t0

print(f"file={path} faces_in_model={len(model.faces)} sampled_faces={len(selected)} "
      f"loops={len(loops)} repeat={repeat}")
print(f"read_s={t_read:.2f} inversion_s={t_inv:.2f} points={points} "
      f"us_per_point={1e6 * t_inv / max(points, 1):.0f} worst_residual={worst:.3e} "
      f"tolerance={tolerance:.3e}")