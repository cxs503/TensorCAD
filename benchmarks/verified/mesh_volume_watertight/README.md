# Verified case `mesh_volume_watertight`

Anchor **`CASE-CAD-VOLUME-CLOSED`** (case pool, round 1): does a closed TriangleCAD mesh
bound the volume of the body it claims to be, and is it actually closed?

**Reference = analytic closed forms, evaluated in code (no transcribed decimals).**
The full report with every number injected from `result.json` is
[`docs/mesh_volume_watertight_bench.md`](../../../docs/mesh_volume_watertight_bench.md).

| quantity | reference | where it comes from |
|---|---|---|
| box volume / area | `a b c` / `2(ab+bc+ca)` | elementary solid geometry, `tensorcad.mesh_volume` |
| sphere volume | `4/3 pi r^3`, cross-checked by `(2/3)(pi r^2)(2r)` | Euclid *Elements* XII.18 (sphere = 2/3 of the circumscribed cylinder) |
| sphere area | `4 pi r^2` | same |
| icosahedron volume | `(5/12)(3+sqrt5) a^3`, `a = 4R/sqrt(10+2 sqrt5)` | regular-polyhedron closed form; this is the `sub0` mesh itself |
| closed-body topology | `open_edges = 0`, `euler = 2`, no non-manifold edge / duplicate facet / winding conflict, `signed_volume > 0`, audit `passed` | `tensorcad.geometry_audit.audit_watertightness` |

## Gates (predeclared, in `tensorcad/mesh_volume.py`)

| gate | value | why |
|---|---|---|
| `EXACT_TOLERANCE` | `1e-12` relative | a triangle mesh represents a box exactly; only round-off may appear |
| `RELATIVE_TOLERANCE` | `1e-3` (0.1 %) | the case's acceptance line, applied to the finest ladder rung (sub6) and to the r=2.5 m scale case |
| `ROUTE_TOLERANCE` | `1e-10` | divergence theorem vs centroid-tet summation over the same triangles |
| `FLOAT32_TOLERANCE` | `1e-6` | binary STL stores float32; the file, not the integrator, is the error source |
| observed order band | `[1.5, 2.5]` | piecewise-linear approximation of a smooth surface is second order |
| **falsification gate** | sub4 **must** miss 0.1 % | if a coarse inscribed polyhedron "gave" the sphere volume, the measurement would be detached from the geometry |

## Run

```bash
python benchmarks/verified/mesh_volume_watertight/run_benchmark.py            # ~20 s, exits 1 on a broken gate
python benchmarks/verified/mesh_volume_watertight/run_benchmark.py --out runs/mesh_volume_watertight
python -m pytest tests/test_mesh_volume.py -q                                  # 12 tests, includes one full run
python benchmarks/verified/mesh_volume_watertight/publish_report.py           # regenerate the doc + report payload
```

Artifacts written here: `result.json` (all gates, bodies, ladder, controls, environment),
`convergence.csv`, `convergence.png` (log-log error vs h with a slope-2 reference and the
0.1 % line), `gates.png` (every gate's |relative error|), `box_stl.stl` (float32 round-trip
sample).

## What the case established

* Box family: volume and area **exact** (0.0 absolute error; `8.327e-15` relative after a
  rotation + translation), i.e. the divergence-theorem integrator is verified on a body where
  the answer is a closed form.
* Icosphere: the *inscribed* tessellation is **below** the sphere and converges at the mesh's
  own order — observed `p = 1.989` (volume) / `2.005` (area) over sub0…sub6, crossing the
  0.1 % line between sub4 (−0.216 %) and sub5 (−0.054 %), reaching `-1.35e-4` at sub6.
* Watertightness: four injected defect families (hole, flipped facets, duplicate facets,
  unwelded triangle soup) are all named, and the audit stays green on every clean rung.
* The two checks are **not interchangeable**: an unwelded soup has *zero* volume error yet 960
  open edges; conversely a single duplicated or flipped facet at sub6 moves the volume by only
  `+1.15e-5` / `-2.31e-5` (inside the 0.1 % gate) while the audit still names it.