# Machine-gated, resolution-controlled export

Every number below was produced by `examples/hull_stl_export/run_example.py` on
this checkout (`L = 10 m`, `B = 2 m`, `T = 0.8 m`); the JSON files it writes are
the source of the tables, not a transcription.

## The two failures this closes

A geometry artifact can be wrong in ways that a viewer cannot show:

1. **The half model.** A starboard half is a surface with no inside: every edge on
   the symmetry plane is used once, the enclosed volume is meaningless, and a
   volume mesher will either refuse it or close it with an arbitrary lid. This
   repository exported such models before the topology checks existed.
2. **The under-resolved model.** A closed, manifold, outward-oriented mesh can
   still integrate the wrong displacement. The historical 304-facet export was
   topologically perfect and several percent off in volume.

`tensorcad.validation.inspect_geometry` reports both conditions, but a report is
easy to ignore. `tensorcad.gate.gate_mesh` turns the report into a gate, and
`enforce_gate` raises `MeshGateError` instead of handing the mesh on. Every
reference value it uses comes from outside the mesh, so a self-consistent wrong
export cannot satisfy it.

## Gate contract

| check | what it proves | reference |
|---|---|---|
| `finite_vertex_coordinates`, `non_empty_mesh`, `face_indices_in_range` | the file parsed into geometry at all | numeric |
| `watertight` | no boundary edge: the surface is closed | edge/face incidence |
| `edge_manifold` | no edge shared by more than two faces | edge/face incidence |
| `vertex_manifold` | every vertex link is a single cycle (no pinches) | link traversal |
| `no_degenerate_faces` | no zero-area triangle | cross product |
| `outward_orientation` | signed volume is positive: the normals point out | divergence theorem |
| `volume_within_bounding_box` | volume is physically possible for the bounds | bounding box |
| `component_count` | the requested number of shells (default 1) | union-find over faces |
| `euler_characteristic` | `V - E + F = 2·components` (closed genus-0 shell) | topology, not the mesh's own output |
| `reference_volume` (opt-in) | the integrated displacement matches an **external** reference | `tensorcad.hull.analytic_volume` |
| `facet_count` (opt-in) | the consumer's minimum resolution is met | caller's budget |
| `facet_aspect` (opt-in) | no sliver beyond `longest_edge²/area` | caller's limit |

The gate assumes a closed genus-0 shell per component, which is what a hull
without through-holes is. A torus-like or open surface must pass
`expected_euler`/`expected_components` or simply not be gated this way — the
gate refuses to guess.

## Resolution control

`resolve_resolution(parameters, target_volume_error=0.002)` doubles `nx` and `nz`
(the counts, so the facet aspect ratio stays roughly constant) until the
mesh-integrated displacement is within the requested relative error of
`analytic_volume`, or until `MAX_LONGITUDINAL_STATIONS`/`MAX_VERTICAL_LEVELS`/
`max_triangles` is reached. The returned `ResolutionPlan` carries the grid, the
measured error, the triangle/vertex counts, the refinement count and
`accepted`. An already adequate grid is returned unchanged (`refinements = 0`);
a target the budget cannot reach returns `accepted = False` with the best
measured grid instead of silently stopping.

Acceptance is against the closed-form integral
`V = L·B·T·(2/π)·√π·Γ(0.85)/(2·Γ(1.35))·(2/3)`, an integral of the analytic
generating functions. Refining against a finer mesh of the *same* generator shows
that the numbers move, not that they are right; the analytic integral is what
catches a uniformly under-resolved export.

`export_hull(..., target_volume_error=...)` wires the two together: it resolves
the grid, regenerates the model, gates it at the same tolerance, writes the
artifacts, and records the plan under `resolution` and the gate evidence under
`gate` in `metadata.json`. Passing a pre-built `metadata` document for a grid
that then has to change is rejected rather than silently overwritten.

## Measured ladder

| nx | nz | triangles | volume [m³] | error vs analytic |
|---:|---:|---:|---:|---:|
| 11 | 5 | 178 | 7.115242 | −5.290 % |
| 21 | 9 | 678 | 7.374345 | −1.842 % |
| 41 | 13 | 1998 | 7.446554 | −0.880 % |
| 81 | 21 | 6558 | 7.483649 | −0.387 % |
| 161 | 41 | 25918 | 7.502488 | −0.136 % |

Analytic volume: 7.512698 m³. The gate at the default 2 % tolerance rejects the
11×5 grid (`reference_volume`) and accepts 21×9 and finer.

With `target_volume_error = 0.002` from an 11×5 request the plan lands on
161×65 (41 278 triangles, −0.0759 % vs analytic, 4 refinements). The doubling
search overshoots: the error target is met with margin because the next grid is
4× denser, not because the tolerance was tuned. That is deliberate — a smaller
growth factor costs more generations for a tighter landing — and it is recorded
in the plan so a consumer can see it.

## The STL artifact

`write_stl` writes binary STL; `read_stl` rebuilds a `Mesh` by welding identical
float32 coordinates back into shared vertices. Two properties follow:

* **The topology is recoverable.** Two faces that shared a corner in double
  precision still share it after the deterministic double→float32 conversion, so
  the recovered soup is watertight and has `V - E + F = 2`. `export_hull(stl=True)`
  re-runs the whole gate on the file it just wrote and records the result under
  `stl` in `metadata.json`; a read-back failure is an error, not a warning.
* **Float32 is the floor.** Coordinates are quantised to 24-bit significands, so
  the recovered volume differs from the double-precision mesh by ~2.4e-8 relative
  on the reference hull. `STL_VOLUME_TOLERANCE = 1e-4` is that floor plus two
  orders of margin; the measured value is reported, not hidden.

The writer is byte-deterministic for a given mesh, ASCII STL is rejected with an
explicit message instead of being half-parsed, and a truncated or non-finite file
is refused.

## API

```python
from tensorcad import (HullParameters, export_hull, gate_mesh, generate_hull,
                       analytic_volume, resolve_resolution, read_stl)

model = generate_hull(HullParameters(nx=41, nz=13))
report = gate_mesh(model.mesh, reference_volume=analytic_volume(10.0, 2.0, 0.8))
report.passed                      # False for a coarse grid at the default 2 %
report.measurements["relative_volume_error"]

plan = resolve_resolution(HullParameters(nx=11, nz=5), target_volume_error=0.002)
plan.parameters                    # the grid that reaches 0.2 %
plan.accepted

artifacts = export_hull(model, "results/hull", stl=True, target_volume_error=0.002)
metadata = artifacts["metadata"]   # carries "resolution", "gate" and "stl" records
mesh = read_stl(artifacts["stl"])  # the file, re-welded and re-checkable
```

CLI:

```bash
python -m tensorcad.cli --nx 11 --nz 5 --target-volume-error 0.002 --stl --gate --output out
```

## Limitations

- The gate is a *geometry* gate: it proves the artifact is closed, single-shell,
  correctly oriented and resolved to the requested error. It knows nothing about
  fairness, resistance, stability or structural suitability — those are not
  implemented in this repository.
- Resolution control is analytic-reference based and therefore hull-family
  specific; imported meshes and the SUBOFF field can be gated structurally, but
  not against an analytic displacement.
- Only binary STL is read. ASCII STL is refused rather than converted.
- The gate checks topology and the integrated volume; it does not detect
  self-intersection, which would need a spatial index and is absent by design
  (no third-party dependencies in the core).