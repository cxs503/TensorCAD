# STEP export (optional) and the TensorFEM handoff

Two ways out of the geometry core, for two different consumers:

| Path | Module | Extra | Consumer | What it is |
|---|---|---|---|---|
| surface exchange | `tensorcad.adapters.fem` | none | TensorFEM, any FEM preprocessor | node/element tables + named boundary patches, gmsh 2.2 ASCII and JSON |
| solid exchange | `tensorcad.adapters.ocp` | `pip install "tensorcad[ocp]"` | CAD/CAE tools that read STEP | the same triangles sewn into a closed solid, checked and written by OpenCASCADE |

Neither is imported by `tensorcad` itself. `import tensorcad` still loads nothing but the
standard library, and `tests/test_public_surface.py` runs that check in a fresh
interpreter.

Everything below is measured by `python examples/step_export/run_example.py`, which runs
the handoff part always and the STEP part only when OCP imports. `tests/test_adapters_fem.py`
covers the handoff in the normal test lane; `tests/test_adapters_ocp.py` skips itself
unless the wheels are installed.

## The handoff, measured

Input: `to_mesh(SuboffFullField(), (72, 9, 13), margin=2)` — 7,386 vertices, 14,764
triangles, watertight, sampled in 0.41 s.

| Patch | triangles | area m² | share |
|---|---|---|---|
| hull | 14,188 | 5.90651 | 96.749 % |
| sail | 281 | 0.10520 | 1.723 % |
| fins | 295 | 0.09326 | 1.528 % |

Volume 0.694708 m³, wetted area 6.104978 m² (both lower than the analytic
0.70754064 m³ / 6.350646 m² because 72 cells along a 4.36 m hull gives a 60 mm cell).

Round trips:

- gmsh 2.2: 0.89 MB, 7,386 nodes, 14,764 triangles, three physical names; **worst node
  displacement 0.0 m** (floats are written with `repr`) and worst per-patch area
  difference 2.2 × 10⁻¹⁴ m².
- JSON: 1.49 MB, digest identical to the surface it came from — the table digest and
  `SampledField.mesh_digest()` use the same SHA-256 convention, so a consumer holding only
  the file can check it against the producer's log.

The patches are assigned geometrically: each triangle centroid goes to the analytic
component whose signed distance is closest to zero, and the worst centroid distance to
its own component is reported (30.5 / 31.7 / 30.2 mm for hull / sail / fins at this
resolution, against a 60 mm cell). That number is the honest measure of "the union and
its components agree": if it ever approached a cell diagonal, the tags would be noise.

What the format does **not** carry: a size field, a boundary-layer prism stack, or a
volume mesh. It is a surface for a solver that meshes the volume itself, or for a BEM
code that never needs one.

## STEP, measured

Same surface, `write_step(mesh, path, length_scale=1000.0, round_trip=True)`:

| Quantity | Measured |
|---|---|
| sewing | 12.8 s for 14,764 triangles |
| writing | 4.4 s, 37.14 MB |
| free edges / multiple edges | 0 / 0 |
| closed shells | 2 (see below) |
| `BRepCheck_Analyzer` | valid for both solids |
| kernel volume | 0.694708267 m³ |
| mesh signed volume | 0.694708267 m³ (**relative difference +5.3 × 10⁻¹⁵**) |
| re-imported volume | 0.694708267 m³ (relative difference −3.1 × 10⁻¹³) |
| faces in the file | 14,764 — one planar face per triangle |

The volume parity is the point. OpenCASCADE integrates a B-rep; TensorCAD integrates a
triangle soup; both are exact for a closed surface, so agreement to round-off is the
check that the sewn shell is the same solid. A hole, a fold or a dropped face shows up
as a discrepancy at once, which is why the exporter measures instead of asserting.

**The two shells are a real observation.** At 72 cells along the hull, marching tets
leaves a 48-triangle island behind a fin root (volume 1.6 × 10⁻⁵ m³), and sewing turns
each closed island into its own shell. The exporter builds one solid per shell, writes
them as a compound, and reports `shells` and `solids`, so nothing is silently dropped —
dropping the island would have changed the volume in the file. At (100, 12, 18) and
(150, 18, 26) the same field produces a single component, so this is a sampling
artefact of the coarse grid, not a geometric feature.

## What STEP export here is not

A faceted solid is not a B-rep model of a submarine. The file contains 14,764 planar
faces and no analytic surfaces, no fillets, no feature tree, no parametric history. That
is what a solver-facing exchange file needs and what a designer does not want. Getting a
feature model means rebuilding the hull, sail and fins as swept/lofted surfaces in a CAD
kernel — a different project, tracked as [issue #4](https://github.com/cxs503/TensorCAD/issues/4).

The exporter is also deliberately one-directional in practice: `read_step` exists to
close the loop (it parses the file back and measures it), not to import arbitrary STEP
into the field protocol. A STEP file of an unknown ship is not a
`GeometryField`.

## Reading the numbers the exporter returns

```python
from tensorcad.adapters.ocp import have_ocp, read_step, write_step

record = write_step(mesh, "results/aff8.step", length_scale=1000.0, round_trip=True)
record.passed                       # free edges 0, shells == solids, closed, valid, V > 0
record.volume_m3                    # integrated by the kernel, in m³
record.round_trip_volume_error      # re-imported volume / written volume - 1
record.shells, record.solids        # one solid per closed shell, none dropped
record.to_dict()                    # JSON-ready, with timings and the scale used
```

`length_scale` multiplies the coordinates on the way out (default 1000, i.e. millimetres)
and divides them back on re-import; the value used is recorded, never assumed. The sewing
tolerance defaults to `1e-9 · length_scale`, a nanometre in model units, far below any
cell in this repository, so it can only join vertices that are meant to be identical.

When OCP is missing, `write_step` raises `OcpUnavailable` with the pip command in the
message — it does not fall back to writing a mesh file with a `.step` extension, which
would be worse than failing.
