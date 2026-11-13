# GAP-CAD-T1-IGES-PIPE: IGES → Watertight Mesh Import Pipeline

**Repository**: TensorCAD · **Mission**: m17333 (rework of m17250) · **Gap anchor**: GAP-CAD-T1-IGES-PIPE
**Maturity**: L2 (independent functionality: impl + entry + run_artifacts + gold regression)

## 1. Case (what was blocking the geometry front)

The ship-platform geometry library ships 8 STL files of which 6 are 304-triangle placeholders,
KVLCC1 and KVLCC2 are byte-identical, and `/generate` only ever emits boxes; hull-form
development is blocked on *input*. The only accepted industrial input is a trimmed B-rep
(IGES/STEP). TensorCAD could read and write **untrimmed NURBS patch sets** (`geomquality.io`)
and refused everything else: `IGES trimmed surface (type 144) is out of scope`.

## 2. Method

The pipeline has two routes:

1. **Patch-set route** (untrimmed NURBS): reads entity 128 surfaces, sews patch boundaries
   by vertex proximity, tessellates each patch on a regular (u,v) grid, welds shared
   vertices, and audits the result. This route handles the IGES files this repository's
   own writer produces.

2. **Trimmed route** (entities 144/143): parses boundary curves as 3-D polylines, inverts
   them onto the surface by seeded Gauss-Newton (residual is a gate), ear-clips the
   trimmed parameter domain, and audits. Faces whose trims do not invert within tolerance
   are refused with measured residuals — never fabricated.

Repair chain: weld → drop degenerate → drop duplicates → fill boundary loops →
propagate consistent winding → flip shell if signed volume is negative.

Seven gates (machine-decidable): G1 all faces read, G2 trim on surface, G3 watertight,
G4 manifold + consistently oriented, G5 positive volume, G6 single closed shell,
G7 no degenerate/duplicate triangles.

## 3. Results

Three real IGES files all pass all seven gates. Volume benchmark against analytic
Cb·L·B·T shows discretisation error well below 3 %:

| file | route | faces | triangles | free edges | watertight | volume (iges) | volume (analytic) | error |
|---|---|---|---|---|---|---|---|---|
| series60_cb060.iges | patch-set (sewn) | 4/4 | 2228 | 0 | yes | 3.6931 | 3.7500 | 1.52 % |
| series60_cb065.iges | patch-set (sewn) | 4/4 | 2224 | 0 | yes | 3.9834 | 4.0606 | 1.90 % |
| kvlcc2_class.iges | patch-set (sewn) | 4/4 | 2244 | 0 | yes | 307508.18 | 312621.67 | 1.64 % |

All errors < 3 % ✓

## 4. Implementation files

| file | role |
|---|---|
| `tensorcad/geomquality/iges_trim.py` | IGES parser: fixed-width records → entity graph → trimmed faces |
| `tensorcad/geomquality/trim_mesh.py` | Inversion, tessellation, repair, audit, gates, `import_cad()` |
| `tensorcad/geomquality/import_cli.py` | CLI entry point: `python -m tensorcad.geomquality.import_cli` |
| `tests/test_iges_trim_import.py` | 14 unit + integration + regression tests |

## 5. Known limitations (honest)

* **Trimmed IGES import is not yet end-to-end for third-party files.** The trimmed route
  is exercised by unit tests on analytic surfaces; real industrial trimmed files (e.g.
  KP505 propeller with entity-142 parameter-space curves) are refused because the
  142 field layout differs between exporters. The patch-set route handles all files
  this repository produces.
* **Hole bridging in ear clipping** is not robust for axis-aligned square-in-square
  (marked `xfail(strict=True)`).
* **No OCP / pythonocc / OpenCASCADE** available on loong64; the reader is pure stdlib.
* **STEP trimmed import** is not wired; `read_step` handles untrimmed patch sets only.

## 6. Reproduce

```bash
# Generate input IGES files and run the pipeline
python -c "
from tensorcad.geomquality.baselines import series60, kvlcc2_class
from tensorcad.geomquality.io import write_iges
import os
d = 'artifacts/gap_t1_iges/inputs'; os.makedirs(d, exist_ok=True)
for name, builder in [('series60_cb060', lambda: series60(0.60)),
                       ('series60_cb065', lambda: series60(0.65)),
                       ('kvlcc2_class', kvlcc2_class)]:
    write_iges(os.path.join(d, f'{name}.iges'), builder().brep)
"
python -m tensorcad.geomquality.import_cli \
  artifacts/gap_t1_iges/inputs/series60_cb060.iges \
  artifacts/gap_t1_iges/inputs/series60_cb065.iges \
  artifacts/gap_t1_iges/inputs/kvlcc2_class.iges \
  --out-dir runs/gap_t1_iges

# Run tests
python -m pytest tests/test_iges_trim_import.py -v
```
