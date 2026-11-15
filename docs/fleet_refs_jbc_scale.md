# JBC model scale: derived from the CAD bytes vs a published statement

Anchor **REPRO-GEOM-FLEET-REFS** (TensorCAD). Machine artefact:
`benchmarks/verified/fleet_refs/model_scale.json` (regenerate:
`python tools/baseline_fleet/model_scale_refs.py build`). Every number below is rendered
from that artefact; nothing is transcribed by hand.

## 1. The claim under test

`fleet_specs.json` turns the admitted JBC geometry into full-scale metres with a **derived**
factor `λ = 39.7285`, because the Tokyo 2015 JBC geometry page publishes no scale row.
The bench compares displacement volumes, so that factor enters cubed: a 0.7 % scale error
is a 2.1 % volume error, i.e. twice the 1 % acceptance band of the baseline benchmark.

## 2. Documents read

| key | source | licence | sha256 | claims read |
|---|---|---|---|---|
| jbc_scale_pr | Research Square preprint (Tarafder, Hasanat, Hasan, Sayeed, Hannan), posted 2025-09-30, 'Scale Effects on Bare-Hull Resi | CC BY 4.0 (stated on page 1 of the document) | `4016a43201d350af…` | beam, cb, declared_scale_ratio, displacement_volume, draft, full_length_pp, model_length_pp, model_length_statement |

## 3. Quoted statements

| ship | claim | value | unit | page | quote |
|---|---|---|---|---|---|
| JBC | full-scale B (m) | 45.0 | m | 8 | `Parameter Full Scale Model Scale Unit Length, Lpp 280.0 7.00 m Breadth, B 45.0 1.125 m Depth, H 25.0 0.625 m Draft, T 16.50 0.4125 m Displacement, 178369.9 2.7` |
| JBC | Cb | 0.858 | - | 8 | `m Draft, T 16.50 0.4125 m Displacement, 178369.9 2.7872 m3 Wetted Surface, S 19556.1 12.201 m2 Block coefficient, Cb 0.858 0.858 Design Speed, V 7.459 1.179 m/` |
| JBC | declared scale ratio λ | 40 | - | 7 | `, is also another benchmark bulk carrier with a full-scale design length of 280 m and the model scale 7.00 m with a scale ratio of 40. The geometry file for the` |
| JBC | full-scale ∇ (m³) | 178369.9 | m3 | 8 | `ength, Lpp 280.0 7.00 m Breadth, B 45.0 1.125 m Depth, H 25.0 0.625 m Draft, T 16.50 0.4125 m Displacement, 178369.9 2.7872 m3 Wetted Surface, S 19556.1 12.201` |
| JBC | full-scale T (m) | 16.50 | m | 8 | `Scale Model Scale Unit Length, Lpp 280.0 7.00 m Breadth, B 45.0 1.125 m Depth, H 25.0 0.625 m Draft, T 16.50 0.4125 m Displacement, 178369.9 2.7872 m3 Wetted S` |
| JBC | full-scale Lpp (m) | 280.0 | m | 8 | `Parameter Full Scale Model Scale Unit Length, Lpp 280.0 7.00 m Breadth, B 45.0 1.125 m Depth, H 25.0 0.625 m Draft, T 16.50 0.4125 m Displacement, 178369.9` |
| JBC | model Lpp (m) | 7.00 | m | 8 | `Parameter Full Scale Model Scale Unit Length, Lpp 280.0 7.00 m Breadth, B 45.0 1.125 m Depth, H 25.0 0.625 m Draft, T 16.50 0.4125 m Displacement, 178369.9` |
| JBC | model length from the same sentence (m) | 7.00 | m | 7 | `, is also another benchmark bulk carrier with a full-scale design length of 280 m and the model scale 7.00 m with a scale ratio of 40. The geometry file for the` |

## 4. Derived vs published scale

### JBC

* derived scale (fleet_specs): **39.728496** — model scale derived from the publisher CAD (no scale row is published): B=45 m / D=25 m over the measured model spans, half-spread 0.18%; record artifacts/baseline_fleet/jbc_model_scale.json (page sha256 fc2ce9054a8b7842…, archive sha256 e19800ee954f4960…)
* published scale: **40.000000** — sentence λ=40 and Table 3 280 m / 7 m = 40 agree to 0.000 %
* discrepancy: **+0.683 %** in λ = **1.02064×** in ∇

| quantity | at derived scale | at published scale |
|---|---|---|
| reconstructed ∇ (m³) | 145769.0 | 148778.1 |
| error vs published 178369.9 m³ | -18.28 % | -16.59 % |

The scale statement moves the comparison by **+1.69 percentage points**; the benchmark verdict stays **WITHHELD** and the residual -16.59 % is geometric, not a unit error.

## 5. What this does and does not establish

* It **does** show that the derived λ is not the published λ: the admitted JBC CAD spans imply
  1.1327 m where a λ=40 model would measure 1.125 m across the published B = 45 m, and the
  difference is quoted, digests pinned, and priced in the benchmark's own currency.
* It does **not** decide which number is right. The CAD spans are measured; the λ=40 statement
  is published. Both stay in the record (`geometry_scale` untouched, the statement quoted in
  `scale_cross_check`), and no volume comparison is silently rescaled.
* It does **not** upgrade the JBC verdict. The withheld −18 % is dominated by the IGES
  control-net layout guard, not by the scale; the scale explains only the number of points
  tabulated in section 4.

## 6. Limitations

* The published statement is a 2025 preprint, not the workshop's own page; its full-scale
  column reproduces the Tokyo 2015 table to the digit, which is why it is treated as a
  transcription witness for the particulars and as the source of the *scale* statement only.
* The model column of that table (Lpp 7.00 m) is itself consistent with λ = 280/7 = 40; the
  tool requires the sentence and the table to agree within 0.5 % before using either.
* No published scale exists for the remaining hulls (KCS/KVLCC2/ONRT use published ratios;
  DTC/Series 60 have no admitted geometry in this fleet), so the comparison is one-hull.
