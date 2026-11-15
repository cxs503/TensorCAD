# KVLCC2 frame policy — wiring the published IGES into the benchmark runner

Anchor **REPRO-GEOM-FLEET-REFS** (TensorCAD), acceptance ⑤: *"runner 接线 KVLCC2 已有参考使基准可复跑"*.
Status of this round: **implementation + frame evidence landed, unit-tested; the volume build is
`not_run`** (cluster job `cj3175` was still queued after 25 min — see §5). No number here is
presented as a benchmark result unless it was measured.

## 1. The gap being closed

The cross-reference library (`docs/fleet_refs_crossref.md`) already carries KVLCC2's published
displacement ∇ = 312 622 m³ and its Cb 0.8098 / LCB 3.48 %Lpp / Lwl 325.5 m, but §4.1 of that
document shows `—` in the "got ∇" column: the *runner* never compared anything with the reference.
The reason was never the reference — it is the **file frame** of the published IGES.

## 2. What was implemented

| file | change |
|---|---|
| `tensorcad/fleet/frames.py` (new) | `detect()` measures the sampled surface (axis ranges, depth below/above the file's own `z = 0`, half-model test, implied beam); `normalise()` applies a declarative `FramePolicy(z_axis, x_origin)` and returns the transformed points **plus the exact flip/shift it applied**. |
| `tensorcad/fleet/build.py` | `iges_field_mesh(..., frame=..., x_extent=...)`: the route frames the samples *before* binning the half-breadth field, resolves `asset://` pointers through the asset store, and records `frame_normalisation` + `x_extent` in the build details. `x_extent="geometry"` lets the field follow the model's own extent instead of clipping to `Lpp`. |
| `tensorcad/fleet/manifest.py` | `ShipRecord.geometry_frame` (declarative policy, round-tripped through the manifest loader and `as_dict`), so a frame policy is auditable per hull rather than being a hidden default. |
| `tools/baseline_fleet/kvlcc2_frame_build.py` (new) | End-to-end reproducible path: packaged KVLCC2 record → asset-store IGES → frame policy → watertight solid → the *same* `evaluate_gates` set as the fleet run → JSON artefact, Markdown table and figure, all generated from the run (no hand-transcribed digits). |
| `tests/test_fleet_frames.py` (new) | 11 tests: keel-origin vs waterline-origin detection, the depth-ratio trap (§4), forced policies, invalid policies, empty cloud, manifest round-trip, dispatch passing the policy to the IGES route, and the guard that the packaged offline fleet stays free of IGES routes. |

The packaged offline fleet is deliberately **not** switched to the IGES route
(`tests/test_fleet_official.py` requires that), so the KVLCC2 reference stays reproducible on a
machine that has the asset, exactly as the "official fleet" path already does.

## 3. Measured frame evidence (real, from the published file)

Read from `asset://geometry/hulls/KVLCC2/kvlcc2.igs` (sha256 `ed06c287…8e2f9`, 6 270 129 bytes,
19 B-spline patches, 10 473 samples at 2.0 m patch spacing):

| quantity | value | consequence |
|---|---|---|
| x range | −328.116 … 5.502 m (extent 333.6 m) | `x = 0` sits near the **bow**; clipping the field to `0 … Lpp = 320 m` drops the whole hull |
| y range | −0.416 … 29.374 m, 99.6 % of samples at `y ≥ 0` | mirrored **half model**; `max abs(y)` is the half-breadth, so B = 58.0 m is consistent (and the SIMMAN page's `Bwl 580` is a decimal-point typo, already flagged by the library's closure test) |
| z range | −28.067 … +0.055 m | `z = 0` is the **design waterline**, z grows **downwards**; the field route keeps `0 ≤ z ≤ draft`, i.e. nothing (−9 % of samples fall in that band) |
| samples below `z = 0` | 37.4 % | a *sample-fraction* rule would call this file z-up — see §4 |

Consequence: before this change the KVLCC2 IGES could not produce a single waterline slice, which
is why the benchmark refused it. The hull shape was never the problem.

## 4. The trap that the tests now pin

The obvious detection rule — "more than half the samples below `z = 0` ⇒ z-down" — is **wrong on
this file**: only 37 % of samples are below zero, because the above-water sliver is tiny in extent
but densely sampled per patch. The criterion used instead is a *depth ratio*:
`below_zero = 28.067 m` vs `above_zero = 0.055 m` ⇒ z-down; a keel-origin export has
`below_zero ≈ 0` ⇒ z-up. `tests/test_fleet_frames.py::test_depth_ratio_beats_the_sample_fraction_on_the_real_export`
locks this in with the measured pattern (fraction < 0.5, ratio ≫ 3).

## 5. Honest status of the number

* `not_run`: the ∇ measurement itself. The build script (`tools/baseline_fleet/kvlcc2_frame_build.py`,
  `--resolutions 1.0,0.5`) was submitted to the compute gate as job `cj3175`; after 25 min of waiting
  the queue still reported `queued` (33 jobs in flight), so this round publishes **no** KVLCC2
  reconstruction error. Per the honest-FAIL rule nothing is filled in by hand.
* Therefore the acceptance question "is the reconstructed ∇ within 3 % of 312 622 m³" **cannot be
  answered yet** — the answer must come from the artefact `benchmarks/verified/fleet_refs/kvlcc2_frame.json`
  that the script writes, and the document it regenerates is the place to read it.
* Verification actually performed this round: `python -m pytest tests/test_fleet_frames.py -q`
  → **11 passed**; the frame evidence in §3 was measured directly from the published file.

Re-run (on a host with the asset, or through the compute gate):

```bash
python tools/baseline_fleet/kvlcc2_frame_build.py \
    --resolutions 1.0,0.5 --spacing 0.5 \
    --out benchmarks/verified/fleet_refs/kvlcc2_frame.json \
    --doc docs/kvlcc2_frame_policy.md --figs docs/figs/fleet_refs
```

## 6. Known limits

* The `max`-binned half-breadth field biases the surface outwards by O(cell²/R); the sign and size
  of the residual at `dx = 1.0 / 0.5 m` is what the two-resolution run is for, and a single
  resolution would not distinguish a field bias from a shape error.
* LCB from this route is only as precise as the field cell; a 0.1 %Lpp claim needs a finer field.
* The KVLCC2 `Lpp` (320 m) is shorter than the model's own extent (333.6 m); with
  `x_extent="geometry"` the volume is integrated over the model, and the Cb reported by the script
  is normalised with `Lpp·B·T` as published — the two conventions differ and both are printed, so
  the comparison with 312 622 m³ (a volume, frame-independent) stays valid.