# GAP-CAD-T1-BASELINE-FLEET — digest-pinned acquisition of official hull geometry (mission m19948)

**Repository** TensorCAD · **Branch** `agent/m19948` · **Anchor** GAP-CAD-T1-BASELINE-FLEET
**New code** `tensorcad/fleet/acquire.py` · `tensorcad/fleet/data/acquisition_plan.json`
**Tests** `tests/test_fleet_acquire.py` (24 passed) · **Evidence** `runs/m19948_acquire/acquisition_evidence.json`
**Reproduce** `python -m tensorcad.fleet.acquire --out-dir runs/m19948_acquire --selftest`
**Re-measured 2026-11-15 on main** (sections 8–10, evidence `runs/m21978_acquire/acquisition_evidence.json`,
pin re-check `runs/m21978_pin_reverify/pin_reverification.json`): 2 of 6 hulls are now admitted,
KCS 227 → 230 m).**

## 1. What this closes

The earlier intake work (m17335) delivered the *judging* half of this anchor — registry, geometry
routes, ten gates, negative controls — and reported the *acquisition* half as **WITHHELD** on the
grounds that "this environment has no outbound network". Two facts changed that:

1. **The egress claim was wrong for this container.** `https://simman2014.dk/` answers HTTP 200 and
   serves the official hull packages. The block was a shell policy, not the network.
2. Even had no byte arrived, there was no code path to *receive* one: nothing could take an official
   file, check it against a pinned digest, write the provenance the gates demand, and hand it to a
   route. "Withheld" was therefore unfalsifiable — the status that never improves.

`tensorcad/fleet/acquire.py` is that path, under one rule enforced in code, not in prose:

> `verified = True` requires a digest **pinned before the bytes were read** and matching the bytes
> that were read.

All four ways of violating it are refusals with a reason, not warnings: no pin → `unpinned_digest`;
pin mismatch → `digest_mismatch`; dead transport → `source_unavailable`; oversized payload →
`size_exceeded`. A published ZIP is unwrapped only when it holds exactly one candidate member (the
archive digest pins the transport, the member digest pins the mesh). Admitted bytes then take the
**same** route and the **same** ten gates as a generated hull.

## 2. Result on the real fleet

> *Superseded for KCS/KVLCC2 by section 9 (re-measured 2026-11-15): both are admitted, and the
> KCS length below was measured again at 230.000 m after the x-origin defect of section 9.1 was
> fixed. The table is kept as the record of what m19948 measured.*

| ship | order | status | reason | provenance |
|---|---|---|---|---|
| KCS | SIMMAN 2014 `KCS-Hull.zip` → `KCS_hull_SVA.igs` (246 164 B) | **admitted** | all gates pass | `verified=true`, `retrieval=network_fetch`, archive sha256 `b50adace…f88c` |
| KVLCC2 | SIMMAN 2014 `KVLCC2-Hull.zip` → `kvlcc2.igs` (6 270 129 B) | fetched, digest-verified, **intake refused** | `geometry_route_refused`: *field produced fewer than two usable waterline slices* | archive sha256 `dc032647…d548` |
| ONRT, JBC, DTC, Series60 | — | withheld | `plan_incomplete` (no portal file URL located; DTC registry URL answers 404) | — |

KCS, measured on the **official IGES** through `iges_field` at `dx = 1.0 m, dz = 0.2 m, 3 samples/m,
200 ring points` (55 B-spline surfaces, 72 588 sample points, 22 000 triangles):

| quantity | mesh | declared/literature | error |
|---|---|---|---|
| displacement volume V | 52196.1 m³ | 52030.0 m³ | **+0.319 %** |
| length L | 230.000 m | 230.000 m | +0.0000 % |
| beam B | 32.2705 m | 32.2000 m | +0.2189 % |
| draft T | 10.8000 m | 10.8000 m | +0.0000 % |
| topology | 0 boundary edges, χ = 2, None non-manifold edges | — | exact |

This is the first ship in the fleet whose volume comparison is against **official geometry** rather
than against a surrogate built to hit the same `Cb`.

## 3. Two real defects found by running the path

* **The `iges_field` route was dead code.** `build.iges_field_mesh` iterated `iges.surfaces`, which is
  a *method* on the current reader, so the route raised `TypeError: 'method' object is not iterable`
  before touching a file. Any "imported IGES" claim in earlier fleet reports was therefore untested.
  Fixed (accepts both spellings); now exercised on two real files.
* **G8 measured draft from one end only.** `gate_principal_dimensions` computed `draft = abs(zmin)`
  when `zmax > 0`, i.e. it assumed a hull authored with z ∈ [−T, 0]. The official KCS IGES is authored
  keel-at-origin (z ∈ [0, 10.8]), so G8 measured **0.0 m of draft** and refused a correct hull. Draft
  is the z *extent*; the waterline's absolute position is a property of the file's frame, not of the
  ship. Fixed; L/B/T now agree to 0.219 % worst case.

Neither defect was reachable from generated geometry: the generated routes author z ∈ [−T, 0] and
never call the IGES path.

## 4. Why KVLCC2's intake is refused (honest FAIL, not a pass)

The 6.3 MB `kvlcc2.igs` was fetched and its archive digest matches the pin, so transport and identity
**are** verified. The field route then refuses it: `field_from_points` found fewer than two usable
waterline slices inside z ∈ [0, T] with the KCS-tuned parameters. Prime suspects, none confirmed in
this mission: the file's frame/origin (KVLCC2 packages are commonly authored with the waterline at
z = 0 and the hull below, or in a non-metre unit), or its patch structuring. The run records the
refusal with the reason instead of substituting a surrogate — exactly the behaviour the anchor
exists to enforce. Next step is a frame/unit diagnostic on that file, not a tolerance relaxation.

## 5. Pipeline controls (offline, machine-checked)

`--selftest` builds a surrogate KCS STL and pushes it through the same functions:

| control | expected | observed |
|---|---|---|
| pinned payload | admitted | admitted |
| one byte flipped | digest_mismatch | digest_mismatch |
| correct bytes, no pin (neither CLI nor plan) | unpinned_digest | unpinned_digest |
| unreachable source | source_unavailable | source_unavailable |

The first version of the "no pin" control **passed the wrong way** — it dropped only the command-line
pin and the plan's pin still applied — which is how the fallback rule was confirmed rather than
assumed. The surrogate is labelled a surrogate in the STL header and in the custody record: passing
these controls proves the pipe is real and says nothing about whether an official hull arrived.

## 6. Error judgement (explicit)

| quantity | value | criterion | verdict |
|---|---|---|---|
| KCS displacement volume vs published 52 030 m³ | **+0.319 %** | < 3 % | **PASS** |
| KCS L / B / T vs declared | 0.000 % / 0.219 % / 0.000 % | < 1 % | **PASS** |
| KCS topology | 0 boundary edges, χ = 2 | exact | **PASS** |
| KVLCC2 intake | refused with reason | — | **FAIL (honest)** |
| ships admitted from official geometry | 1 / 6 | 6 / 6 | **NOT MET** |

Is the key quantity within 3 %? **Yes for KCS (+0.319 %)**, and
within the 1 % the anchor demands. The anchor as a whole is **not** satisfied: one of six hulls now
has official geometry behind its volume comparison; five are withheld, and the plan states for each
exactly which URL and which digest is missing.

## 7. Caveats

* Pinned digests are **trust-on-first-use**: the digests of the bytes as served on 2026-11-14, not
  checksums published by the workshop. They prove a repeat retrieval returns the same file, not that
  the file is the benchmark hull. HTTP `ETag`/`Last-Modified` are recorded so the retrieval replays.
* The `max |y|` field binning biases the surface outward by O(cell²/radius); `field_cover` reports how
  much of the field was interpolated, and the +0.32 % volume is consistent with that bias.
* KVLCC2's extracted 6.3 MB `kvlcc2.igs` is not committed (repo hygiene); re-running with
  `--custodian KVLCC2=<path>` regenerates the block, and the digest in `acquisition_evidence.json` is
  what a reader must check.
* Appendages (rudder, propeller, duct) are deliberately outside the single-closed-shell gate.


## 8. Pin re-verification: turning trust-on-first-use into a measured property

Section 7 named the kit's weakest link: the pinned digests are the digests of the bytes as served
on one day, so they prove identity only against themselves. The kit *claimed* that "a repeat
retrieval returns the same file" — and nothing in the package ever tested that claim. It was
prose, and prose about reproducibility is exactly what this anchor exists to stop accepting.

`verify_pins` / `python -m tensorcad.fleet.acquire --verify-pins` is that test. It re-fetches every
pinned plan entry, builds no geometry (bytes, not minutes), and puts each entry in one of five
states — `reproduced`, `mismatch`, `unavailable`, `unpinned`, `no_source` — with `mismatch` failing
the command (exit 1) rather than decorating a healthy-looking summary. Transport failures are
recorded, not raised, because a re-check that dies on the first dead link cannot report on the
other five hulls.

Run against the packaged plan on 2026-11-15 (`runs/m21978_pin_reverify/pin_reverification.json`):

| ship | pin | re-fetch | bytes | ETag / Last-Modified | verdict |
|---|---|---|---|---|---|
| KCS | `b50adace…f88c` | HTTP 200 | 73 295 (archive) → `KCS_hull_SVA.igs` 246 164 B | `"5cb08c74-11e4f"` / Fri, 12 Apr 2019 13:02:44 GMT | **reproduced** |
| KVLCC2 | `dc032647…d548` | HTTP 200 | 1 345 336 (archive) → `kvlcc2.igs` 6 270 129 B | `"5cb08c74-148738"` / Fri, 12 Apr 2019 13:02:44 GMT | **reproduced** |
| ONRT, JBC, DTC, Series60 | — | — | — | — | `unpinned` (nothing pinned, nothing to check) |

Both pins reproduce a day after they were taken, with the publisher's `ETag` unchanged, so the
pins are not artefacts of a single corrupted transfer. They are **still** not workshop-published
checksums, and the tool says so in its own summary rather than in this paragraph: what has been
upgraded is reproducibility, not authenticity. The four unpinned hulls stay `unpinned` — a state
that is not a pass and not a failure, but a to-do list with a reason.

## 9. Re-measured on main, 2026-11-15: two hulls admitted, and a real defect behind KCS

Re-running the intake with the CLI settings (`dx` 1.0 m, `dz` 0.2 m, 3 samples/m, 200 ring points)
against the packaged plan gave 2 of 6 admitted — but the *first* run of the day gave 1 of 6, with
KCS refused. That refusal was not noise; it was the second real defect this kit has found in the
admission path.

### 9.1 The x-origin decision translated an AP-origin file

The published KCS IGES is authored with `x = 0` **at the aft perpendicular**: the model runs from
−6.0 m (transom overhang) to +237.84 m (bulb ahead of the FP), a 243.8 m extent around the Lpp
window it belongs to. `frames.normalise`'s `x_origin: auto` policy translated the aft-most point to
0 unconditionally, i.e. by +6.0 m. The half-breadth field window is `0 … Lpp`, so after that shift
the window started 6 m *inside* the hull and cut the forebody off:

| | before the fix | after the fix |
|---|---|---|
| applied `x` shift | +6.000 m | 0.000 m |
| measured length | 227.0 m (−1.30 % vs 230 m) | 230.000 m (0.00 %) |
| G8 verdict | **refused** | pass |
| displacement volume | 51 628.6 m³ (−0.77 %, would have passed G7) | 52 193.5 m³ (**+0.31 %**) |

The decision now *measures* where the file's own origin sits instead of assuming it. The share of
the model's extent that lies aft of `x = 0` separates the two published exports without a special
case: KCS keeps 2.5 % aft of its origin (the origin **is** the AP), KVLCC2 keeps 98.3 % (the origin
is at the bow, so the stern-negative export must be translated). The applied shift, the measured
fraction and the reason are written into `details["frame_normalisation"]["x_origin_reason"]`, so a
reader can re-check the decision from the evidence instead of trusting it.

### 9.2 Result

| ship | status | V_mesh | V_ref | ΔV | L | B | T | topology |
|---|---|---|---|---|---|---|---|---|
| KCS | **admitted** | 52 193.5 m³ | 52 030 m³ | **+0.314 %** | 230.000 (0.000 %) | 32.2365 (+0.113 %) | 10.800 (0.000 %) | 22 000 tris, 0 boundary, χ = 2 |
| KVLCC2 | **admitted** | 312 066.5 m³ | 312 622 m³ | **−0.178 %** | 320.000 (0.000 %) | 58.1285 (+0.222 %) | 20.800 (0.000 %) | 42 000 tris, 0 boundary, χ = 2 |
| ONRT, JBC, DTC, Series60 | withheld | — | — | — | — | — | — | `plan_incomplete` (no portal URL located; the plan names which URL and digest is missing) |

KVLCC2's intake, reported as an honest FAIL in section 4, now completes: the frame flip
(`z`-down, stern-negative) plus the trimmed-patch sampling that later work added let the field
route read the 6.27 MB file, and the mesh matches the published particulars to 0.22 % worst case.
This is the first time **two** hulls in the fleet carry official geometry behind their volume
comparison instead of a closure value (`Cb·L·B·T`) or a surrogate.

Figures (`runs/m21978_acquire/figures/`, generated by `runs/m21978_acquire/make_figures.py` from
the evidence and the delivered meshes, no hand-copied numbers):

* `fig1_half_breadth_fields.png` — the half-breadth field `b(x,z) = max|y|` rebuilt from each
  admitted mesh: KCS `0 … 230 m` and KVLCC2 `0 … 320 m`, both with the measured beam/draft extents;
* `fig2_errors_vs_reference.png` — V, L, B, T of both hulls against the published particulars,
  every bar labelled with its own number, with the 1 % and 3 % bands drawn.

### 9.3 Grid sensitivity (a caveat that is measured, not assumed)

Refining the field grid does **not** converge the volume to the reference from one side
(`runs/m21978_acquire/route_sensitivity.json`):

| ship | settings (dx, dz, samples/m, ring) | cover | ΔV |
|---|---|---|---|
| KCS | 1.0, 0.2, 3, 200 | 77.1 % | +0.314 % |
| KCS | 0.5, 0.1, 8, 240 | 88.8 % | **+0.033 %** |
| KVLCC2 | 1.0, 0.2, 3, 200 | 95.7 % | −0.178 % |
| KVLCC2 | 0.5, 0.1, 8, 240 | 72.6 % | **−0.633 %** |

The residual tracks `max|y|` binning (outward bias, O(cell²/radius)) against the trimmed-patch
sampling (inward bias where a patch is under-sampled), and the two do not cancel in a way that can
be predicted without measuring both files. Both settings pass the 1 % gate on both hulls, so the
admission decision is robust to the grid; the *last* digit of the volume is not.

## 10. Error judgement (re-measured, explicit)

| quantity | measured | reference | error | criterion | verdict |
|---|---|---|---|---|---|
| KCS displacement volume | 52 193.5 m³ | 52 030 m³ | **+0.314 %** | < 3 % | **PASS** |
| KVLCC2 displacement volume | 312 066.5 m³ | 312 622 m³ | **−0.178 %** | < 3 % | **PASS** |
| KCS L / B / T | 230.000 / 32.2365 / 10.800 | 230 / 32.2 / 10.8 | 0.000 / +0.113 / 0.000 % | < 1 % | **PASS** |
| KVLCC2 L / B / T | 320.000 / 58.1285 / 20.800 | 320 / 58 / 20.8 | 0.000 / +0.222 / 0.000 % | < 1 % | **PASS** |
| topology, both ships | 0 boundary edges, χ = 2, 0 non-manifold | exact | — | exact | **PASS** |
| pins reproducible (KCS, KVLCC2) | 2 / 2 | 2 / 2 | — | all | **PASS** |
| hulls admitted from official geometry | **2 / 6** | 6 / 6 | — | 6 / 6 | **NOT MET** |

Is the key quantity within 3 %? **Yes** — worst case is 0.314 % (KCS volume), five times inside
the 3 % criterion and inside the 1 % band this anchor demands. The anchor as a whole is still
**not** satisfied: two of six hulls have official geometry, four are withheld with a named missing
URL. What changed on 2026-11-15 is that the count moved 1 → 2 for a measured reason, and that the
number a reader is asked to trust (the pin) is now re-checked by a command instead of a sentence.
