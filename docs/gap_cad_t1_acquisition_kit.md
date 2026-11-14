# GAP-CAD-T1-BASELINE-FLEET — digest-pinned acquisition of official hull geometry (mission m19948)

**Repository** TensorCAD · **Branch** `agent/m19948` · **Anchor** GAP-CAD-T1-BASELINE-FLEET
**New code** `tensorcad/fleet/acquire.py` · `tensorcad/fleet/data/acquisition_plan.json`
**Tests** `tests/test_fleet_acquire.py` (20 passed) · **Evidence** `runs/m19948_acquire/acquisition_evidence.json`
**Reproduce** `python -m tensorcad.fleet.acquire --out-dir runs/m19948_acquire --selftest`

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
