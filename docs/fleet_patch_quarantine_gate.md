# G9: the patch quarantine has to reach the verdict

**Defect class:** fail-open in the baseline-fleet acceptance chain.
**Status:** fixed; read side unchanged, verdict side gated.
**Evidence:** `runs/m22653/quarantine_gate_ab.json`,
`docs/figs/patch_quarantine_ab.png`, `tests/test_fleet_patch_quarantine_gate.py`.

## 1. What the reader already did, and where it stopped

The type-128 reader fails closed. When a B-spline patch's control-point tail can be
split into control points two different ways and the evidence does not decide it — the
canonical case is a degenerate `2x2` (or `2xn`) net, where every second difference is
identically zero, so control-net roughness scores `0.0` for *both* candidate layouts —
the reader does not guess. It returns the patch with `layout_ambiguous = True`, and
`tensorcad.fleet_baseline.partition_surfaces` quarantines it
(`docs/fleet_geometry_kernel_diagnostics.md`).

That is where the chain stopped, in two places:

1. `tensorcad.fleet.build.iges_field_mesh` — the route the **baseline fleet acceptance
   actually runs** (`fleet/run.py` → `build_ship_geometry` → `iges_field_mesh`) — never
   called `partition_surfaces`. It sampled every type-128 surface, including the ones the
   reader had refused to decide, so the quarantine was honoured by the legacy
   `fleet_baseline.build_ship` and by the diagnostic tools, and bypassed by the route
   that produces the fleet verdicts.
2. The gate chain (`tensorcad.fleet.gates.evaluate_gates`, G0–G8) had no parameter that
   could carry a patch census at all, and no committed verdict carries one.

The consequence is silent by construction. The IGES route does not mesh patches: it
samples them into a half-breadth field `b(x, z) = max |y|` and lofts that field into a
watertight mesh. A dropped patch is therefore **not** a hole in the mesh — the field
interpolates across it from its neighbours, the mesh stays watertight, the Euler
characteristic stays 2, and G7's volume is the volume of a hull that is missing surface.
Measured: dropping an entire patch changes the volume by exactly `0.0 m3`
(`runs/m22653/quarantine_gate_ab.json`, `volume_difference_m3`).

## 2. The repair

* `fleet/build.py` — `iges_field_mesh` partitions the surfaces before sampling, refuses
  the file when every patch is quarantined, excludes quarantined patches (and, on the
  trimmed route, trimmed patches whose *basis* surface is ambiguous) from the field, and
  reports the census in `details["patch_quarantine"]`
  (`patches_total` / `patches_used` / `patches_quarantined`).
* `fleet/gates.py` — new gate **G9 `patch_quarantine`**: a hull with more dropped patches
  than `QUARANTINE_TOLERANCE` (default **0**) does not pass, whatever its volume error.
  `quarantine_census(details)` reads the census from either builder generation (the new
  `patch_quarantine` block or the legacy `source` block), and returns `None` — not `[]` —
  when a route reports no census, so the verdict can say `geometry_completeness =
  "unknown"` instead of pretending "complete".
* `fleet/run.py`, `fleet/acquire.py`, `fleet/replace.py` — every call site that has a
  `BuiltGeometry` forwards its census into `evaluate_gates`.

G9 is a *completeness* gate, not an accuracy gate: it does not say the volume is wrong,
it says the geometry the other gates graded is not the geometry the file describes. That
distinction is why it must be a hard gate and not a note: a note in a report is not read
by the next consumer, a failed gate is.

## 3. Measurements

### 3.1. Controlled A/B — same mesh, same volume, census is the only change

A synthetic barge (L = 40 m, B = 10 m, T = 4 m, analytic 1600 m3) described by three
type-128 patches: two clean bilinear patches and one deliberately undecidable patch
(`2x2`, weights `1, 1.06, 0.94, 1`, weight test silent) that overlaps surface the other
two already carry.

| build | patches | used | quarantined | V_mesh (m3) |
|---|---|---|---|---|
| with the undecidable patch | 3 | 2 | 1 | 1519.8418 |
| without it | 2 | 2 | 0 | 1519.8418 |

`volume_difference = 0.0 m3` — the lost patch is invisible to G7, as predicted.

Grading that one mesh twice (declared volume set to the measured volume, `closure`
basis, so that only the census can differ):

| verdict | passed | failed gates | `geometry_completeness` |
|---|---|---|---|
| without census (G1–G8) | **True** | — | `unknown` |
| with census (+G9) | **False** | `G9_patch_quarantine` | `patches_dropped` |

The barge's `1519.84 m3` against the analytic `1600 m3` is `-5.0 %`, the field/loft
discretisation of the square ends at `dx = 2 m`; it is not part of the claim above (the
two builds agree to `1e-12`).

### 3.2. Is the hole live on today's four-hull fleet? Cross-tabulation

Cross-tabulating the quarantined-patch counts recorded in
`benchmarks/verified/baseline_fleet/geometry_diagnostics.json` with the verdicts in
`benchmarks/verified/baseline_fleet/result.json` (evidence that predates G9):

| hull | quarantined / patches | verdict |
|---|---|---|
| KCS | 0 / 55 | VERIFIED |
| KVLCC2 | 1 / 19 | WITHHELD |
| JBC | 4 / 45 | WITHHELD |
| ONRT | 26 / 132 | WITHHELD |

**No hull was passed while carrying quarantined patches.** The class is therefore
*latent* on this fleet, not observed: the three hulls with dropped patches were already
refused for other reasons (KVLCC2: fewer than two usable waterline slices; JBC: −18.28 %
volume; ONRT: model-scale bytes, extent 42.5 m vs `L_PP` 154 m). Saying "we fixed a
defect that was firing" would be false; the honest statement is that the chain had no
mechanism to notice, and one of the four hulls (KCS) passes only because it happens to
have zero quarantined patches.

## 4. Caveats and scope

* **The fleet numbers cannot be re-measured here.** The publisher IGES bytes are not in
  this repository (`assets` are resolved through `asset://`, absent in a clean
  checkout), so the effect of excluding JBC's 4 and ONRT's 26 patches from their fields
  is *unmeasured*. Those hulls are `WITHHELD` either way, so no committed verdict flips,
  but their recorded volume errors will change when the fleet is next run against real
  bytes. That is the correct direction: they were computed from fields that included
  mis-read patches.
* **The undecidable patch in the A/B is synthetic**, by design: it isolates the plumbing.
  It is not evidence that any specific publisher patch is mis-read.
* **Identity matching on the trimmed route** is by basis-surface object handed out by
  `IgesFile.trimmed_surfaces()`, not by object identity from `IgesFile.surfaces()` —
  the reader re-parses each entity on demand, so the two lists contain distinct objects.
  The census for the trimmed route is therefore taken over the trimmed patches
  themselves.
* **`quarantine_tolerance` exists and defaults to 0.** A caller that has measured the
  impact of a dropped patch (e.g. proved the field is unaffected) may raise it; it cannot
  happen by accident.
* **Not covered here:** `tensorcad/fleet/official_geometry.py` builds its own verdict
  (it uses `fleet_baseline.build_ship`, which quarantines, but writes its own
  `verdict` field rather than running `evaluate_gates`). Threading G9 through that path
  is a separate, smaller change and is deliberately not bundled into this one.

## 5. Reproduction

```bash
python -m pytest tests/test_fleet_patch_quarantine_gate.py -q     # 14 tests
python tools/baseline_fleet/quarantine_gate_ab.py --out runs/m22653
```