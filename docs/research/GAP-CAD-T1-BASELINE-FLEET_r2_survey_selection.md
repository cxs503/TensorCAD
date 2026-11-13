# GAP-CAD-T1-BASELINE-FLEET — round 2 survey: open-source landscape, method lineage, selection

Mission m18641 (round 2 of the two-round protocol; round 1, m18602, ended honest-FAIL with no
stage recorded). This document is the survey half: what the open-source landscape offers for a
*validated baseline hull fleet*, what the methodological literature prescribes, and which pieces
this repository adopts. The migration plan and the numbers of the implementation are in
`docs/gap_cad_t1_baseline_fleet_r2.md`; the machine-readable selection record is §4 here
(`oss_solutions`-shaped, `feature_slug = GAP-CAD-T1-BASELINE-FLEET`).

## 1. Open-source horizontal survey

### 1.1 What was actually available to read

The mission brief points at the cluster `~/oss-cache` inventory (`oss_inventory` where
`feature_slug ILIKE '%FLEET%'`). **Result: nothing to learn from there, and that is recorded
rather than papered over.** `~/oss-cache` on this box holds three source trees — `pylbm`,
`Maarten-vd-Sande_lbm`, `hschwane` — all lattice-Boltzmann; there is no CAD/geometry-baseline
project cached, and no `oss_inventory`/`oss_solutions` table is reachable from this worktree
(the platform DB files present locally contain no `oss%` tables, and the sandbox blocks the
egress that would let me query the platform). So the survey below is a *documented reading* of
each project's published behaviour, and where I could not execute a project I say so instead of
claiming experience with it. This is the main limit of round 2: no third-party hull-baseline
implementation was executed here.

### 1.2 Projects whose structure is worth copying, and what exactly to copy

| project | licence | what it does that we need | what we adopted |
|---|---|---|---|
| **trimesh** (mikedh) | MIT | `is_watertight`, `is_winding_consistent`, `euler_number`, volume via divergence theorem; the everyday vocabulary of "is this mesh usable" | the *gate semantics*: watertight = 0 boundary edges, manifold + winding-consistent checked separately, Euler 2 for a single closed shell. Our `gates.py` measures the same four facts with its own code path (no dependency added: the repo must audit meshes without a third-party install). |
| **manifold3d** (Emmett Lalish) | Apache-2.0 | guarantees a manifold output solid so that boolean ops downstream cannot produce non-manifold input | the ordering principle: *the intake is where manifoldness is decided*, not the solver. Nothing imported. |
| **OpenCASCADE / OCCT** (with FreeCAD as consumer) | LGPL-2.1 + exception | STEP/IGES reading, `BRepCheck_Analyzer`, `GProp` volume properties — the reference implementation for "load the official hull file and measure it" | the *route design*: an `imported_stl` route and an IGES field route already exist in `build.py`; OCCT is not a dependency, so an official file is parsed by our own reader and the gates take responsibility for the result. |
| **Gmsh** | GPL-2.0 | `-check` style topology diagnostics before meshing | the habit of reporting *numbers with tolerances* in the diagnostic, which our `GateResult.detail` follows. |
| **SU2 `TestCases` / OpenFOAM tutorial regression harnesses** | LGPL / GPL | a *fleet* of cases with per-case expected values and a driver that refuses a case whose result drifted | the fleet-level unit of analysis: per-ship manifests plus negative controls that must fail, so a green run cannot be green by accident. |
| **Ship-D (Bagazinski & Ahmed, IDETC-CIE 2023) dataset** | public research dataset | 30 000 parametric hulls with validity filters (the parent form of this repo's `series60`/parametric routes) | the principle that a *parametric* hull is only admissible with its parameter window declared; our `series60_form` route refuses `Cb` outside its verified window instead of clipping. |
| **Software Heritage / FAIR / W3C PROV / RO-Crate** | — (standards) | content-addressed identity (sha256 / SWHID) and provenance graphs | every delivered artefact is identified by the sha256 of its bytes and every manifest digit by a source anchor; the ledger is a small PROV-like graph (file → sha → record → source anchor). |

### 1.3 What no open-source project provides (the gap this anchor fills)

None of the above answers the question the owner's defect poses: *"is this collection of hull
files internally consistent and externally anchored?"* Mesh libraries validate one file; test
suites validate one solver case. A **fleet** gate needs cross-file analysis (byte-identical
geometry across two ships), record-vs-geometry analysis (mesh volume vs the published
displacement), and provenance analysis (does every digit have a named source?). That triple is
what `tensorcad/fleet/` implements, and round 2 adds the two pieces that were missing: the
reference-triple cross-check and the replacement delivery.

## 2. Method lineage (what the literature prescribes, and where it is used here)

| decision | standard practice | source (external anchor) |
|---|---|---|
| benchmark hulls | methodical-series / workshop test cases are the accepted baselines; a hull without a published displacement and particulars is not a baseline | Todd (1963) DTMB Report 1710 (Series 60); Gothenburg 2010 workshop (eds. Larsson, Stern, Visonneau, Springer 2014); Tokyo 2015 workshop (eds. Hino, Stern, Springer 2020); SIMMAN 2008/2014; El Moctar, Shigunov & Zorn (2012) *Ship Technology Research* 59(3) 50–64 (DTC) |
| volume from a triangulated surface | divergence theorem on the closed polyhedral boundary; independent recomputation on a second code path | Mirtich (1996) *J. Graphics Tools* 1(2) 31–50 (*Fast and accurate computation of polyhedral mass properties*) |
| discretisation error must be shown, not assumed | grid-refinement study with observed order and Richardson extrapolation | Richardson (1911); Roache (1994) *J. Fluids Eng.* 116 405–413; Stern et al. (2001) *J. Fluids Eng.* 123 793–802; ASME V&V 20-2009 |
| "validated" has a defined meaning | verification (is the equation solved right?) is separated from validation (is the right equation solved?); a self-consistent closure is *verification only* | ASME V&V 20-2009; Stern et al. (2001) — hence the manifest's `basis` column: `closure` is never printed as an external validation |
| provenance and identity | FAIR principles; PROV-DM; content-addressed artefacts | Wilkinson et al. (2016) *Sci. Data* 3 160018; W3C PROV-DM (2013); RO-Crate (Soiland-Reyes et al. 2022, *Data Science* 5(2) 97–138) |
| solid-model integrity | a valid solid is closed, 2-manifold, consistently oriented, genus 0 for a hull shell | Mäntylä (1988) *An Introduction to Solid Modeling*; STL format spec (3D Systems, 1989) |
| ship geometry exchange | AP216 (ship moulded forms) / AP218 (ship structures) are the exchange formats an official hull would arrive in | ISO 10303-216 / ISO 10303-218 |

Two consequences that are *enforced in code* rather than promised in prose:

1. **A cited displacement over-determines the particulars.** If a record carries published
   `L, B, T, Cb` *and* a published displacement volume `∇`, the four numbers must satisfy
   `∇ ≈ Cb·L·B·T`. Round 2 makes that a gate (`G6b_reference_consistency`, tolerance 0.5 %);
   the band is far above the round-off of published digits (7 significant figures give ≲0.01 %)
   and far below any transcription error. **This is how a live defect was found**: the DTC
   record carried `B = 45.8 m` next to a cited `∇ = 173 467 m³`, i.e. 11.3 % inconsistent. The
   published quadruple `(355.0, 51.0, 14.5, 0.661)` reproduces the cited volume to 0.0005 % and
   gives `L/B = 6.96`, `B/T = 3.52` (post-Panamax container ship); it is the correction now in
   the manifest, and the mis-transcribed record is kept as a negative control.
2. **A closure reference is not evidence.** `∇ = Cb·L·B·T` can never disagree with the record,
   so `G6b` reports itself *vacuous* for such a record instead of scoring a pass, and `G7`
   labels the basis `closure (NOT an external validation)`.

## 3. Why this environment cannot close the provenance gap (honest limits)

* No egress: `https://simman2014.dk`, `https://tokyo2015.k.u-tokyo.ac.jp`,
  `https://www.uni-due.de/ship/dtc.php` and `raw.githubusercontent.com` are unreachable from the
  sandbox, so no workshop page and no official IGES/offset file can be re-read here. Every
  source-registry entry therefore carries `retrieval = offline_recollection`, `verified = false`.
* No official geometry on the box either: `/root/ship-data-official` contains only a DTMB 5415
  `.mat` database (not one of the six fleet hulls).
* Disk: the filesystem is 100 % full (≈28 KB free), which caps this round's artefact production
  at the coarse grid; the finer-grid numbers quoted in the report are the *previously recorded*
  production-grid run (`docs/gap_cad_t1_baseline_fleet.md`) plus this round's Richardson study.

## 4. Selection record (`oss_solutions`, `feature_slug = GAP-CAD-T1-BASELINE-FLEET`)

```json
{
  "feature_slug": "GAP-CAD-T1-BASELINE-FLEET",
  "round": 2,
  "mission": "m18641",
  "problem": "A hull geometry library cannot be trusted as a solver baseline unless every file is watertight+manifold, every ship's volume matches a published displacement, every digit has a named source, and placeholders can be replaced by gated geometry.",
  "candidates": [
    {"option": "adopt trimesh/manifold3d/OCCT as the validator", "verdict": "rejected",
     "reason": "adds a heavyweight dependency whose version drift becomes the audit trail; the four facts we need (boundary edges, non-manifold edges, winding consistency, Euler, signed volume) are ~120 lines and are already measured by an independent code path (tensorcad.hydrostatics)"},
    {"option": "trust the vendor STL and validate downstream in the solver", "verdict": "rejected",
     "reason": "the owner's defect (6 of 8 files are 304-triangle placeholders, one byte-identical pair) is invisible to a solver and contaminates every result computed from it"},
    {"option": "gate the fleet + anchor each digit to a source + deliver gated replacements", "verdict": "selected",
     "reason": "attacks all four acceptance criteria with machine-decidable gates, needs no network and no new dependency, and its failure modes are themselves gated by negative controls"}
  ],
  "selected": {
    "intake": "tensorcad/fleet/gates.py: G0 provenance, G0b source anchor, G1 not-placeholder, G2 no duplicate bytes, G3 watertight, G4 manifold+oriented, G5 single closed shell, G6 particulars coherent, G6b cited-volume vs Cb*L*B*T, G7 volume vs reference, G8 principal dimensions",
    "provenance": "tensorcad/fleet/data/source_registry.json (citation, URL, licence, locator, retrieval, verified) joined to the manifest and gated by G0b",
    "delivery": "tensorcad/fleet/replace.py — audit a library, pair defect files with manifest records, build through the same routes, gate, write, re-ingest, and emit ledger.json with before/after sha256",
    "artefact identity": "build.write_binary_stl + run._artifact_block: every admitted hull is written to manifests/<ship>.stl, sha256'd, re-read and re-gated (float32 round-trip tolerance 1e-6)",
    "verification": "per-ship verification_conclusion composed from the gates (no hand-written verdicts) + 5 negative controls that must be rejected"
  },
  "migration_plan": [
    {"step": 1, "action": "run the fleet gate over the ship-platform geometry library", "cmd": "python -m tensorcad.fleet.replace --source-dir <library> --out-dir <delivery>", "done": true},
    {"step": 2, "action": "ship the delivery directory back as a drop-in replacement (same file names, ledger.json alongside)", "done": "mechanism ready; needs a run against the real library, which is not in this worktree"},
    {"step": 3, "action": "route every new geometry through tensorcad.fleet before it reaches a solver", "done": "gate exists; wiring into the ship platform's /generate endpoint is the platform-side follow-up"},
    {"step": 4, "action": "replace offline-recollected anchors with fetched pages/IGES + sha256 once egress exists", "done": false, "blocker": "no outbound network; source registry records retrieval=offline_recollection, verified=false"}
  ],
  "reusable_for": ["TensorFEM", "TensorFVM", "TensorBEM", "TensorLBM"],
  "reuse_note": "The pattern (content-addressed artefacts + per-case manifest + must-fail controls + verification_conclusion) is domain-neutral; only data/source_registry.json and data/fleet_manifest.json are ship-specific."
}
```

## 5. What round 2 deliberately did *not* do

* It did not fabricate geometry for `KVLCC1` / `DTMB5415`: the replacement ledger reports them
  `WITHHELD` with the reason "no fleet-manifest record", because a hull invented to make a
  directory look complete is the same class of defect as the placeholder.
* It did not raise the volume tolerance to make the coarse grid pass: the 121 × 49 run is
  executed at a 2 % setting and its 1.19–1.23 % errors are printed as such, while the 1 %
  criterion is reported against the production grid and the Richardson study.
* It did not claim a fetched source for any digit. The registry's `verified=false` is the
  machine-readable form of that admission.