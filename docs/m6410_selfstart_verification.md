# m6410 selfstart: full-suite audit, artifact regeneration, intersect defect

**Date**: 2026-11-08 · **Branch**: agent/m6410 · **Scope**: whole-repo test audit +
defect evidence. No solver code was modified in this mission.

## 1. Case

The suite had grown to 60 test files / 1533 collected cases since the 2026-10-24
audit baseline (262 cases). A full run on this worktree (loong64, 128 cores,
load ~125) measured:

| suite | collected | passed | skipped | failed | wall time |
|---|---|---|---|---|---|
| `pytest tests/ -q` | 1533 | 1513 | 15 | **5** | 614.87 s |

## 2. What was done

1. **Three failures were missing committed artifacts**, filtered out of the
   merge pipeline ("N 个文件被门禁滤除" in the merge commits). Regenerated with
   the committed generators, no code change:
   * `artifacts/m4873/box_gz_areas.json` (+ `.png`) via `tools/verify_box_gz_areas.py`
     -- fixes `test_box_gz_area.py::test_the_committed_evidence_still_makes_its_claims`
     and `test_the_artifact_is_what_the_code_computes_today`.
   * `artifacts/m5479_geometry_audit/gates_verified.json` via
     `tools/generate_gates_m5479.py` -- fixes
     `test_geometry_audit_m5479.py::TestC5GatesArtifact::test_gates_file_exists_and_valid`.
   * Re-run after regeneration: both files green (see report).
2. **Two failures are a real defect in `tensorcad/intersect`** (merged by m6300,
   ~24 h before this mission). `TestExactDifferential` requires the two exact
   routes to agree with each other; they do not:
   * `test_three_exact_routes_agree_on_dyadic_pairs`: **275 / 800** random dyadic
     triangle pairs disagree between the line route and the clip3d route.
     First disagreeing pair (reproduction):
     `T1 = ((-0.5,1.0,-0.5),(0.0,0.5,0.0),(-1.0,-1.0,0.5))`,
     `T2 = ((0.5,0.5,0.5),(-0.5,0.5,-0.5),(-1.0,1.0,0.5))`.
   * `test_thirds_lattice_matches_halfspace_clipping`: **12** thirds-lattice
     pairs disagree with the independent half-space clipping oracle.
   * Dyadic coordinates are exact in binary, so this is a logic defect, not
     an architecture/rounding effect. It was **not** fixed in this mission
     (root-causing two 3-D exact routes is beyond a single lightweight
     session); a mission was dispatched instead. Honest FAIL, not covered.
3. **NURBS closed-surface verification gap** (verification debt, not closed):
   `docs/nurbs-surface.md` anchors an open tube and a capped cylinder only.
   A closed NURBS surface (pole rings welding to single vertices, degenerate
   pole quads, winding survival) has no committed anchor. A first attempt
   (rational B-sphere, radius anchor `||x|| = R`) was started and withdrawn
   because the control net could not be verified within this session's budget;
   the withdrawn file was removed rather than committed red. Follow-up
   candidate for the next selfstart.

## 3. Results

| check | before | after |
|---|---|---|
| `test_box_gz_area.py` | 2 failed | 0 failed |
| `test_geometry_audit_m5479.py` | 1 failed | 0 failed |
| `test_intersect.py::TestExactDifferential` | 2 failed | 2 failed (defect, mission dispatched) |
| full suite | 5 failed | 2 failed (the intersect defect only) |

## 4. Notes / limitations

* The full-suite number above is from this worktree on this architecture
  (loong64); the intersect defect is coordinate-exact and architecture-
  independent by construction (dyadic inputs).
* No physics quantity is computed here, so there is no literature benchmark
  to compare against; the "benchmark" is the suite itself plus the committed
  generators.
