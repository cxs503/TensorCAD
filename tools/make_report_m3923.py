"""Build the PR #40 cure report payload (markdown + figure) for the platform endpoint."""
from __future__ import annotations

import base64
import json
import pathlib

ROOT = pathlib.Path(__file__).resolve().parents[1]
png = (ROOT / "artifacts/m3923/merge_verification.png").read_bytes()
d = json.loads((ROOT / "artifacts/m3923/merge_verification.json").read_text())
p, b, t = d["preservation"], d["m3169_route"], d["main_route_transom"]
i, c, f = d["main_route_inertia"], d["composition"], d["small_angle_cross_module_finding"]
r0 = f["rows"][0]

md = f"""## 1. Case

PR #40 (`agent/m3169` into `main`) in `TensorCAD` was `mergeable=false`. Both sides moved
after merge base `df31744`: `main` sits at `2ff80b2` (PR #39 landed, 12 commits ahead), the
branch at `fcf15ad` (1 commit) adding the large-angle righting-arm curve `GZ(phi)`. `main`
meanwhile added the transom-stern hull family, the waterplane second moments on the
`Hydrostatics` record, an asymmetric hull and a curvature rewrite.

Geometry used for the checks: wall-sided barge L=10 m, B=2 m, T=0.8 m, freeboard 0.4 m,
KG=0.55 m above keel; unit sphere at tessellation resolution 32; transom hull family L=10 m,
B=2 m, T=0.8 m, tau in {{0, 0.2, 0.4}}, breadth exponent 0.7, 41x13 grid.

## 2. Method

`git merge --no-ff origin/main` from `fcf15ad`. Two files conflicted (`CHANGELOG.md`,
`tests/test_public_surface.py`) and were resolved as unions: both sides' lines kept verbatim,
`main`'s `### Added` heading kept once. `tensorcad/__init__.py` auto-merged and was then
*verified*: both sides' `__all__` parsed and diffed against the merged list (23 names added by
the two sides, 0 missing). Preservation was re-derived from git plumbing by
`tools/verify_merge_m3169.py`: per-path existence, per-side *added-line* presence in the merged
file, and blob identity for one-sided paths. Capability liveness was measured against
independent references (wall-sided closed form, `GM sin(phi)` for a sphere, closed-form taper
integrals, `L B^3/12`), and the two routes were composed on one object: main's transom mesh
solved by m3169's heeled integrator.

## 3. Results

| quantity | value | reference |
|---|---|---|
| paths touched by either side | {p['paths_touched']} | - |
| files missing from merged tree | {len(p['files_missing_from_result'])} | 0 |
| added lines checked (both sides) | {p['added_lines_checked']} | - |
| added lines lost | {len(p['added_line_losses'])} | 0 |
| one-sided paths blob-identical | {p['files_touched_by_one_side']}/{p['files_touched_by_one_side']} | all |
| barge GZ, max abs error | {b['barge']['max_abs_error_m']:.2e} m | closed form, heels < {b['barge']['applied_validity_deg']:.1f} deg |
| sphere GZ, max abs error | {b['sphere']['max_abs_error_m']:.2e} m | GM sin(phi), GM = {b['sphere']['gm']} m |
| transom volume, max abs rel. error | {t['max_abs_volume_relative_error']*100:.2f} % | closed-form taper integrals |
| waterplane inertias, max rel. error | {i['max_abs_relative_error']:.1e} | L B^3/12, B L^3/12 |
| V record vs heeled (upright) | {c['max_volume_relative_difference']:.1e} | two independent integrations |
| KB record vs heeled (upright) | {c['max_kb_abs_difference_m']:.1e} m | same |
| displacement residual at every heel | {c['max_volume_residual']:.1e} | constant-displacement solve |
| kg lever, max deviation from -sin(phi) | {c['max_kg_lever_error_m']:.1e} m | exact linearity |

Tests, identical command on both trees (the two neural modules that cannot import here,
missing `safetensors`, excluded identically): `main` baseline **818 passed, 2 failed, 11
skipped (472 s)**; merged **859 passed, 2 failed, 11 skipped (823 s)**. The two failures are the
*same* two environment-related `tests/test_geometry_field.py` tests on both trees, so the merged
tree is **not below baseline** (+41 passes, identical failure set). Merge-critical modules:
`test_righting_arm + test_transom + test_public_surface` = 143 passed in 40 s.

![verification](figure)

## 4. Figures

Panel 1: barge `GZ(phi)` from the merged tree over the exact closed form (indistinguishable at
1e-14 m). Panel 2: `GZ` of transom hulls tau = 0, 0.2, 0.4 - main's geometry, m3169's solver.
Panel 3: error budget on a log axis with the 3 % and 15 % acceptance lines, plus the preservation
count.

## 5. Notes and limits

* **Finding left standing.** On a fine hull the small-angle slope from the heeled solver is
  `KB - KG` while the `Hydrostatics` record predicts `KB + I_T/nabla - KG`; the gap is *exactly*
  `BM` (parametric hull: {r0['slope_from_heeled_solver_m']:.5f} vs {r0['kb_plus_bm_minus_kg_record_m']:.5f},
  gap {abs(r0['absolute_difference_m']):.5f} m = BM {abs(r0['bm_m']):.5f} m). On the barge - whose waterline
  does not land on a tessellation ring - the solver matches the closed form to 1e-14, so the
  likely cause is the cap/waterline ring coincidence that m3169's own documentation already
  flags, **but that diagnosis is a hypothesis, not verified**. The merge modified neither module,
  so the effect predates the union: verdict WITHHELD, published with both numbers instead of a
  weakened assertion.
* `git push` of the cured branch is **refused by this sandbox's security policy**, so the head
  could not be published from here: old head `fcf15ad`, new head `1eaf438` (merge commit
  `2e5e934` plus the record commit `1eaf438`). The same commit is the head of `agent/m3923`, so
  the platform's automatic push carries it; publishing it on `agent/m3169` makes PR #40
  mergeable, and the merge into `main` remains the merge-governor's step.
* Accepted scope: the merge, its record and the verification tool. No solver physics was changed,
  so the finding above is reported rather than fixed.
"""

payload = {
    "title": ("PR #40 cure (m3923): semantic union of main and agent/m3169 - 0 lines lost, "
              "both capabilities verified"),
    "summary": ("PR #40 (agent/m3169) was mergeable=false against main 2ff80b2. Merged main into the "
                "branch (2 conflicted files resolved as unions), audited preservation with git plumbing "
                f"(47 paths, {p['added_lines_checked']} added lines, 0 lost, 43 one-sided paths "
                "blob-identical), and measured both sides on the merged tree: barge GZ vs closed form "
                "1.2e-14 m, transom family 0.88 %, waterplane inertias 1.7e-16, cross-route agreement "
                "6.9e-13. Tests 859 passed vs main baseline 818, same 2 pre-existing environment "
                "failures. New head 1eaf438 (old fcf15ad); the sandbox refused the push from here, so "
                "the platform push carries the branch."),
    "scenario_key": "tensorcad_pr40_merge_cure",
    "mission_id": 3923,
    "tags": ["TensorCAD", "merge-cure", "PR40", "righting-arm", "transom", "preservation-audit"],
    "markdown": md,
    "figures": [{
        "name": "merge_verification.png",
        "caption": ("Left: merged-tree barge GZ over the exact closed form. Middle: GZ of transom hulls "
                    "tau=0, 0.2, 0.4 (main's geometry, m3169's solver). Right: error budget against the "
                    "3 % / 15 % acceptance lines and the preservation count."),
        "png_base64": base64.b64encode(png).decode(),
    }],
}

text = json.dumps(payload)
(ROOT / "artifacts/m3923/report_payload.json").write_text(text)
(ROOT / "artifacts/m3923/report_body.md").write_text(md)
pathlib.Path("/tmp/report.json").write_text(text)
print("payload bytes:", len(text))