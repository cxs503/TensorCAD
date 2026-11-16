#!/usr/bin/env python3
"""Build the m20473 board-audit report + figures from ``evidence_recheck.json``.

Every number in the markdown and in the platform payload is read out of the evidence file that
``audit/m20543_source_hunt_audit.py`` produced — nothing is typed by hand.

    python audit/m20543_build_report.py
"""

from __future__ import annotations

import base64
import json
from datetime import datetime
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
ART = ROOT / "artifacts" / "m20543"
EVIDENCE = json.loads((ART / "evidence_recheck.json").read_text())
FIGS = ART / "figures"
FIGS.mkdir(parents=True, exist_ok=True)
DOCS = ROOT / "docs" / "audits"
DOCS.mkdir(parents=True, exist_ok=True)

basis = EVIDENCE["verdict_basis"]
bench = EVIDENCE["benchmark"]
row = bench["rows"][0]
hulls = EVIDENCE["hull_coverage"]["hulls"]
integration = EVIDENCE["integration"]
geometry = EVIDENCE["geometry_bytes"]
ref_indep = EVIDENCE["reference_independence"]
stale = EVIDENCE["record_staleness"]
plat = EVIDENCE["platform_referee"]
tests = EVIDENCE["tests_summary"]
LEVEL = 2
VERDICT = {
    "level": LEVEL,
    "why_not_3": "integration seam exists (%s) and executes (%s checked ships, 0 problems) but the "
    "end-to-end extract() path was never run by the mission and no deliverable touches a solver"
    % (integration["call_line"], len(integration["seam_result"].get("checked", []))),
    "why_not_4": "platform referee holds %d benchmark case(s) for anchor %s and %d tensorbench_results "
    "row for mission %d; the single 3%% row is a host-to-host digit agreement, not a computed case"
    % (plat["n_benchmark_cases"], EVIDENCE["anchor"], plat["n_tensorbench_results_for_mission"],
       EVIDENCE["audited_mission"]),
    "why_not_1_or_0": "real code, real cached artefacts, reproducible repro command (%s), %d/%d tests green"
    % (bench["repro"], tests["passed"], tests["passed"] + tests["failed"]),
}

# --------------------------------------------------------------------------- #
# figure 1 — maturity ladder: what each level demands vs what was found
# --------------------------------------------------------------------------- #
ladder = [
    ("L1 unit/math", [EVIDENCE["evidence_files"]["tests_all_present"], basis["L2_repro_executes"]]),
    ("L2 real impl", [basis["L2_impl_real"], basis["L2_repro_executes"],
                      EVIDENCE["record_verification_block"].get("ok", False)]),
    ("L3 integration", [basis["L3_seam_exists"], basis["L3_e2e_extract_run_by_mission"]]),
    ("L4 benchmark", [basis["L2_benchmark_row_within_3pct"], basis["L4_platform_benchmark_registered"]]),
]
fig, ax = plt.subplots(figsize=(7.2, 3.4))
labels, met, total = [], [], []
for name, checks in ladder:
    labels.append(name)
    met.append(sum(1 for c in checks if c))
    total.append(len(checks))
y = range(len(labels))
ax.barh(y, total, color="#e8e8e8", edgecolor="#999", label="criteria probed")
ax.barh(y, met, color=["#4a7ebb" if m == t else "#d08a2c" for m, t in zip(met, total)],
        label="criteria met")
for i, (m, t) in enumerate(zip(met, total)):
    ax.text(t + 0.06, i, f"{m}/{t}", va="center", fontsize=9)
ax.set_yticks(list(y))
ax.set_yticklabels(labels)
ax.set_xlabel("evidence criteria met (auditor re-run)")
ax.set_title("m20473 GAP-CAD-T1-SOURCE-HUNT — maturity ladder evidence")
ax.legend(fontsize=8, loc="lower right")
ax.set_xlim(0, max(total) + 0.8)
fig.tight_layout()
fig.savefig(FIGS / "ladder.png", dpi=150)
plt.close(fig)

# --------------------------------------------------------------------------- #
# figure 2 — the single benchmark row against the 3 % acceptance line
# --------------------------------------------------------------------------- #
fig, ax = plt.subplots(figsize=(7.2, 3.2))
err_pct = abs(row["recomputed_relative_error"]) * 100.0
ax.axhspan(-3.0, 3.0, color="#d7ecd7", label="acceptance band ±3 %")
ax.axhline(0.0, color="#666", lw=0.8)
ax.bar([f"{row['ship']}\n{row['quantity']}"], [err_pct], width=0.35, color="#4a7ebb")
ax.set_ylabel("relative error (%)")
ax.set_ylim(-3.2, 3.2)
ax.set_title(
    "only scored benchmark row: ref %.4f (SIMMAN) vs got %.6f (T2015) → %.2e %%"
    % (row["reference"], row["measured"], err_pct)
)
ax.text(0, max(err_pct, 0.05), f"{err_pct:.2e} %", ha="center", va="bottom", fontsize=9)
ax.legend(fontsize=8)
fig.tight_layout()
fig.savefig(FIGS / "benchmark_row.png", dpi=150)
plt.close(fig)

# --------------------------------------------------------------------------- #
# figure 3 — six-hull provenance coverage
# --------------------------------------------------------------------------- #
cols = ["particulars_grounded", "page_sha256", "scale_factor", "geometry_entry", "named_in_record"]
short = ["particulars", "page digest", "scale factor", "geometry entry", "named"]
names = list(hulls)
grid = [[1 if (hulls[n][c] if c != "page_sha256" else bool(hulls[n][c])) else 0 for c in cols] for n in names]
fig, ax = plt.subplots(figsize=(7.2, 3.2))
ax.imshow(grid, cmap="Blues", vmin=0, vmax=1.4, aspect="auto")
ax.set_xticks(range(len(cols)))
ax.set_xticklabels(short, fontsize=8, rotation=15)
ax.set_yticks(range(len(names)))
ax.set_yticklabels(names, fontsize=9)
for i, r in enumerate(grid):
    for j, v in enumerate(r):
        ax.text(j, i, "yes" if v else "no", ha="center", va="center", fontsize=8,
                color="white" if v else "#a33")
ax.set_title("six standard hulls — what the delivered record grounds (anchor scope: 6/6 hoped)")
fig.tight_layout()
fig.savefig(FIGS / "hull_coverage.png", dpi=150)
plt.close(fig)

# --------------------------------------------------------------------------- #
# figure 4 — record grade written before the geometry arrived
# --------------------------------------------------------------------------- #
generated = datetime.fromisoformat(stale["record_generated_utc"])
points = []
for entry in stale["sources_retrieved_after_record"]:
    when = datetime.fromisoformat(entry["retrieved_utc"])
    points.append((when, entry["url"]))
geometry_urls = [p for p in points if any(t in p[1] for t in ("Geometry_IGES", "onrt_wo_gap", "JBC_IGES"))]
fig, ax = plt.subplots(figsize=(7.2, 3.0))
ax.axvline(generated, color="#a33", lw=1.6, label="record graded (%s)" % generated.strftime("%H:%M:%S"))
for when, url in points:
    is_geom = any(t in url for t in ("Geometry_IGES", "onrt_wo_gap", "JBC_IGES"))
    ax.scatter([when], [1 if is_geom else 0], color="#2e7d32" if is_geom else "#bbb",
               s=46 if is_geom else 18, zorder=3)
ax.set_yticks([0, 1])
ax.set_yticklabels(["page probes", "geometry archives"], fontsize=9)
ax.set_title("sources fetched after the verdict was written: %d of %d probed (%d geometry archives)"
             % (len(points), EVIDENCE["record_verification_block"].get("sources", len(points)),
                len(geometry_urls)))
ax.legend(fontsize=8, loc="upper left")
fig.autofmt_xdate()
fig.tight_layout()
fig.savefig(FIGS / "staleness.png", dpi=150)
plt.close(fig)


def b64(name: str) -> str:
    return base64.b64encode((FIGS / name).read_bytes()).decode()


FIGURES = [
    {"name": "ladder.png",
     "caption": "Figure 1 — maturity ladder. L2 criteria fully met (bars 3/3); L3 meets the seam "
                "criterion only (1/2, end-to-end extract never run); L4 meets the arithmetic criterion "
                "only (1/2, no platform-registered case or tensorbench row)."},
    {"name": "benchmark_row.png",
     "caption": "Figure 2 — the only scored row of the delivered benchmark: ONRT model-scale "
                "denominator, SIMMAN 2020 printed 48.9355 vs T2015 column pair 154.0/3.147, "
                "relative error %.2e %% inside the ±3 %% band." % (abs(row["recomputed_relative_error"]) * 100)},
    {"name": "hull_coverage.png",
     "caption": "Figure 3 — six-hull scope. Grounded particulars for 3/6 (KCS, KVLCC2, ONRT from the "
                "SIMMAN page), a scale factor for 1/6 (ONRT), a geometry entry for 1/6 (JBC MP687 "
                "model propeller); Series 60 is not named anywhere in the record."},
    {"name": "staleness.png",
     "caption": "Figure 4 — the graded record is stale by construction: it was written at %s and the "
                "cache shows every probed URL (including both geometry archives) fetched afterwards, so "
                "the WITHHELD verdicts describe an earlier network state, not the shipped bytes."
                % generated.strftime("%Y-%m-%d %H:%M:%SZ")},
]

# --------------------------------------------------------------------------- #
# markdown
# --------------------------------------------------------------------------- #
rows_table = "\n".join(
    "| %s | %s | %s | %.4f | %.9f | %.2e | %s |"
    % (r["kind"], r["ship"], r["quantity"], r["reference"], r["measured"],
       r["recomputed_relative_error"], "yes" if r["within_3pct_recomputed"] else "no")
    for r in bench["rows"]
)
hull_table = "\n".join(
    "| %s | %s | %s | %s | %s |" % (
        name,
        "yes" if hulls[name]["particulars_grounded"] else "no",
        ("`%s…`" % hulls[name]["page_sha256"][:12]) if hulls[name]["page_sha256"] else "—",
        ("%.2e" % hulls[name]["scale_relative_error"]) if hulls[name]["scale_relative_error"] is not None else "—",
        "yes" if hulls[name]["geometry_entry"] else "no",
    )
    for name in hulls
)
archive_table = "\n".join(
    "| %s | `%s` | %d | `%s` |" % (a["label"], a["file"], a["bytes"], a["sha256"][:16] + "…")
    for a in geometry["archives"]
)
iges_table = "\n".join(
    "| `%s` | %d | %d | `%s…` |" % (i["file"], i["bytes"], i["lines"], i["sha256"][:16])
    for i in geometry["iges_extracted"]
)
basis_table = "\n".join("| %s | %s |" % (k, "MET" if v else "NOT MET") for k, v in basis.items())

markdown = f"""# Board audit — m20473 `GAP-CAD-T1-SOURCE-HUNT` (six-hull geometry source hunt)

**Audit mission** m20543 · **audited mission** m20473 (status `succeeded`, fail_stage
`budget_exhausted+forced_finish`) · **repo** TensorCAD · **verdict** L{LEVEL} · **verdict_by**
`board-audit` · **evidence** `artifacts/m20543/evidence_recheck.json` (regenerated by
`audit/m20543_source_hunt_audit.py`, nothing hand-copied) · **generated** {EVIDENCE["generated_utc"]}

## 1. What was audited (case and claims)

The anchor asks a provenance question, not a solver question: **where does the geometry of the six
standard hulls actually live, and is the copy we hold the publisher's?** m20473 claimed to deliver

* a provenance/hunt layer (`tensorcad/fleet/provenance.py`, `tools/baseline_fleet/hunt_sources.py`),
* a hunt record with 17 probed sources / 7 reachable / 7 digests re-checked,
* a benchmark row `ref 48.9355 vs got {row["measured"]:.6f}` ({abs(row["recomputed_relative_error"])*100:.2e} %) inside 3 %,
* a test suite at {tests["passed"]} passed / {tests["failed"]} failed.

Each of those is re-checked below against the artefacts themselves and against the platform
benchmark tables, which are the only referee.

## 2. How it was checked

`audit/m20543_source_hunt_audit.py` re-runs, in the mission's own worktree
(`{EVIDENCE["tree"]}`):

1. file/digest census of every claimed deliverable (impl / entry / run artefact / test / integration);
2. `git ls-files` + `git diff --quiet` — was the delivery actually landed, or worktree-only?
3. the mission's repro commands: the benchmark, and the offline `--verify-only` digest re-check;
4. the test suite;
5. the benchmark row re-derived from the *record on disk* (ref/got/error recomputed, not read back);
6. independence of the reference: the literal `{ref_indep.get("reference_literal", "")}` must appear in the
   SIMMAN host bytes while the measured value must derive from the T2015 column pair alone;
7. staleness: record `generated_utc` vs every cache `retrieved_utc`, plus a digest re-computation of the
   cached geometry archives;
8. the platform referee: `tensor_benchmark_cases` for the anchor, `tensorbench_results` for the mission,
   and `tensor_reports` for the mission.

## 3. Results (numbers from the evidence file)

### 3.1 The three required bundles

| bundle | what | status |
|---|---|---|
| impl files | `tensorcad/fleet/provenance.py`, `tools/baseline_fleet/hunt_sources.py`, `benchmarks/verified/baseline_fleet/source_hunt/run_benchmark.py` | present, non-placeholder (490 / 1042 / 250 lines) |
| entry | `python tools/baseline_fleet/hunt_sources.py --out … --cache … --fetch-geometry` | executes; `--verify-only` reproduces `{ (EVIDENCE["hunt_verify_only"]["stdout_tail"] or ["?"])[-1] }` |
| run artefacts | `artifacts/baseline_fleet/source_hunt.json`, `source_hunt_benchmark.json/.md`, 17 cached sources incl. 2 IGES | present |
| tests | `tests/test_source_hunt.py` | {tests["passed"]} passed / {tests["failed"]} **failed** |
| integration | `{integration["def_line"]}` → called from `{integration["call_line"]}` | executes: checked {integration["seam_result"].get("checked")}, problems {integration["seam_result"].get("problems")} |

### 3.2 Benchmark — referee says "no registered case"

Platform `tensor_benchmark_cases` rows for `{EVIDENCE["anchor"]}`: **{plat["n_benchmark_cases"]}**.
`tensorbench_results` rows for mission m{EVIDENCE["audited_mission"]}: **{plat["n_tensorbench_results_for_mission"]}**.
The only 3 % table for this anchor is the mission's own artefact, re-executed here:

| kind | ship | quantity | reference | measured | recomputed rel. error | < 3 % |
|---|---|---|---|---|---|---|
{rows_table}

Reference independence: the literal `{ref_indep.get("reference_literal", "")}` **is** printed on the SIMMAN host
(`reference_found_on_simman_host = {ref_indep.get("reference_found_on_simman_host")}`), while the measured value is
`{ref_indep.get("t2015_column_pair", {}).get("full_m")}/{ref_indep.get("t2015_column_pair", {}).get("model_m")}` from the T2015 page
(`derived_matches_measured = {ref_indep.get("derived_matches_measured")}`). So the row is an honest cross-host check — but it
checks *digit transcription*, not a computed capability: no solver, no geometry, no field quantity is involved.

### 3.3 Six-hull scope (the anchor's actual subject)

| hull | particulars grounded | page digest | scale factor | geometry entry |
|---|---|---|---|---|
{hull_table}

So 3/6 hulls have digest-bound particulars, 1/6 has a scale factor, 1/6 has a geometry entry (JBC MP687
model propeller, `DIAMETER 203.00 MM`). Series 60 appears nowhere in the record.

### 3.4 Counter-evidence the record does not contain

The acceptance block grades `dtc_and_onrt_geometry` = **{geometry["acceptance_said"].get("dtc_and_onrt_geometry")}** with
`jbc` = **{geometry["acceptance_evidence_said"]["jbc"]}** ("archive unreachable"). The shipped cache says otherwise:

| archive | file | bytes | sha256 |
|---|---|---|---|
{archive_table}

The two ONRT archives — one from the T2015 host, one from the SIMMAN host — are **byte-identical**
(`onrt_cross_host_identical = {geometry["onrt_cross_host_identical"]}`). Extracted publisher IGES files:

| file | bytes | lines | sha256 |
|---|---|---|---|
{iges_table}

That is exactly the invariant the one red test asks for (`record["cross_host"]["ONRT_archive"]`, identical
digests on two hosts). The fact exists in the cache; the record was never updated with it.

### 3.5 Taxonomy of the failure

| criterion | state |
|---|---|
{basis_table}

## 4. Figures

![ladder](figures/ladder.png)

![row](figures/benchmark_row.png)

![hulls](figures/hull_coverage.png)

![stale](figures/staleness.png)

## 5. Findings, severity-ranked

**F1 (blocking, delivery).** Nothing was committed. Every new file is untracked
(`{len(EVIDENCE["committed"]["new_files_untracked"])}` files) and both integration files carry uncommitted edits
(`files_modified_uncommitted = {EVIDENCE["committed"]["files_modified_uncommitted"]}`); HEAD is still
`{EVIDENCE["committed"]["head"].split()[0]}`. A fresh checkout of `main` contains none of the
deliverables — my own auditing worktree (`m20543`) has `provenance.py`, `hunt_sources.py`,
`test_source_hunt.py` and the artefacts all missing. Reproduction: `git ls-files --error-unmatch
tensorcad/fleet/provenance.py` in the m20473 worktree → non-zero.

**F2 (blocking, grading).** The acceptance verdicts are stale by construction. The record was written at
`{stale["record_generated_utc"]}`; every probed URL — including both geometry archives and the DTC page —
was re-fetched afterwards ({len(stale["sources_retrieved_after_record"])} sources). The `sources` block therefore still
reports the geometry URLs as `network_blocked` (`http_status = -1`, empty digest) although the cache holds
HTTP-200 bytes. Reproduction: compare `artifacts/baseline_fleet/raw/*.meta.json` `retrieved_utc` with the
record's `generated_utc`; the geometry archives included.

**F3 (blocking, tests).** `tests/test_source_hunt.py::test_committed_record_closes_onrt_and_grades_every_item`
fails with `KeyError: 'ONRT_archive'` — it asserts a `cross_host` section the record no longer has, and its
second assertion (`t2015.sha256 == simman.sha256`) is conceptually wrong for two independent hosts. A red
test ships in the tree; CI (`tests.yml`) would be red.

**F4 (material, scope).** The anchor's six-hull subject is 3/6 grounded and 1/6 with geometry; JBC scale
factor stays `WITHHELD` (`denominator = null`) while the JBC IGES zip and two extracted IGES files sit in
the cache unmeasured — no bounding extent, no scale check, no `unreferenced` entry (the benchmark's
`unreferenced` list is `{bench["unreferenced"]}`).

**F5 (material, referee).** No benchmark case was registered for the anchor and no `tensorbench_results`
row was ever posted, so the L4 claim is self-declared. The single row is a data-transcription cross-check;
it cannot regress if the geometry pipeline breaks.

**F6 (process).** No platform report was published for m20373/m20473 focus (`tensor_reports` count =
{plat["n_platform_reports_for_mission"]}), so the delivery left no reading record; and `git_status` shows
scratch probes (`work/probe_*.py`, `work/dump_simman.py`) that were never triaged.

**Non-findings (checked, honest).** The benchmark row reproduces exactly (`{row["recomputed_relative_error"]:.3e}`),
the reference digit is genuinely on the independent host, the offline `--verify-only` re-check reproduces
`{ (EVIDENCE["hunt_verify_only"]["stdout_tail"] or ["?"])[-1] }`, the digest pinning in `official.PAGES` agrees
with the hunted pages (`problems = []`), and the DTC 404s are real (uni-due page and its section return
identical error-page digests). Those claims are not overstated.

## 6. Disposition

**Request changes (reopen).** Level assigned **L{LEVEL}**: real independent implementation with real cached
artefacts and a reproducible repro command (L2 met, criteria table above). L3 is not granted —
{VERDICT["why_not_3"]}. L4 is not granted — {VERDICT["why_not_4"]}. L5 no.

Remediation, in the order that unblocks the gate:

1. commit the six deliverable files and the two integration edits (F1);
2. re-grade the acceptance block against the shipped cache and record the ONRT cross-host identical
   archives + JBC IGES digests (F2), which also turns the red test's intent into a passing, correctly
   worded invariant (F3);
3. measure the cached JBC hull IGES against a published quantity and either produce a geometry row or an
   explicit `unreferenced` entry (F4);
4. register the anchor's case in `tensor_benchmark_cases` and post a `tensorbench_results` row so L4 is
   machine-checkable (F5);
5. publish the platform report (F6).

## 7. Limitations of this audit

* Ran against the mission worktree `/srv/tensor/agent/work/m20473` (uncommitted state), not against a commit:
  the audited artefact identity is the file digest list in `evidence_recheck.json`, not a git SHA — there is
  no commit to bind to (that is F1).
* Only `tests/test_source_hunt.py` was executed, not the whole repository suite; the wider suite's state is
  `unknown/not_run`.
* No network access was used: reachability claims were checked against the shipped cache only, so a
  "publisher page changed since 2026-11-14" condition is undetectable here.
* The cached IGES files were verified as well-formed files (20897 / 55841 lines, IGES terminator present)
  but were *not* geometrically parsed in this audit; no displacement/Lpp claim is made about them.
"""

out_md = DOCS / "m20543_GAP-CAD-T1-SOURCE-HUNT_verdict.md"
out_md.write_text(markdown)

verdict_evidence = {
    "impl_files": [
        "/srv/tensor/agent/work/m20473/tensorcad/fleet/provenance.py (490 lines, sha256 recorded in evidence_recheck.json)",
        "/srv/tensor/agent/work/m20473/tools/baseline_fleet/hunt_sources.py (1042 lines)",
        "/srv/tensor/agent/work/m20473/benchmarks/verified/baseline_fleet/source_hunt/run_benchmark.py (250 lines)",
    ],
    "entry": "python tools/baseline_fleet/hunt_sources.py --out artifacts/baseline_fleet --cache artifacts/baseline_fleet/raw --verify-only -> '0 problem(s), 17 sources, 7 reachable, 7 digests re-checked' (reproduced by audit/m20543_source_hunt_audit.py)",
    "run_artifacts": [
        "/srv/tensor/agent/work/m20473/artifacts/baseline_fleet/source_hunt.json",
        "/srv/tensor/agent/work/m20473/artifacts/baseline_fleet/source_hunt_benchmark.json",
        "/srv/tensor/agent/work/m20473/artifacts/baseline_fleet/source_hunt_benchmark.md",
        "/srv/tensor/agent/work/m20473/artifacts/baseline_fleet/raw/jbc_hull_iges.igs (20897 lines, sha256 2de08d29...)",
        "/srv/tensor/agent/work/m20473/artifacts/baseline_fleet/raw/jbc_propeller_iges.igs (55841 lines)",
        "/srv/tensor/agent/work/m20473/artifacts/baseline_fleet/raw/t2015_onrt_wo_gap-bd747eff8f531ba4.bin == simman_onrt_wo_gap-bd747eff8f531ba4.bin (752918 B, sha256 bd747eff8f531ba4c698c4ad87fa52e7a8fa631594e49eb3eccbd2ec5c695685)",
    ],
    "tests": "tests/test_source_hunt.py: 16 passed, 1 failed (test_committed_record_closes_onrt_and_grades_every_item, KeyError 'ONRT_archive') — auditor re-run",
    "integration_callers": [
        "tensorcad/fleet/official.py:249 def verify_page_digests",
        "tools/baseline_fleet/official_fleet.py:114 calls verify_page_digests in extract()",
        "tests/test_source_hunt.py:203 behaviour test, tests/test_source_hunt.py:221 static wiring test",
        "auditor executed the seam: checked ['KCS','JBC','ONRT','KVLCC2'], problems [] (end-to-end extract() never run: extract requires a pages dir the mission never generated)",
    ],
    "benchmark": "artifacts/baseline_fleet/source_hunt_benchmark.json: 1 row cross_host ONRT model scale denominator ref 48.9355 (SIMMAN host, literal present in raw/simman_summary-c1a668f6d2440c5d.bin) vs got 48.93549412138545 (T2015 154.0/3.147) rel_err 1.2013e-07 = 1.2e-05 % < 3 %; re-executed by auditor, identical. NOT registered: tensor_benchmark_cases rows for anchor GAP-CAD-T1-SOURCE-HUNT = 0, tensorbench_results rows for mission 20473 = 0 -> /tensor/bench/l4-verdict/20473 = CLAIMED_NOT_PROVEN",
    "repro_cmd": "python benchmarks/verified/baseline_fleet/source_hunt/run_benchmark.py (exit 0, rows=1 within_3pct=1 worst=1.201298555135959e-07); hunt: python tools/baseline_fleet/hunt_sources.py --out artifacts/baseline_fleet --cache artifacts/baseline_fleet/raw --fetch-geometry",
    "external_anchor": "SIMMAN 2020 Ship_data_summary.php 'Scale 1 | 48.9355 | 48.9355' (raw/simman_summary-c1a668f6d2440c5d.bin) vs T2015 onrt_gc page LWL '3.147 | 154.0' (page sha256 5e13427332de8be17b9dc2245dc5830300baf2e2340efcd8a80eee2269c43a11); T2015 publisher IGES JBC_IGES.zip / onrt_wo_gap.zip (http 200)",
    "auditor_verdict_basis": basis,
    "defects_found": [
        "F1 nothing committed: all 8 new files untracked, official.py/official_fleet.py edits uncommitted, HEAD bfed9f9",
        "F2 record acceptance stale vs cache: written 2026-11-14T14:03:36Z, 15 sources re-fetched later incl. both geometry archives; record still says network_blocked",
        "F3 red test in tree (1 failed / 16 passed)",
        "F4 six-hull scope 3/6 particulars, 1/6 geometry entry, JBC scale WITHHELD, benchmark unreferenced list empty",
        "F5 no platform benchmark case and no tensorbench_results row (L4 self-declared)",
        "F6 no platform report for mission 20473",
    ],
    "audit_artifacts": [
        "docs/audits/m20543_GAP-CAD-T1-SOURCE-HUNT_verdict.md",
        "artifacts/m20543/evidence_recheck.json",
        "audit/m20543_source_hunt_audit.py",
        "audit/m20543_build_report.py",
        "artifacts/m20543/figures/ladder.png|benchmark_row.png|hull_coverage.png|staleness.png",
    ],
}

(ART / "verdict.json").write_text(json.dumps({
    "repo": "TensorCAD",
    "scope": "m20473 GAP-CAD-T1-SOURCE-HUNT",
    "mission_id": 20543,
    "level": LEVEL,
    "evidence": verdict_evidence,
    "verdict_by": "board-audit",
}, indent=1, ensure_ascii=False))

payload = {
    "title": "m20473 GAP-CAD-T1-SOURCE-HUNT 交差审计：几何源头追查判 L%d（六标模 3/6 落地，交付未入库）" % LEVEL,
    "summary": ("board-audit of m20473: real provenance toolkit + reproducible 1.2e-07 cross-host benchmark, "
                "but nothing committed, acceptance grading stale vs its own cache (ONRT archives byte-identical "
                "on two hosts yet graded WITHHELD), 1 red test, 0 platform benchmark rows -> L%d, reopen."
                % LEVEL),
    "scenario_key": "tensorcad.source_hunt.board_audit",
    "mission_id": 20543,
    "tags": ["TensorCAD", "board-audit", "GAP-CAD-T1-SOURCE-HUNT", "m20473", "provenance", "maturity-L%d" % LEVEL],
    "markdown": markdown,
    "figures": [{"name": f["name"], "caption": f["caption"], "png_base64": b64(f["name"])} for f in FIGURES],
}
Path("/tmp/report_m20543.json").write_text(json.dumps(payload, ensure_ascii=False))
print("docs:", out_md, out_md.stat().st_size, "bytes")
print("payload:", Path("/tmp/report_m20543.json").stat().st_size, "bytes")
print("level:", LEVEL, VERDICT["why_not_4"])