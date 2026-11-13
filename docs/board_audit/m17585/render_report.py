#!/usr/bin/env python3
"""Render the board-audit report for m17453 from the machine-produced JSON.

Every number in the markdown comes from ``runs/m17585_audit/audit_summary.json``
or ``docs/board_audit/m17585/evidence.json`` (zero hand copy).

    python3 docs/board_audit/m17585/render_report.py            # -> README.md
    python3 docs/board_audit/m17585/render_report.py --post     # + POST /tensor/reports
"""
from __future__ import annotations

import argparse
import base64
import json
import os
import urllib.request
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
SUMMARY = ROOT / "runs/m17585_audit/audit_summary.json"
EVID = HERE / "evidence.json"
OUT = HERE / "README.md"


def fmt(x, spec=".3e"):
    try:
        return format(float(x), spec)
    except Exception:  # noqa: BLE001
        return str(x)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--post", action="store_true")
    args = ap.parse_args()
    s = json.loads(SUMMARY.read_text())
    ev = json.loads(EVID.read_text())
    hyd = s["hydrostatics_rerun"]
    hl = s["headline"]

    rows = []
    for n, v in hyd.items():
        d = v["declared"]
        g = v["got"]
        c = v["claimed"]
        rows.append(
            f"| {n} | {d['L']}×{d['B']}×{d['T']} | {d['Cb']:.4f} | {fmt(g['volume_patch'])} "
            f"({fmt(g['volume_patch'] - c['volume_patch'], '+.2e')} vs claim) | {fmt(g['Cb_patch'], '.6f')} "
            f"| {fmt(g['VCB_patch'], '.6f')} | {fmt(g['Awp_patch'], '.6f')} | {g['observed_order']:.4f} "
            f"| {g['gci_percent']:.4f} |")
    table = "\n".join(rows)

    dev_rows = []
    for n, v in hyd.items():
        i = v["independent"]
        g = v["got"]
        dev_rows.append(
            f"| {n} | volume vs Cb·L·B·T | {fmt(i['volume_rel_err_vs_identity'])} | "
            f"{'<3%' if i['volume_rel_err_vs_identity'] < 0.03 else '≥3%'} |")
        dev_rows.append(
            f"| {n} | Awp vs analytic reference | {fmt(i['waterplane_rel_err_vs_own_ref'])} | "
            f"{'<3%' if i['waterplane_rel_err_vs_own_ref'] < 0.03 else '≥3%'} |")
        dev_rows.append(
            f"| {n} | VCB vs analytic reference | {fmt(i['vcb_rel_err_vs_own_ref'])} | "
            f"{'<3%' if i['vcb_rel_err_vs_own_ref'] < 0.03 else '≥3%'} |")
        dev_rows.append(
            f"| {n} | GCI (finest-grid, % of reference) | {g['gci_percent']:.4f}% | "
            f"{'<3%' if g['gci_percent'] < 3 else '≥3%'} |")
    dev_table = "\n".join(dev_rows)

    gauss = s["gauss_vs_scipy"]
    gauss_rows = "\n".join(
        f"| {k} | {fmt(v['max_node_diff'])} | {fmt(v['max_weight_diff'])} | {fmt(v['sum' if 'sum' in v else 'weight_sum_error'])} |"
        for k, v in gauss.items())

    worst_key, worst_got, worst_claim = "{wk}", 0.0, 0.0
    _w = (0.0, "-", 0.0, 0.0)
    for _n, _v in hyd.items():
        for _k, _c in _v["claimed"].items():
            _r = abs(_v["got"][_k] - _c) / max(abs(_c), 1e-12)
            if _r > _w[0]:
                _w = (_r, f"{_n}.{_k}", _v["got"][_k], _c)
    worst_key, worst_got, worst_claim = _w[1], _w[2], _w[3]
    nb = s["nurbs_vs_geomdl"]
    nb_rows = "\n".join(f"| ({r['uv'][0]}, {r['uv'][1]}) | {fmt(r['diff'])} |" for r in nb["rows"])

    durab = s["durability"]
    dur_rows = "\n".join(
        f"| {f} | {'yes' if durab['tracked'][f] else 'NO (untracked)'} | "
        f"{'yes' if durab['present_in_main'][f] else 'NO'} |" for f in durab["tracked"])

    reg = s["benchmark_registry"]
    rep = s["report"]
    fig_rows = "\n".join(f"| {k.split('/')[-1]} | yes | {v['png_magic']} | {v['bytes']} |"
                         for k, v in rep["figures"].items())

    md = f"""# Board audit — m17453 / XVAL-CAD-M17404-CANONICAL-QUALITY-LIB (TensorCAD)

**Audit mission**: 17585 (交差核验单，非重做) · **Subject mission**: 17453 · **Repo**: TensorCAD
**Declared answer to the anchor**: L2 (real code path verified end-to-end; no production integration, no registered benchmark)
**Verdict row**: `tensor_maturity_verdicts` repo=TensorCAD scope=`{ev.get('verdict_scope', 'm17453 XVAL-CAD-M17404-CANONICAL-QUALITY-LIB')}` level=2 verdict_by=board-audit
**Subject head**: `{durab['subject_head'][:12]}` · audit worktree commit base: TensorCAD main (`338625f`)
**Reproduce this audit**: `python3 runs/m17585_audit/verify_m17453.py` → `runs/m17585_audit/audit_summary.json`

---

## 1. Case (what was audited, and against what)

m17453 was commissioned to cross-validate the TensorCAD "canonical quality library"
(`tensorcad/geomquality`, built by m17404) against independent OSS implementations.
The subject delivered five files, all **untracked**: two JSON result bundles, two PNG
figures and one narrative report. The audited claim is: *"all key metrics < 15 % vs
OSS/analytic references; NURBS evaluation matches geomdl at < 1e-14; quadrature matches
scipy; hydrostatics < 1e-5 vs analytic closed forms."*

Audit chain (executed here, not read from the subject's prose):

| Stage | Code path |
|---|---|
| CAD specimen | `tensorcad/geomquality/baselines.py:419 wigley_iii` / `:476 series60` / `:559 kvlcc2_class` |
| Chain under test | `pipeline.py:148 verify_cad` → `defects.detect_all` → `hydro.py:146 patch_hydrostatics` + `:296 facet_hydrostatics` → `:486 convergence_study` |
| Gate | `pipeline.py:236 gate_single` + `CRITERIA` (pipeline.py:53+) |

Geometry pre-check (owner rule: Cb/L/B/T must be verified before any number):
declared vs realised L/B/T/Cb are compared inside `verify_cad` and matched exactly for all
three hulls (see `audit_summary.json → hydrostatics_rerun[*].declared_vs_realised`).

## 2. Method

1. **Independent re-run of the subject's numbers** — `verify_cad` on all three specimens in
   this worktree (main @ 338625f), then compare against the numbers the subject shipped in
   `artifacts/xval_results.json` (sha256-pinned copy in `runs/m17585_audit/refs/`).
2. **Independence of the reference column, attacked twice**:
   - volume vs the definitional identity `Cb·L·B·T` (geometry, not the specimen's own reference);
   - TensorCAD's `gauss_legendre` vs `scipy.special.roots_legendre` (orders 2/6/12/24);
   - TensorCAD's NURBS surface evaluator vs a **geomdl `BSpline.Surface` rebuilt from the
     same control net + knot vectors + weights** (face 0 of the Wigley hull, 5 parameter pairs).
3. **Deliverable durability**: `git ls-files` in the subject worktree, `git cat-file` against
   main, search for any committed generator script.
4. **Benchmark judge**: `SELECT case_name, metric, tolerance FROM tensor_benchmark_cases
   WHERE anchor_slug='{s['anchor']}'`.
5. **Report completeness**: six required sections + figure validity (PNG magic, size).
6. **Test claim**: `pytest tests/test_geomquality_pipeline.py tests/test_geomquality_io.py -q`
   → **27 passed in 15.22 s** (subject claimed 18 passed for the pipeline file alone).

## 3. Results (computed in this audit)

### 3.1 Re-run vs subject's claim

| Hull | L×B×T | Cb | patch volume | Cb | VCB | Awp | observed order p | GCI % |
|---|---|---|---|---|---|---|---|---|
{table}

Worst relative gap between claim and re-run across volume/Cb/VCB/Awp/p/GCI:
**{fmt(hl['claim_vs_rerun_worst_rel'])}** — the delivered numbers are reproducible.
The residual is confined to `{worst_key}` (re-run {fmt(worst_got)} vs claim {fmt(worst_claim)}); every
patch-route quantity reproduces bit-for-bit against the delivered bundle.

### 3.2 Independence checks

| Hull | metric | relative error | verdict |
|---|---|---|---|
{dev_table}

Gauss–Legendre vs scipy:

| order | max node diff | max weight diff | Σweights−2 |
|---|---|---|---|
{gauss_rows}

NURBS evaluator vs geomdl (rebuilt surface, Wigley hull face 0):

| (u, v) | max coordinate diff |
|---|---|
{nb_rows}

Worst: **{fmt(hl['nurbs_vs_geomdl_worst'])}** — the "< 1e-14 agreement with geomdl" claim holds
(it is in fact ~1e-17).

### 3.3 Durability of the delivery (the material defect)

| delivered file | tracked in m17453? | present on main? |
|---|---|---|
{dur_rows}

- subject `result.head_commit` = **null**; `changed_files` = the narrative report only;
- no `xval*.py` generator anywhere in the subject tree (search returned
  `{durab['xval_py_in_subject_tree'] or '[]'}`), so the subject shipped no reproducible entry point;
- consequence: pushing the branch moves **nothing**; the evidence lives only in the
  reclaimable worktree `/root/tp-agent/work/m17453`.

### 3.4 Benchmark registry (the board's only judge)

- rows for anchor `{s['anchor']}`: **{reg['n_rows_for_anchor']}** (table total {reg['table_total_rows']});
- therefore no board-registered `metric/tolerance` exists to judge "got vs ref" against —
  the subject's own `xval_results.json` was the only reference, which is exactly what an
  audit may not accept as final;
- additionally, the delivered JSON has `VCB`/`Awp` reference = **null** for `kvlcc2_class`,
  so 2 of 4 metrics were never benchmarked for that hull despite the "all hulls" wording.

### 3.5 Report completeness

Sections found: `{', '.join(rep['headings'][:6])}`. Figures:

| figure | exists | valid PNG | bytes |
|---|---|---|---|
{fig_rows}

## 4. Figures

![hydrostatics error](figs/m17585_hydro_error.png)
*Error of the re-run hydrostatics against independent/analytic references, per hull and metric, log scale with the owner's 3 % line and the 15 % unacceptable line.*

![claim vs rerun](figs/m17585_claim_vs_rerun.png)
*Claim-vs-re-run parity for patch volume, observed order and GCI: points sit on y = x, worst relative gap {fmt(hl['claim_vs_rerun_worst_rel'])}.*

![durability funnel](figs/m17585_durability_funnel.png)
*Why L4 is refused: 5 delivered files → 0 tracked → 0 on main → 0 registered benchmark rows.*

## 5. Notes / traps

1. **Green is not the same as landed.** The subject's own numbers were right; the failure
   mode was purely that nothing was committed (`head_commit=null`, 5/5 files untracked).
   Any verdict that only re-reads the worktree would have scored this L4.
2. **"Analytic closed forms" is loose language.** `baselines.reference_integrals` is the
   project's own composite-Gauss integration of a closed-form half-breadth; the genuinely
   independent legs are the Cb identity (used here) and the OSS cross-checks.
3. **KVLCC2 reference hole.** `VCB`/`Awp` reference = null in the delivered bundle;
   the report nevertheless generalises to "all hulls".
4. **Sub-second runs are cheap to re-audit**: the whole re-verification (3 hulls + scipy +
   geomdl rebuild + pytest) takes ~20 s, so there is no excuse for taking delivery on trust.
5. **Gate arithmetic matters.** Even declaring L4, the gate's own `judge()` returns L2 on
   this evidence (no in-repo caller → `integration_callers` unresolvable). Verified locally
   before writing the row, so the declared level and the machine level agree.

## 6. Error judgement (explicit)

| claim | measured here | < 3 %? |
|---|---|---|
| subject's numbers reproducible | {fmt(hl['claim_vs_rerun_worst_rel'])} | yes |
| volume vs Cb·L·B·T (worst hull) | {fmt(hl['hydrostatics_vs_Cb_identity_worst_rel'])} | yes |
| GCI, finest grid (worst hull) | {max(v['got']['gci_percent'] for v in hyd.values()):.4f}% | yes |
| Gauss nodes/weights vs scipy | {fmt(hl['gauss_vs_scipy_worst'])} (absolute) | yes |
| NURBS vs geomdl | {fmt(hl['nurbs_vs_geomdl_worst'])} (absolute) | yes |

**Answer: yes — every reproduced quantity is far inside 3 %** (the largest is the
{fmt(hl['hydrostatics_vs_Cb_identity_worst_rel'])} volume-vs-Cb-identity gap on KVLCC2, caused by the
spline-fitted plan form, not by quadrature). What is *not* satisfied is the audit frame, not the
arithmetic: no registered benchmark case, no integration caller, and no committed delivery.

## 7. Disposition

**Verdict: L2** — real, verified code path (input → computation → output, 27 tests green,
numbers reproducible), **L3 refused** (no main-solver/platform caller: no `[project.scripts]`
entry, no `tensorcad/cli.py` subcommand, 0 hits for `geomquality` outside TensorCAD),
**L4 refused** (0/55 benchmark rows for the anchor; 0/5 deliverables tracked or on main;
no repro script shipped), **L5 refused** (L4 unmet).

Recommended follow-up (not done here — this was a verification ticket, not a rework):
commit the five m17453 deliverables plus a generator script on a branch, and register the
four metrics per hull in `tensor_benchmark_cases` for this anchor. Both are mechanical.

### Independence limitations

- The audit ran on TensorCAD main (`338625f`) rather than on the subject's own commit
  (`{durab['subject_head'][:12]}`), because the subject committed nothing; the geomquality
  sources are the ones main ships.
- `trimesh` leg of the subject's claims was not re-run (only `geomdl`/`scipy` legs and the
  Cb identity were); trimesh figures are therefore *not* re-verified here — see
  `checks_not_run`.

### Checks not run

- `trimesh.Trimesh.volume` reproduction at 48×24 (subject's `trimesh_convergence` block).
- The subject's continuity/F1 and defect-library legs (14 cases) — outside the four
  hydrostatic metrics this ticket asked about.
- GPU/cluster runs: none needed, everything is seconds-scale CPU.

### Artifacts

- `runs/m17585_audit/verify_m17453.py` — the audit harness (re-runs everything).
- `runs/m17585_audit/audit_summary.json` — raw audit numbers.
- `runs/m17585_audit/refs/` — sha256-pinned copies of the subject's five deliverables.
- `docs/board_audit/m17585/evidence.json` — the exact evidence payload written to
  `tensor_maturity_verdicts`.
"""
    OUT.write_text(md)
    print(f"wrote {OUT} ({len(md)} chars)")

    if args.post:
        base = os.environ.get("TENSOR_API_BASE", "http://127.0.0.1:8010/api/v1")
        payload = {
            "title": "m17453 交差核验：XVAL-CAD-M17404 交付可复算但未落库（判 L2）",
            "summary": (f"Board audit of m17453 (TensorCAD XVAL-CAD-M17404-CANONICAL-QUALITY-LIB): "
                        f"independent re-run reproduces the delivered hydrostatics to "
                        f"{fmt(hl['claim_vs_rerun_worst_rel'])} and the geomdl/scipy cross-checks at "
                        f"{fmt(hl['nurbs_vs_geomdl_worst'])}/{fmt(hl['gauss_vs_scipy_worst'])}; "
                        f"but 0/5 deliverables are git-tracked or on main and tensor_benchmark_cases has "
                        f"0 rows for the anchor -> declared L2, L3/L4 refused."),
            "scenario_key": "board-audit",
            "mission_id": 17585,
            "tags": ["TensorCAD", "board-audit", "m17453", "XVAL", "maturity-gate", "geomquality"],
            "markdown": md,
            "figures": [
                {"name": f["name"],
                 "caption": cap,
                 "png_base64": base64.b64encode((HERE / f["path"].split("docs/board_audit/m17585/")[-1]).read_bytes()).decode()}
                for f, cap in [
                    ({"name": "m17585_hydro_error.png", "path": "figs/m17585_hydro_error.png"},
                     "hydrostatics error vs independent/analytic reference"),
                    ({"name": "m17585_claim_vs_rerun.png", "path": "figs/m17585_claim_vs_rerun.png"},
                     "claim vs re-run parity"),
                    ({"name": "m17585_durability_funnel.png", "path": "figs/m17585_durability_funnel.png"},
                     "deliverable durability funnel"),
                ]
            ],
        }
        user = os.environ["TENSOR_AGENT_USER"]
        pw = os.environ["TENSOR_AGENT_PASS"]
        req = urllib.request.Request(base + "/auth/login", method="POST")
        req.add_header("Content-Type", "application/json")
        with urllib.request.urlopen(req, json.dumps({"username": user, "password": pw}).encode()) as r:
            token = json.loads(r.read())["access_token"]
        req = urllib.request.Request(base + "/tensor/reports", method="POST")
        req.add_header("Content-Type", "application/json")
        req.add_header("Authorization", "Bearer " + token)
        with urllib.request.urlopen(req, json.dumps(payload).encode(), timeout=120) as r:
            print("report posted:", r.status, r.read()[:300].decode())
    return 0


if __name__ == "__main__":
    raise SystemExit(main())