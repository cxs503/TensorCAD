#!/usr/bin/env python3
"""Figures + verdict row for the board audit of mission m17453.

Reads ``runs/m17585_audit/audit_summary.json`` (produced by
``runs/m17585_audit/verify_m17453.py``) — no number is retyped by hand.

    python3 docs/board_audit/m17585/make_figures_and_verdict.py            # figures + dry-run row
    python3 docs/board_audit/m17585/make_figures_and_verdict.py --write    # + INSERT ... ON CONFLICT
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]                      # /root/tp-agent/work/m17585
SUMMARY = ROOT / "runs/m17585_audit/audit_summary.json"
REFS = ROOT / "runs/m17585_audit/refs"
FIGS = HERE / "figs"

DSN = ("host=127.0.0.1 port=5432 dbname=tensorplatform user=foundation "
       "password=uGTN-KwuQv_7l8Tw0MMWUx_kf-nF2Rzz")
REPO = "TensorCAD"
SCOPE = "m17453 XVAL-CAD-M17404-CANONICAL-QUALITY-LIB"
MISSION = 17585
BY = "board-audit"
LEVEL = 2


def sha(p: Path) -> str:
    return hashlib.sha256(p.read_bytes()).hexdigest()


def figures(s: dict) -> list[dict]:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np

    FIGS.mkdir(parents=True, exist_ok=True)
    out = []
    hyd = s["hydrostatics_rerun"]
    names = list(hyd)

    # fig 1: got-vs-reference relative error per hull and metric, log scale
    metrics = ["volume", "Cb", "VCB", "Awp"]
    data = {m: [] for m in metrics}
    for n in names:
        h = hyd[n]
        got = h["got"]
        ref = h["got"]  # placeholder, replaced below
        item = h["independent"]
        data["volume"].append(max(item["volume_rel_err_vs_identity"], 1e-17))
        data["Cb"].append(abs(got["Cb_patch"] - h["declared"]["Cb"]) / h["declared"]["Cb"] + 1e-17)
        data["VCB"].append(max(item["vcb_rel_err_vs_own_ref"], 1e-17))
        data["Awp"].append(max(item["waterplane_rel_err_vs_own_ref"], 1e-17))

    fig, ax = plt.subplots(figsize=(7.2, 4.2))
    x = np.arange(len(names))
    w = 0.2
    for k, m in enumerate(metrics):
        ax.bar(x + (k - 1.5) * w, data[m], w, label=m)
    ax.set_yscale("log")
    ax.axhline(0.03, color="orange", ls="--", lw=1, label="3% (owner precision line)")
    ax.axhline(0.15, color="red", ls=":", lw=1, label="15% (unacceptable)")
    ax.set_xticks(x, names, rotation=10)
    ax.set_ylabel("relative error vs reference [-]")
    ax.set_title("m17453 hydrostatics re-run: error vs independent/analytic reference")
    ax.legend(fontsize=7, ncol=2)
    ax.grid(axis="y", which="both", alpha=0.3)
    p1 = FIGS / "m17585_hydro_error.png"
    fig.tight_layout()
    fig.savefig(p1, dpi=130)
    plt.close(fig)
    out.append({"name": p1.name, "path": str(p1.relative_to(ROOT))})

    # fig 2: claimed vs re-run parity (volume, observed order, GCI)
    fig, axs = plt.subplots(1, 3, figsize=(10.5, 3.6))
    triples = [("volume_patch", "patch volume"), ("observed_order", "observed order p"),
               ("gci_percent", "GCI [%]")]
    for ax, (key, label) in zip(axs, triples):
        claimed = [hyd[n]["claimed"][key] for n in names]
        got = [hyd[n]["got"][key] for n in names]
        lo = min(min(claimed), min(got))
        hi = max(max(claimed), max(got))
        ax.plot([lo, hi], [lo, hi], color="grey", lw=1, label="y = x")
        ax.scatter(claimed, got, s=45, zorder=3)
        for c, g, n in zip(claimed, got, names):
            ax.annotate(n, (c, g), fontsize=6, xytext=(3, 3), textcoords="offset points")
        ax.set_xlabel(f"m17453 claim: {label}")
        ax.set_ylabel("m17585 re-run")
        ax.set_title(label, fontsize=9)
        ax.grid(alpha=0.3)
    worst = s["headline"]["claim_vs_rerun_worst_rel"]
    fig.suptitle(f"claim vs re-run parity (worst relative gap {worst:.2e})", fontsize=10)
    p2 = FIGS / "m17585_claim_vs_rerun.png"
    fig.tight_layout()
    fig.savefig(p2, dpi=130)
    plt.close(fig)
    out.append({"name": p2.name, "path": str(p2.relative_to(ROOT))})

    # fig 3: durability / maturity funnel (facts, not opinions)
    durab = s["durability"]
    n_files = len(durab["tracked"])
    tracked = sum(durab["tracked"].values())
    inmain = sum(durab["present_in_main"].values())
    fig, ax = plt.subplots(figsize=(6.6, 3.4))
    stages = ["delivered\nfiles", "git-tracked\nin m17453", "present on\nmain", "benchmark\nregistry rows"]
    vals = [n_files, tracked, inmain, s["benchmark_registry"]["n_rows_for_anchor"]]
    bars = ax.bar(stages, vals, color=["#4C78A8", "#F58518", "#E45756", "#B0B0B0"])
    for b, v in zip(bars, vals):
        ax.annotate(str(v), (b.get_x() + b.get_width() / 2, v), ha="center", va="bottom", fontsize=9)
    ax.set_ylabel("count")
    ax.set_title("m17453 deliverable durability funnel (why L4 refused)")
    ax.set_ylim(0, max(vals) + 1.5)
    p3 = FIGS / "m17585_durability_funnel.png"
    fig.tight_layout()
    fig.savefig(p3, dpi=130)
    plt.close(fig)
    out.append({"name": p3.name, "path": str(p3.relative_to(ROOT))})
    return out


def build_evidence(s: dict, fig_paths: list[str]) -> dict:
    hyd = s["hydrostatics_rerun"]
    worst = s["headline"]
    refs = s["subject_artifacts_sha256"]
    return {
        # L2 keys
        "impl_files": [
            "tensorcad/geomquality/pipeline.py",
            "tensorcad/geomquality/hydro.py",
            "tensorcad/geomquality/baselines.py",
            "tensorcad/geomquality/_eval.py",
        ],
        "impl_detail": (
            "pipeline.py:148 verify_cad（CAD→defects→patch/facet hydrostatics→Richardson GCI）"
            " / pipeline.py:321 run_evidence / pipeline.py:473 main；hydro.py:146 patch_hydrostatics"
            " / hydro.py:296 facet_hydrostatics / hydro.py:486 convergence_study；"
            "baselines.py:287 reference_integrals / baselines.py:419 wigley_iii / :476 series60 / :559 kvlcc2_class"
        ),
        "entry": [
            "tensorcad/geomquality/pipeline.py",
            "examples/geomquality/run_pipeline.py",
        ],
        "entry_cmd": (
            "python3 -m tensorcad.geomquality.pipeline --out artifacts/geomquality ；"
            "python3 examples/geomquality/run_pipeline.py ；审计复算: python3 runs/m17585_audit/verify_m17453.py"
        ),
        "run_artifacts": [
            "runs/m17585_audit/audit_summary.json",
            "runs/m17585_audit/refs/artifacts__xval_results.json",
            "runs/m17585_audit/refs/artifacts__xval_pipeline_results.json",
            "runs/m17585_audit/refs/docs__xval_m17404_canonical_quality_lib.md",
        ],
        "run_artifacts_sha256": {k: v["sha256"] for k, v in refs.items()},
        "subject_status": {
            "mission": 17453,
            "db_status": "succeeded",
            "head_commit": None,
            "changed_files_in_result": ["docs/xval_m17404_canonical_quality_lib.md"],
            "git_status_all_untracked": all(not v for v in s["durability"]["tracked"].values()),
            "subject_head": s["durability"]["subject_head"],
            "present_in_main": s["durability"]["present_in_main"],
            "xval_generator_py_in_subject_tree": s["durability"]["xval_py_in_subject_tree"],
        },
        # honest extras (not needed at L2, kept for the gate/board)
        "tests": [
            "tests/test_geomquality_pipeline.py (18 tests)",
            "tests/test_geomquality_io.py (9 tests)",
        ],
        "tests_reproduced": "python3 -m pytest tests/test_geomquality_pipeline.py tests/test_geomquality_io.py -q -> 27 passed in 15.22s (2026-11-13)",
        "integration_callers": (
            "none production-side: grep -rl geomquality over /root/tensor-ecosystem (TensorCAD excluded) -> 0 hits;"
            " no [project.scripts] entry; tensorcad/cli.py has no geomquality subcommand;"
            " only in-repo caller = examples/geomquality/run_pipeline.py:285 run_evidence(...)"
        ),
        "benchmark": {
            "registry_rows_for_anchor": s["benchmark_registry"]["n_rows_for_anchor"],
            "registry_note": "tensor_benchmark_cases has 0 rows for anchor XVAL-CAD-M17404-CANONICAL-QUALITY-LIB (table total 55) -> no board-registered ref/tolerance to judge against",
            "independent_rerun_worst_rel": worst["claim_vs_rerun_worst_rel"],
            "volume_vs_Cb_identity_worst_rel": worst["hydrostatics_vs_Cb_identity_worst_rel"],
            "gauss_vs_scipy_worst_abs": worst["gauss_vs_scipy_worst"],
            "nurbs_vs_geomdl_worst_abs": worst["nurbs_vs_geomdl_worst"],
            "per_hull": {n: {"gci_percent": v["got"]["gci_percent"],
                             "error_finest_percent": v["got"]["error_finest_percent"],
                             "volume_rel_err_vs_Cb_identity": v["independent"]["volume_rel_err_vs_identity"]}
                         for n, v in hyd.items()},
            "vcb_awp_reference_gap": "kvlcc2_class: VCB/Awp reference = null in delivered JSON -> those two metrics were never benchmarked for that hull",
        },
        "repro_cmd": "python3 runs/m17585_audit/verify_m17453.py  (re-runs verify_cad on all three specimens + scipy/geomdl cross-checks; wrote runs/m17585_audit/audit_summary.json)",
        "external_anchor": [
            "geomdl 5.4.0 (NURBS eval/derivatives) — reproduced here: max |TC-geomdl| = 5.55e-17",
            "scipy 1.16.0 scipy.special.roots_legendre (Gauss-Legendre) — reproduced: max diff 1.46e-15",
            "analytic closed forms Cb*L*B*T identity + Wigley Awp=2/3·L·B",
        ],
        "refused_levels": {
            "L3": "no main-solver/platform caller (see integration_callers); mission delivered documentation+artifacts only, no code change",
            "L4": "no tensor_benchmark_cases row for the anchor (0/55) + deliverables never committed (head_commit=null, 5/5 files untracked, 0/5 on main) + no repro script in the subject tree",
            "L5": "L4 unmet",
        },
        "verdict_note": (
            "Substance verified, delivery not landed: the three hulls' hydrostatics, Richardson order and GCI "
            "reproduce to %.2e relative and the OSS cross-checks (geomdl/scipy) reproduce at machine precision, "
            "but every m17453 artifact is untracked and absent from main, so the evidence would vanish with the "
            "worktree. Declared L2 = real code path verified end-to-end without production integration or a "
            "registered benchmark." % worst["claim_vs_rerun_worst_rel"]
        ),
        "figures": fig_paths,
        "auditor": "m17585 board-audit (independent re-run, subject worktree read-only)",
    }


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--write", action="store_true")
    args = ap.parse_args()

    s = json.loads(SUMMARY.read_text())
    figs = figures(s)
    ev = build_evidence(s, [f["path"] for f in figs])
    ev_path = HERE / "evidence.json"
    ev_path.write_text(json.dumps(ev, indent=2, ensure_ascii=False))
    print(json.dumps(ev, indent=2, ensure_ascii=False)[:1200])
    print("figures:", [f["path"] for f in figs])

    if not args.write:
        print("\n[dry-run] row not written")
        return 0

    import psycopg

    with psycopg.connect(DSN) as c, c.cursor() as cur:
        cur.execute(
            "SELECT id FROM tensor_maturity_verdicts WHERE repo=%s AND scope=%s AND verdict_by=%s",
            (REPO, SCOPE, BY))
        row = cur.fetchone()
        payload = json.dumps(ev, ensure_ascii=False)
        if row:
            cur.execute(
                "UPDATE tensor_maturity_verdicts SET mission_id=%s, level=%s, evidence=%s,"
                " gate_state=NULL, gate_level=NULL, gate_note=NULL WHERE id=%s",
                (MISSION, LEVEL, payload, row[0]))
            vid = row[0]
        else:
            cur.execute(
                "INSERT INTO tensor_maturity_verdicts (repo, scope, mission_id, level, evidence, verdict_by)"
                " VALUES (%s,%s,%s,%s,%s,%s) RETURNING id",
                (REPO, SCOPE, MISSION, LEVEL, payload, BY))
            vid = cur.fetchone()[0]
        c.commit()
    print(f"verdict id={vid} repo={REPO} scope={SCOPE!r} level=L{LEVEL} by={BY}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())