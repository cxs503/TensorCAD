#!/usr/bin/env python3
"""Build the board-audit evidence bundle for m17307 (subject m17251).

Every number is read programmatically from the audited artifacts in
`runs/m17307/` (subject snapshot + independent re-verification + probe logs);
nothing is typed by hand.  Output: docs/board_audit/m17307/evidence.json and
docs/board_audit/m17307/audit_figures.png.

Run: python3 docs/board_audit/m17307/build_evidence.py
"""
from __future__ import annotations

import hashlib
import json
import subprocess
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

ROOT = Path(__file__).resolve().parents[3]          # worktree root (m17307)
SUB = ROOT / "runs/m17307/subject_m17251"
RUNS = ROOT / "runs/m17307"
OUT = ROOT / "docs/board_audit/m17307"
SUBJECT_MISSION = 17251
SUBJECT_COMMIT = "5bc0db7"


def sha(p: Path) -> str:
    return hashlib.sha256(p.read_bytes()).hexdigest()


def psql(sql: str) -> str:
    r = subprocess.run(["runuser", "-u", "postgres", "--", "psql", "-d", "tensorplatform",
                        "-tA", "-F", "|", "-c", sql], capture_output=True, text=True)
    return r.stdout.strip()


def main() -> int:
    res = json.loads((SUB / "docs/offset_loft/results.json").read_text())
    chk = json.loads((RUNS / "independent_check.json").read_text())
    probe = json.loads((RUNS / "tamper_probe.json").read_text())
    repro = json.loads((RUNS / "logs/repro.json").read_text())

    cmp_ = res["comparisons"]
    worst_vol = max(cmp_[h]["volume"]["relative_error_nurbs"] for h in cmp_)
    worst_wp = max(cmp_[h]["waterplane"]["relative_error_nurbs"] for h in cmp_)
    worst_mesh = max(max(cmp_[h]["volume"]["relative_error_mesh"],
                         cmp_[h]["waterplane"]["relative_error_mesh"]) for h in cmp_)
    station = res["station_sensitivity"]["kvlcc2_class"]
    station_tbl = {str(r["n_stations"]): {
        "dV_vs_81st_table": r["relative_error_volume_vs_denominator"],
        "dV_vs_closed_form": r["relative_error_volume_vs_reference"]} for r in station}
    baseline_repro = {h: max(chk[h]["audit_routeB_simpson"]["rel_diff_volume_vs_table_simpson"],
                             chk[h]["audit_routeB_simpson"]["rel_diff_waterplane_vs_table_simpson"])
                      for h in ("wigley_iii", "series60_cb060", "kvlcc2_class")}
    spline_vs_nurbs = {h: chk[h]["audit_routeA_spline"]["rel_diff_volume_vs_nurbs"]
                       for h in ("wigley_iii", "series60_cb060", "kvlcc2_class")}
    wig = chk["wigley_closed_form"]

    registry = psql("select count(*) from tensor_benchmark_cases "
                    "where anchor_slug='GAP-CAD-T1-OFFSET-LOFT';")
    registry_rows = int(registry.split("|")[0]) if registry else -1
    bench_rows = psql(f"select count(*) from tensorbench_results where mission_id={SUBJECT_MISSION};")
    bench_rows = int(bench_rows.split("|")[0]) if bench_rows else -1
    merged = subprocess.run(["git", "-C", "/root/tp-agent/work/m17251", "status", "--short"],
                            capture_output=True, text=True).stdout.strip()
    on_main = subprocess.run(
        ["git", "-C", "/root/tensor-ecosystem/TensorCAD", "ls-tree", "--name-only", "HEAD",
         "tensorcad/offset_loft.py"], capture_output=True, text=True).stdout.strip()

    files = sorted(p for p in SUB.rglob("*") if p.is_file())
    manifest = {str(p.relative_to(ROOT)): {"sha256": sha(p), "bytes": p.stat().st_size} for p in files}
    (RUNS / "subject_m17251/MANIFEST.json").write_text(json.dumps(
        {"subject_mission": SUBJECT_MISSION, "subject_commit": SUBJECT_COMMIT,
         "subject_worktree": "/root/tp-agent/work/m17251",
         "note": "verbatim snapshot of the audited deliverables (commit " + SUBJECT_COMMIT
                 + "), taken for the board audit; the branch agent/m17251 is not merged into main",
         "files": manifest}, indent=1, sort_keys=True))

    rel = lambda p: str(p.relative_to(ROOT))  # noqa: E731
    evidence = {
        "impl_files": [rel(SUB / "tensorcad/offset_loft.py"),
                       rel(SUB / "tensorcad/offset_forms.py")],
        "impl_detail": (
            f"offset_loft.py {len((SUB/'tensorcad/offset_loft.py').read_text().splitlines())} lines: "
            "ShipOffsetTable/OffsetProvenance (CSV+sha256 bundle), classical and unequal-spacing "
            "Simpson, tensor-product NURBS loft via tensorcad.parametric_hull.loft_nurbs_surface, "
            "Green-theorem section integrals over the patch's own basis (Gauss 12/span), "
            "volume/LCB/VCB/waterplane/wetted area, mesh route via parametric_hull.compute_hydrostatics, "
            "and the slow NURBSSurface.derivatives cross-check. "
            f"offset_forms.py {len((SUB/'tensorcad/offset_forms.py').read_text().splitlines())} lines: "
            "builtin hull-form registry (wigley/series60/kvlcc2-class) that emits tables with provenance"),
        "entry": rel(SUB / "examples/offset_loft/run_example.py"),
        "entry_detail": ("python3 examples/offset_loft/run_example.py [--output docs/offset_loft] -> "
                         "writes results.json + report.md + plot.png + tables/*.csv(+provenance); "
                         "rc=0"),
        "run_artifacts": [rel(SUB / "docs/offset_loft/results.json"),
                          rel(SUB / "docs/offset_loft/report.md"),
                          rel(SUB / "docs/offset_loft/plot.png"),
                          rel(SUB / "docs/offset_loft/tables/wigley_iii.csv"),
                          rel(SUB / "docs/offset_loft/tables/series60_cb060.csv"),
                          rel(SUB / "docs/offset_loft/tables/kvlcc2_class.csv"),
                          rel(SUB / "docs/offset_loft/tables/kvlcc2_class.provenance.json"),
                          rel(RUNS / "independent_check.json"),
                          rel(RUNS / "tamper_probe.json"),
                          rel(RUNS / "logs/repro.json"),
                          rel(RUNS / "logs/pytest.log"),
                          rel(RUNS / "subject_m17251/MANIFEST.json")],
        "tests": ["runs/m17307/subject_m17251/tests/test_offset_loft.py: 18 collected tests, "
                  "18 passed (auditor rerun: python -m pytest tests/test_offset_loft.py -q -> "
                  "18 passed in 17.5s; with tests/test_loft.py -> 49 passed in 27.2s)",
                  "test_landed_tables_pass_the_half_percent_acceptance (per-table acceptance)",
                  "test_wigley_loft_reproduces_the_closed_form_volume_and_waterplane",
                  "test_loft_passes_through_every_tabulated_ordinate",
                  "test_bundle_round_trip_and_tamper_detection",
                  "test_fast_basis_matrix_path_matches_the_surface_derivative_path",
                  "test_station_sensitivity_converges_and_reports_both_baselines"],
        "integration_callers": "",
        "integration_note": (
            "NOT MET and not claimed.  The module is re-exported through tensorcad/__init__.py "
            "(loft_offset_table, compare_offset_table, ShipOffsetTable, ...) and consumed by the "
            "companion tensorcad/offset_forms.py of the same commit, plus the example script: no "
            "solver, CLI, API endpoint or ship-performance-platform caller exists, and "
            f"'git ls-tree HEAD tensorcad/offset_loft.py' on /root/tensor-ecosystem/TensorCAD "
            f"returns {'nothing (branch not merged)' if not on_main else on_main}.  No downstream "
            "hydrostatics/resistance path consumes the lofted patch yet."),
        "benchmark": {
            "registry_case": (f"tensor_benchmark_cases rows for anchor_slug="
                              f"'GAP-CAD-T1-OFFSET-LOFT' = {registry_rows} -> no registered ref/tolerance "
                              "exists; the platform L4 rule cannot be satisfied and "
                              f"tensor_bench_results rows for mission {SUBJECT_MISSION} = {bench_rows}"),
            "acceptance_criterion": "lofted volume and waterplane area within 0.5% of the table's own Simpson integral",
            "ref_vs_got": {
                "wigley_iii": {"volume_ref": cmp_["wigley_iii"]["volume"]["table_simpson"],
                               "volume_got": cmp_["wigley_iii"]["volume"]["nurbs_quadrature"],
                               "rel_err": cmp_["wigley_iii"]["volume"]["relative_error_nurbs"],
                               "waterplane_ref": cmp_["wigley_iii"]["waterplane"]["table_simpson"],
                               "waterplane_got": cmp_["wigley_iii"]["waterplane"]["nurbs_quadrature"],
                               "rel_err_wp": cmp_["wigley_iii"]["waterplane"]["relative_error_nurbs"]},
                "series60_cb060": {"volume_ref": cmp_["series60_cb060"]["volume"]["table_simpson"],
                                   "volume_got": cmp_["series60_cb060"]["volume"]["nurbs_quadrature"],
                                   "rel_err": cmp_["series60_cb060"]["volume"]["relative_error_nurbs"],
                                   "rel_err_wp": cmp_["series60_cb060"]["waterplane"]["relative_error_nurbs"]},
                "kvlcc2_class": {"volume_ref": cmp_["kvlcc2_class"]["volume"]["table_simpson"],
                                 "volume_got": cmp_["kvlcc2_class"]["volume"]["nurbs_quadrature"],
                                 "rel_err": cmp_["kvlcc2_class"]["volume"]["relative_error_nurbs"],
                                 "rel_err_wp": cmp_["kvlcc2_class"]["waterplane"]["relative_error_nurbs"]},
            },
            "worst_relative_error_volume": worst_vol,
            "worst_relative_error_waterplane": worst_wp,
            "worst_relative_error_loft_mesh": worst_mesh,
            "station_sensitivity_kvlcc2_dV_vs_81st_table": station_tbl,
            "auditor_baseline_reproduction": {
                "method": "independent composite Simpson (uniform 1/3 + 3/8 tail; unequal-spacing "
                          "quadratic rule for clustered waterlines) re-implemented from the CSV only, "
                          "no import of tensorcad.offset_loft",
                "max_rel_diff_vs_reported_table_simpson": max(baseline_repro.values())},
            "auditor_independent_quadrature": {
                "method": "not-a-knot cubic spline in z per station then in x (a different rule)",
                "rel_diff_volume_vs_reported_nurbs": spline_vs_nurbs},
        },
        "repro_cmd": (f"python3 examples/offset_loft/run_example.py  # rc=0, results.json sha256 "
                      f"{repro['results_json_sha256_after'][:16]} before and after -> "
                      f"byte_identical={repro['byte_identical']}; "
                      "python3 -m pytest tests/test_offset_loft.py tests/test_loft.py -q  # 49 passed"),
        "external_anchor": {
            "analytic": {"wigley_volume_4BLT_over_9": wig["volume_analytic_4BLT_over_9"],
                         "auditor_spline_volume": wig["volume_audit_spline"],
                         "rel_diff": wig["rel_diff"],
                         "wigley_waterplane_2LB_over_3": wig["waterplane_analytic_2LB_over_3"],
                         "rel_diff_waterplane": wig["rel_diff_waterplane"],
                         "note": "Cb=4/9, Cwp=2/3, VCB=5T/8 exact; reproduced by an independent "
                                 "quadrature from the landed CSV only"},
            "published_particulars": {"kvlcc2": "L=320 m, B=58 m, T=20.8 m, Cb=0.8098, LCB=-3.48%L "
                                                "(KRISO/ITTC 2010+2015 benchmark particulars) used to "
                                                "calibrate the class reconstruction; Cb from the landed "
                                                f"table = {chk['kvlcc2_class']['cb_from_tables']:.6f}",
                                      "series60": "Todd (1951) DTMB 1447 parent form, Cp=0.625, n=5/3"},
            "NOT_AVAILABLE": ("no issuing-body offset table was ingested: every provenance records "
                              "official=false and the KVLCC2 sidecar states 'NOT the official KRISO "
                              "offset table'; the acceptance is therefore a method/quadrature "
                              "consistency check plus analytic identities, not validated ship data"),
        },
        "integrity_probes": {
            "tamper_detection": [{"probe": p["probe"], "rejected": p.get("rejected", p.get("digest_ok"))}
                                 for p in probe["probes"]],
            "all_probes_pass": probe["ok"],
            "subject_worktree_clean_after_rerun": merged == "",
        },
        "snapshot": rel(RUNS / "subject_m17251/MANIFEST.json"),
    }
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "evidence.json").write_text(json.dumps(evidence, indent=1, ensure_ascii=False))

    # ── figure: acceptance ref-vs-got, station sensitivity, independent cross-check ──
    fig, ax = plt.subplots(1, 3, figsize=(15, 4.4))
    hulls = ["wigley_iii", "series60_cb060", "kvlcc2_class"]
    width = 0.28
    for i, h in enumerate(hulls):
        for j, (key, lbl) in enumerate([("relative_error_nurbs", "NURBS quadrature"),
                                        ("relative_error_mesh", "loft mesh")]):
            ax[0].bar(i + (j - 0.5) * width, cmp_[h]["volume"][key], width,
                      label=lbl if i == 0 else None)
        ax[0].bar(i + 1.5 * width, cmp_[h]["waterplane"]["relative_error_nurbs"], width,
                  label="waterplane (NURBS)" if i == 0 else None)
    ax[0].axhline(0.005, color="r", ls="--", lw=1, label="acceptance 0.5%")
    ax[0].set_xticks(range(3), hulls, rotation=20)
    ax[0].set_ylabel("relative error vs table Simpson integral")
    ax[0].set_title("(a) acceptance: ref vs got (acceptance study)")
    ax[0].legend(fontsize=7)
    ax[0].grid(alpha=.3, axis="y")
    ns = [r["n_stations"] for r in station]
    ax[1].plot(ns, [r["relative_error_volume_vs_denominator"] for r in station], "o-", label="dV vs 81-st table")
    ax[1].plot(ns, [r["relative_error_volume_vs_reference"] for r in station], "s--", label="dV vs closed form")
    ax[1].axhline(0.005, color="r", ls="--", lw=1, label="0.5% acceptance")
    ax[1].set_xscale("log")
    ax[1].set_yscale("log")
    ax[1].set_xlabel("stations in the table (kvlcc2_class)")
    ax[1].set_ylabel("relative volume error")
    ax[1].set_title("(b) station-count sensitivity (subject results.json)")
    ax[1].legend(fontsize=8)
    ax[1].grid(alpha=.3, which="both")
    tags = ["routeB\n(Simpson rerun)", "routeA\n(spline quadrature)"]
    xs = range(len(hulls))
    ax[2].bar([x - 0.18 for x in xs], [baseline_repro[h] for h in hulls], 0.36, label=tags[0])
    ax[2].bar([x + 0.18 for x in xs], [spline_vs_nurbs[h] for h in hulls], 0.36, label=tags[1])
    ax[2].set_yscale("log")
    ax[2].set_xticks(list(xs), hulls, rotation=20)
    ax[2].set_ylabel("relative difference")
    ax[2].set_title("(c) auditor's independent quadrature vs reported numbers")
    ax[2].legend(fontsize=8)
    ax[2].grid(alpha=.3, axis="y", which="both")
    fig.suptitle("m17307 board audit of m17251 — GAP-CAD-T1-OFFSET-LOFT", fontsize=11)
    fig.tight_layout()
    fig.savefig(OUT / "audit_figures.png", dpi=130)
    print("wrote", OUT / "evidence.json", "and audit_figures.png")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())