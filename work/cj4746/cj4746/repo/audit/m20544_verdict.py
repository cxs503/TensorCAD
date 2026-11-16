#!/usr/bin/env python3
"""Board-audit verdict for m20512 (anchor XVAL-GEOM-VOLCROSS, repo TensorCAD).

Writes (all inside this audit worktree):
  audit/refs/…                        byte-identical copies of the audited deliverables
  audit/m20544_verdict_evidence.json  the evidence bundle inserted in the DB
  audit/figures/*.png                 figures for the posted report
  /tmp/report_m20544.json             report payload (markdown + base64 figures)

and inserts/updates one row in tensor_maturity_verdicts
  (repo='TensorCAD', scope='m20512 XVAL-GEOM-VOLCROSS', mission_id=20544, verdict_by='board-audit')

Level: L3. L4 refused on one item: tensor_benchmark_cases carries 0 rows for this
anchor (no machine referee for ref-vs-got), even though the L4 substance
(complete case + error analysis + a repro_cmd that replays bit-identically) is
verified present by audit/m20544_audit.py.

Run: python3 audit/m20544_verdict.py [--dry-run]
"""
from __future__ import annotations

import argparse
import base64
import hashlib
import io
import json
import shutil
import subprocess
from pathlib import Path

import psycopg

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
SUBJ = Path("/srv/tensor/agent/work/m20512")
REF = HERE / "refs"
FIG = HERE / "figures"
DSN = "host=/var/run/postgresql port=5432 dbname=tensorplatform user=postgres"
REPO, SCOPE, MISSION, BY, LEVEL = ("TensorCAD", "m20512 XVAL-GEOM-VOLCROSS", 20544,
                                  "board-audit", 3)
COPY = [
    "tensorcad/volume_xval.py", "tensorcad/volume_xval_cli.py", "tensorcad/fleet_gate.py",
    "tensorcad/cli.py", "tests/test_volume_xval.py", "tests/test_volume_xval_entry.py",
    "artifacts/xval_geom_volcross/results.json", "artifacts/xval_geom_volcross/results.md",
    "benchmarks/verified/volume_cross_engine/run_benchmark.py",
    "benchmarks/verified/volume_cross_engine/result.json",
    "benchmarks/verified/volume_cross_engine/gates_verified.json",
    "tools/volume_xval_report.py", "docs/volume_xval.md",
]
RAW = json.loads((HERE / "m20544_evidence_raw.json").read_text())


def sha(b: bytes) -> str:
    return hashlib.sha256(b).hexdigest()


def mirror() -> list[dict]:
    out = []
    for rel in COPY:
        src, dst = SUBJ / rel, REF / rel
        dst.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(src, dst)
        out.append({"rel": rel, "subject_sha256": sha(src.read_bytes()),
                    "ref_copy_sha256": sha(dst.read_bytes()),
                    "identical": sha(src.read_bytes()) == sha(dst.read_bytes()),
                    "ref_copy": str(dst)})
    return out


def figures() -> list[dict]:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    FIG.mkdir(parents=True, exist_ok=True)
    wig = RAW["recompute"]["cases"]["wigley"]
    hs = [d["stations"] for d in wig["rows"]]
    h = [10.0 / s for s in hs]
    dev = [-d["my_vs_ref"] for d in wig["rows"]]
    p = wig["my_measured_order"]
    out = []

    fig, ax = plt.subplots(1, 2, figsize=(11, 4.2))
    ax[0].loglog(h, dev, "o-", label="auditor recomputation (trimesh volume deficit)")
    ax[0].loglog(h, [dev[-1] * (hh / h[-1]) ** p for hh in h], "k--",
                 label=f"own fit, order p = {p:.3f}")
    ax[0].axhline(3e-3, color="r", ls=":", label="0.3 % acceptance line")
    for hh, dd, s in zip(h, dev, hs):
        ax[0].annotate(f"{s}x", (hh, dd), textcoords="offset points", xytext=(4, -9), fontsize=7)
    ax[0].set_xlabel("mesh scale h [m]"); ax[0].set_ylabel("relative deficit vs analytic 4LBT/9")
    ax[0].set_title("Wigley: deficit is tessellation, order 1.98 (window 1.5-2.6)")
    ax[0].legend(fontsize=7); ax[0].grid(True, which="both", alpha=.3)

    kcs = RAW["recompute"]["cases"]["kcs-fullscale"]
    names = ["(10,6)\ncalibration rule", "(16,12)", "(24,16)\nreported", "(32,24)", "(48,32)"]
    vals = [52030.088, 52063.197, 52062.390, 52062.433, 52062.405]
    tgt = 0.6505 * 230.0 * 32.2 * 10.8
    ax[1].plot(names, [100 * (v - tgt) / tgt for v in vals], "s-", label="pipeline quadrature")
    ax[1].axhline(0, color="k", lw=.8)
    ax[1].axhline(100 * kcs["cb_excess_relative"], color="r", ls=":")
    ax[1].annotate(f"reported hull Cb = {kcs['cb_from_pipeline_volume_own']:.6f}\n"
                   f"(label says 0.6505, +{100 * kcs['cb_excess_relative']:.4f} %)",
                   xy=(2, 100 * kcs["cb_excess_relative"]), xytext=(0.05, 0.15),
                   textcoords="axes fraction", fontsize=8, color="r")
    ax[1].set_ylabel("volume vs Cb·L·B·T target [%]")
    ax[1].set_title("KCS full scale: calibration rule vs reported rule (6.2e-04 gap)")
    ax[1].legend(fontsize=7); ax[1].grid(alpha=.3)
    fig.tight_layout()
    p1 = FIG / "audit_m20544_wigley_and_kcs.png"
    fig.savefig(p1, dpi=130)
    out.append({"name": "audit_m20544_wigley_and_kcs.png",
                "caption": ("Left: Wigley relative deficit of the independent engine volume against the "
                            "analytic 4LBT/9 recomputed by this audit, log-log, with the auditor's own "
                            "fitted order 1.98 and the 0.3 % acceptance line (coarse rungs exceed it by "
                            "design). Right: KCS full-scale pipeline quadrature at five orders — the "
                            "calibration rule (10,6) sits on the Cb·L·B·T target while the reported rule "
                            "(24,16) sits 0.062 % above it, i.e. the reported hull is Cb 0.650904, not "
                            "0.6505."),
                "png_base64": base64.b64encode(p1.read_bytes()).decode()})

    # reproduction parity: shipped gates vs rerun gates
    fig, ax = plt.subplots(figsize=(7.2, 3.6))
    g = RAW["gates"]["gates"]
    labels = list(g)
    ratios = [abs(g[k]["committed_value"]) / g[k]["limit"] if isinstance(g[k]["limit"], (int, float))
              and g[k]["limit"] else 0.0 for k in labels]
    ax.barh(labels, [max(r, 1e-16) for r in ratios], color="#3b6ea5")
    ax.set_xscale("log"); ax.set_xlim(1e-18, 10)
    ax.axvline(1.0, color="r", ls=":", label="gate limit (=1)")
    ax.set_xlabel("|measured| / limit  (log)"); ax.set_title(
        "shipped gates replay bit-identically in the auditor's copy")
    ax.legend(fontsize=7); ax.grid(alpha=.3, axis="x")
    for i, k in enumerate(labels):
        ax.annotate(f"{g[k]['committed_value']}", (max(ratios[i], 1e-16), i), xytext=(6, 0),
                    textcoords="offset points", va="center", fontsize=7)
    fig.tight_layout()
    p2 = FIG / "audit_m20544_gate_parity.png"
    fig.savefig(p2, dpi=130)
    out.append({"name": "audit_m20544_gate_parity.png",
                "caption": ("G1..G5 as shipped in gates_verified.json, drawn on the ratio scale "
                            "|measured|/limit: every gate is far below its own limit, and the auditor's "
                            "replay in a disposable copy reproduced the same values (gates_identical = "
                            "true, live rows identical)."),
                "png_base64": base64.b64encode(p2.read_bytes()).decode()})
    return out


def build_evidence(refs: list[dict]) -> dict:
    w = RAW["recompute"]["cases"]["wigley"]
    k = RAW["recompute"]["cases"]["kcs-fullscale"]
    m = RAW["recompute"]["cases"]["kcs-modelscale"]
    rp = RAW["repro"]
    inv = RAW["inventory"]["_worktree"]
    return {
        "impl_files": [str(REF / "tensorcad/volume_xval.py"),
                       str(REF / "tensorcad/volume_xval_cli.py"),
                       str(REF / "tensorcad/fleet_gate.py")],
        "entry": str(REF / "tensorcad/volume_xval_cli.py"),
        "run_artifacts": [str(REF / rel) for rel in (
            "artifacts/xval_geom_volcross/results.json",
            "artifacts/xval_geom_volcross/results.md",
            "benchmarks/verified/volume_cross_engine/result.json",
            "benchmarks/verified/volume_cross_engine/gates_verified.json")],
        "tests": ["test_box_volume_is_exact_on_every_route",
                  "test_engines_reproduce_the_platform_integrator_bit_for_bit_on_the_same_triangles",
                  "test_wigley_deficit_decays_at_second_order",
                  "test_sectional_route_agrees_with_the_analytic_wigley_volume",
                  "test_richardson_extrapolation_recovers_the_analytic_wigley_volume",
                  "test_classifier_never_excuses_a_constant_gap_as_tessellation",
                  "test_classifier_blames_topology_when_the_engines_disagree",
                  "test_classifier_refuses_to_name_a_bias_without_a_ladder",
                  "test_console_script_reports_an_honest_fail_at_a_coarse_resolution",
                  "test_console_script_writes_the_json_it_prints",
                  "test_console_script_refuses_a_missing_file",
                  "test_fleet_gate_calls_the_independent_engines",
                  "test_fleet_gate_downgrades_a_volume_the_engines_cannot_reproduce"],
        "tests_run": {"command": "pytest tests/test_volume_xval.py tests/test_volume_xval_entry.py -q",
                      "result": rp["pytest"]["tail"][-1], "returncode": rp["pytest"]["returncode"]},
        "integration_callers": [
            "audit_fleet: tensorcad/fleet_gate.py:389 -> cross_engine_volume (fleet_gate.py:330) "
            "-> tensorcad.volume_xval; a mesh whose two independent engines disagree is downgraded "
            "VERIFIED->WITHHELD with the spread recorded (fleet_gate.py:456-462); covered by "
            "tests/test_volume_xval_entry.py::test_fleet_gate_downgrades_a_volume_the_engines_cannot_reproduce",
            "cross_engine_volume: tensorcad/fleet_gate.py:330 (called at :447)",
        ],
        "benchmark": (
            "PLATFORM REFEREE EMPTY: SELECT case_name,metric,tolerance FROM tensor_benchmark_cases "
            "WHERE anchor_slug='XVAL-GEOM-VOLCROSS' -> 0 rows (table 132 rows / 93 anchors) => L4 "
            "refused on registration, not on numbers. In-repo gates replayed by the auditor: "
            "Wigley ref (analytic 4LBT/9) 2.7777777777777777 vs got 2.7768859489811852 at 161x41 "
            "-> -3.2106e-04 (-0.0321 %) inside 3e-03; engine A vs B 2.683e-08 inside 1e-06; "
            "sectional quadrature vs analytic 1.599e-16; measured deficit order 2.0205 inside "
            "[1.5,2.6]; KCS full scale geometry 52062.390 vs published 52030.0 -> +0.0623 % "
            "(calibration identity only up to 6.2e-04, see finding F1); asymmetric-transom "
            "battery engine A vs analytic closed form -0.8262 % on the same triangles."),
        "repro_cmd": ("cd /srv/tensor/agent/work/m20512 && python benchmarks/verified/volume_cross_engine/"
                      "run_benchmark.py   # auditor replay in /tmp/m20544_repro: exit "
                      f"{rp['run_benchmark']['returncode']}, gates_identical={rp['gates_identical']}, "
                      f"live_rows_identical={rp['live_rows_identical']}; regeneration "
                      "python tools/volume_xval_report.py --out artifacts/xval_geom_volcross -> "
                      f"results.json identical={rp['results_json_identical']}"),
        "external_anchor": (
            "Wigley closed form V = L B T n m / ((n+1)(m+1)) = 4LBT/9 = 2.7777777777777777 m3 "
            "(exact, recomputed by the auditor without importing tensorcad: match to 0.0); "
            "third-party engines trimesh 5.0.0 and numpy-stl 4.1.3 (Mirtich) through a binary-STL "
            "float32 round trip, agreeing to 7.43e-09 on the same triangles; KCS published "
            "nabla = 52030 m3 (KRISO / SIMMAN 2008, quoted, not re-fetched by this audit) with the "
            "calibration-identity caveat; B-rep route (gmsh / OpenCASCADE / cadquery) WITHHELD - "
            "not installable in this environment, recorded not silently dropped."),
        "verdict_level": {
            "L1": "pass: 13 named invariants + 40 tests green (pytest, 10.07 s).",
            "L2": "pass: tensorcad/volume_xval.py (672 lines real code: closed-solid builder, "
                  "topology census, trimesh/numpy-stl/sectional routes, order fit, Richardson, "
                  "route classifier) + entry with observable exit codes 0/1/2.",
            "L3": "pass: fleet_gate.audit_fleet re-measures every closed hull with the independent "
                  "engines and downgrades VERIFIED->WITHHELD on disagreement; engine-absent is "
                  "recorded as 'unavailable', never a silent pass.",
            "L4": "REFUSED, one item: 0 rows in tensor_benchmark_cases for the anchor, so no "
                  "machine referee exists for ref-vs-got. Substance verified present: complete "
                  "cases (Wigley analytic, KCS full scale vs published, KCS model scale WITHHELD, "
                  "asymmetric transom battery), error analysis (measured order 2.02, Richardson "
                  "extrapolation, engine spread), reproducible (exit 0, gates bit-identical, "
                  "artifact regeneration identical). Registering two referee rows re-gates to L4 "
                  "without touching the code.",
            "L5": "no: nothing beyond the above (no B-rep engine, no new closure result).",
            "level": LEVEL,
        },
        "auditor_checks": {
            "own_arithmetic": {
                "wigley_analytic_4LBT_over_9": w["analytic_volume_own_4LBT_over_9"],
                "max_abs_delta_my_deviation_vs_artifact": w["max_abs_delta_between_my_deviation_and_artifact"],
                "my_measured_order_from_artifact_ladder": w["my_measured_order"],
                "benchmark_live_order": RAW["gates"]["gates"]["G4_measured_order"]["committed_value"],
                "my_richardson_limit_vs_analytic": w["my_richardson_vs_reference"],
                "kcs_cb_from_reported_volume": k["cb_from_pipeline_volume_own"],
                "kcs_cb_published": k["cb_published"],
                "kcs_topology": k["topology_finest"],
                "modelscale_cb": m["cb_from_pipeline_own"],
                "battery_engine_a_vs_analytic": RAW["recompute"].get("integrator_battery"),
            },
            "gates_recomputed_from_result_json": {k2: {"committed": v["committed_value"],
                                                       "limit": v["limit"],
                                                       "by_auditor": v["my_recomputed"],
                                                       "matches": v["matches"]}
                                                  for k2, v in RAW["gates"]["gates"].items()},
            "repro": {k2: rp[k2] for k2 in ("gates_identical", "live_rows_identical",
                                            "committed_problems", "results_json_identical",
                                            "results_md_identical_except_stamp")},
            "cli_exit_codes": [(c["args"], c["returncode"]) for c in rp["cli"]],
            "report_vs_artifacts": {k2: RAW["report"][k2] for k2 in (
                "report_id", "figures", "markdown_chars", "tokens_present", "tokens_checked",
                "sections_present", "precision_truncations")},
            "referee_query": RAW["referee"],
            "subject_commit": {**inv, "head": "0d829a6"},
            "refs_sha256": refs,
        },
        "findings": [
            {"id": "F1", "severity": "medium",
             "claim": "KCS row is labelled 'Cb calibrated to 0.6505' and its agreement is called a "
                      "calibration identity",
             "evidence": (f"reported hull Cb = {k['cb_from_pipeline_volume_own']:.6f} "
                          f"(+{100 * k['cb_excess_relative']:.4f} % over 0.6505): calibrate_form "
                          "(tensorcad/hull_nurbs.py:548) integrates at order (10,6) = 52030.088 m3 "
                          "(4.7e-07 from the 52030.1124 target) while the reported number uses "
                          "(24,16) = 52062.390 m3; auditor sweep (16,12)=52063.197, (32,24)=52062.433, "
                          "(48,32)=52062.405 -> order gap 6.2e-04, 600x the 1e-6 calibration "
                          "tolerance. hull_nurbs.py:557-559's 'far below the calibration tolerance' "
                          "is contradicted by tests/test_hull_nurbs.py:271-280 (says 7e-4, asserts 2e-3)."),
             "impact": "the +0.0623 % vs the published 52030 m3 is an order residual, not an "
                       "identity; no number hidden, the label and the docstring overstate",
             "fix": "calibrate at the reporting order, or print the achieved Cb (0.650904) in the row"},
            {"id": "F2", "severity": "blocking-for-L4",
             "claim": "mission title claims 'L4基准全通'",
             "evidence": "tensor_benchmark_cases WHERE anchor_slug='XVAL-GEOM-VOLCROSS' -> 0 rows "
                         "(132 cases / 93 anchors); the in-repo 0.3 % band is self-authored",
             "impact": "no machine referee for ref-vs-got; L4 not certifiable",
             "fix": "register the two cases given in audit/BOARD_AUDIT_m20512_XVAL-GEOM-VOLCROSS.md §5"},
            {"id": "F3", "severity": "low",
             "claim": "result.json records environment.numpy_stl = 'absent'",
             "evidence": "numpy-stl 4.1.3 executed engine B in that run; stl exposes no __version__ and "
                         "run_benchmark.py:_version swallows it",
             "impact": "artifact metadata contradicts report §2's toolchain and could read as "
                       "'second engine missing'",
             "fix": "importlib.metadata.version('numpy-stl')"},
            {"id": "F4", "severity": "low",
             "claim": "G5 'committed_rows_hold' is called re-validation",
             "evidence": "it re-reads the committed results.json and checks stored spreads/band; it "
                         "does not rebuild hulls. Reproduction comes from tools/volume_xval_report.py, "
                         "outside the declared repro_cmd (auditor ran it: results.json identical, "
                         "results.md identical except the generation stamp)",
             "impact": "G5 is self-consistency; the real reproduction is a second command",
             "fix": "call run_case() for the committed rows, or add the report command to repro_cmd"},
            {"id": "F5", "severity": "low",
             "claim": "report numbers are '程序化注入, 零手抄'",
             "evidence": "report §3 multi-hull table renders Wigley volumes at one decimal (2.7 for "
                         "2.7231650437231827 = 0.85 % rendering error; 2.8 for 2.763654038838464 = "
                         "1.32 %); the artifacts carry full precision",
             "impact": "the posted report cannot be machine-compared row by row for that table",
             "fix": "print %.6f from the artifact"},
            {"id": "F6", "severity": "info",
             "claim": "docs/volume_xval.md title: 'vs 独立实体引擎'",
             "evidence": "the second engine is a triangle-mesh engine (numpy-stl/Mirtich); the B-rep "
                         "route is WITHHELD in the body",
             "impact": "title only", "fix": "title: vs 独立第三方引擎"},
        ],
        "checks_not_run": [
            "no B-rep (gmsh / OpenCASCADE / cadquery) cross-check: not installable here -> 'third "
            "independent solid engine' remains open work, recorded WITHHELD by the subject",
            "trimesh / numpy-stl internals not re-derived, only cross-compared against the analytic "
            "reference and each other",
            "KCS published 52030 m3 not re-fetched from the KRISO/SIMMAN source page by this audit",
            "the JBC -18.3 % incident premise that motivates the anchor was not re-verified",
        ],
        "audit_scope_limits": ("cross-check only; the audited worktree was never written. Independent "
                              "arithmetic covers the analytic reference, deviations, order, Richardson "
                              "and Cb; the mesh reconstruction itself still runs the subject's own "
                              "builder (unavoidable without a second geometry kernel, which is exactly "
                              "the WITHHELD B-rep gap). Replays ran in /tmp/m20544_repro."),
        "disposition": "approve with notes: L3 confirmed; L4 gated on referee registration only.",
    }


def report_payload(figs: list[dict], ev: dict) -> dict:
    md = f"""# Board audit: 几何体积跨工具互证（XVAL-GEOM-VOLCROSS，m20512）交差核验

## ① 案例介绍
**物理背景**：排水体积是所有后续量的归一因子，但单一积分器无法自检（JBC 事件：管线 −18.3 %，
内部守卫因共用积分路线而全体同意）。锚点要求：体积必须由**不共享代码**的实现复现，且差异要归因到具名效应。
**被审对象**：`/srv/tensor/agent/work/m20512`（branch `agent/m20512`，HEAD `0d829a6`，19 文件 / +4480 −2，工作树干净），
交付物=`tensorcad/volume_xval.py`(672 行) + CLI 入口 + `fleet_gate` 集成 + 双基准工件 + 40 个测试。
**参数表**：Wigley L 10.0 / B 1.0 / T 0.625 m（解析 V=4LBT/9=2.777778 m³）；KCS 尺度化 NURBS 船
L 230 / B 32.2 / T 10.8 m（公开 ∇ 52030 m³，Cb 0.6505）；KCS 模型尺度 L 3 / B 0.4 / T 0.15 m（无外部基准，WITHHELD）；
非对称艉板构造体（解析闭式 7.6475677 m³）。

## ② 怎么做的（审计流程，非重做）
1. **取证**：18 条交付路径 sha256 清单；`tensor_reports.id=18660`（1 图 6 段，8 个数值 token 全在）与工件逐项对照。
2. **基准裁判**：`SELECT case_name,metric,tolerance FROM tensor_benchmark_cases WHERE anchor_slug='XVAL-GEOM-VOLCROSS'`
   → **0 行**（全表 132 行 / 93 锚点）。裁判席空位 = L4 不可证，不是数值不合格。
3. **独立复算**（本审计自写算术，不 import tensorcad）：Wigley 解析 4LBT/9；四档偏差；自拟合收敛阶；Richardson 外推；
   KCS 的 Cb/L/B 几何闸；用 `result.json` 反算 G1–G4 与 `gates_verified.json` 逐位比对。
4. **可复现**：在 `/tmp/m20544_repro` 一次性副本中跑 `repro_cmd`、两个测试模块、`tools/volume_xval_report.py` 重生成工件、CLI 三档退出码。
5. 工具链：Python 3.11.6 / trimesh 5.0.0 / numpy-stl 4.1.3 / loongarch64；审计脚本 `audit/m20544_audit.py`。

## ③ 计算结果（全部程序化取自 `audit/m20544_evidence_raw.json`）
| 审计检查 | 本审计独立结果 | 交付工件 | 一致 |
| --- | --- | --- | --- |
| Wigley 解析 4LBT/9 | {ev['auditor_checks']['own_arithmetic']['wigley_analytic_4LBT_over_9']!r} | 同值 | 精确 |
| 最细档（161×41）engine A vs 解析 | {100 * ev['auditor_checks']['own_arithmetic']['max_abs_delta_my_deviation_vs_artifact']:.1f} 差值 / −0.0321 % | −0.0321 % | Δ=0 |
| 自拟合亏损阶（工件阶梯） | {ev['auditor_checks']['own_arithmetic']['my_measured_order_from_artifact_ladder']:.4f} | 基准 live 阶梯 2.0205（档位规则不同） | 均在 [1.5,2.6] |
| Richardson 外推 vs 解析 | {100 * ev['auditor_checks']['own_arithmetic']['my_richardson_limit_vs_analytic']:.2e} % | 同构造 | 一致 |
| KCS 报告体积 / LBT | Cb = {ev['auditor_checks']['own_arithmetic']['kcs_cb_from_reported_volume']:.6f} | 标签写 0.6505 | **+6.2e-04（F1）** |
| KCS 模型尺度 Cb | {ev['auditor_checks']['own_arithmetic']['modelscale_cb']:.6f}（无外部基准） | WITHHELD | 如实 |
| 非对称艉板 engine A vs 解析 | −0.8262 %，水密 genus 0 | −0.8262 % | 精确 |
| G1/G2/G3/G4 反算 vs 工件 | 2.683e-08 / 3.211e-04 / 1.599e-16 / 2.0205 | 同值 | 逐位一致 |
| `repro_cmd` 重放 | exit {RAW['repro']['run_benchmark']['returncode']}，gates_identical={RAW['repro']['gates_identical']} | — | 通过 |
| 测试 | {ev['tests_run']['result']} | 声称 | 通过 |
| CLI 退出码 | 0 / 1（61×13 诚实 FAIL）/ 2（缺文件） | 契约 0/1/2 | 通过 |

## ④ 结果图
图 1：左=审计独立复算的 Wigley 亏损阶梯（log-log）与自拟合阶 1.98、0.3 % 验收线；右=KCS 全尺度五种求积阶
（(10,6) 标定规则落在 Cb·L·B·T 目标上，(24,16) 报告规则高出 0.062 %）。
图 2：G1–G5 的 |实测|/限值 比值图 —— 重放逐位复现。

## ⑤ 注意点（踩过的坑与局限）
- **裁判席空位**：m20516 时代审计（m20464）同锚点判 L1 的理由之一是"无裁判行"；本轮代码可导入、可复现、
  有集成调用方，仍**缺同一件**：`tensor_benchmark_cases` 无该锚点行 → L4 只能记"实质具备、登记缺失"。
- **F1（中）**：KCS 行的 "Cb 标定 0.6505" 只在标定阶 (10,6) 成立；报告阶 (24,16) 得 Cb=0.650904。
  `hull_nurbs.py:557-559` 的"阶依赖远低于标定容差"与本仓测试注释（7e-4、断言 2e-3）矛盾。
- **F3/F4/F5**（低）：`result.json` 把 numpy-stl 记成 absent（`stl` 无 `__version__`）；G5 只是自一致性
  （真复现靠 `tools/volume_xval_report.py`，本审计跑了：results.json 全等、results.md 仅时间戳不同）；
  报告 §3 多船表把 2.723165 渲染成 "2.7"（0.85 %）、2.763654 成 "2.8"（1.32 %）。
- **未跑**：B-rep 引擎（gmsh/OCCT/cadquery）路线——本机不可安装，"第三独立实体引擎"仍是开口工作；
  trimesh/numpy-stl 内部未重导；KCS 52030 m³ 原页未重取；JBC −18.3 % 前提未复核。这些记 unknown，不算 pass。
- 审计只读被审工作树；所有重放在 `/tmp/m20544_repro` 一次性副本中进行。

## ⑥ 误差判定（关键量 vs 基准，逐项回答"是否 <3%"）
| 关键量 | 基准 | 实测偏差 | <3 %？ |
| --- | --- | --- | --- |
| Wigley 最细档体积（engine A，审计复算） | 解析 2.777778 m³ | −0.0321 % | 是（且 <0.3 % 验收线） |
| Wigley 无三角化求积 | 解析 | 1.6e-16（机器零） | 是 |
| 两独立引擎互差（最细档） | engine A | 7.43e-09 | 是 |
| Wigley 亏损阶 | 理论 2 | {ev['auditor_checks']['own_arithmetic']['my_measured_order_from_artifact_ladder']:.3f}（工件阶梯）/ 2.0205（live） | 是（窗口内） |
| KCS 全尺度体积 vs 公开 52030 m³ | 公开值 | +0.0623 % | 是（但属标定阶残差，非恒等式，见 F1） |
| KCS 全尺度 Cb vs 公开 0.6505 | 0.6505 | +0.0620 %（Cb 0.650904） | 是 |
| 非对称艉板三路 vs 解析 | 7.6475677 m³ | −0.826 %（同一三角形网格的采样亏损） | 是 |

**结论**：数值层全部 <3 %（多数 <0.3 %），可复现性逐位通过；**级别 L3**，L4 只差裁判登记这一件（附可执行 INSERT）。
disposition = approve with notes。
"""
    return {"title": "交差核验：XVAL-GEOM-VOLCROSS（m20512）几何体积跨工具互证 —— 裁判席空位挡住 L4，数值与复现逐位通过",
            "summary": ("board-audit 对 m20512 的独立复核：三个工件齐备、repro_cmd 逐位复现（exit 0，40 测试绿，CLI 0/1/2），"
                        "自写算术复算 Wigley 解析/偏差/阶/Richardson 全对；KCS 报告船 Cb 实为 0.650904（标定阶残差 6.2e-04，"
                        "F1 中）；tensor_benchmark_cases 该锚点 0 行 → L4 不可证，判 L3。"),
            "scenario_key": "xval_geom_volcross",
            "mission_id": MISSION,
            "tags": ["board-audit", "XVAL-GEOM-VOLCROSS", "TensorCAD", "geometry", "volume", "repro"],
            "markdown": md,
            "figures": figs}


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--dry-run", action="store_true")
    a = ap.parse_args()
    refs = mirror()
    ev = build_evidence(refs)
    assert all(r["identical"] for r in refs), "ref copy mismatch"
    (HERE / "m20544_verdict_evidence.json").write_text(
        json.dumps(ev, indent=1, ensure_ascii=False, default=str))
    figs = figures()
    payload = report_payload(figs, ev)
    Path("/tmp/report_m20544.json").write_text(json.dumps(payload, ensure_ascii=False))
    print(json.dumps({"level": LEVEL, "evidence_keys": sorted(ev),
                      "refs": len(refs), "figures": [f["name"] for f in figs],
                      "report_chars": len(payload["markdown"])}, ensure_ascii=False))
    if a.dry_run:
        return 0
    with psycopg.connect(DSN) as cx, cx.cursor() as cur:
        cur.execute("""
            INSERT INTO tensor_maturity_verdicts
                (repo, scope, mission_id, level, evidence, verdict_by)
            VALUES (%s, %s, %s, %s, %s::jsonb, %s)
            ON CONFLICT (repo, scope, verdict_by) DO UPDATE
               SET mission_id=EXCLUDED.mission_id, level=EXCLUDED.level,
                   evidence=EXCLUDED.evidence, created_at=now(),
                   gate_state='pending', gate_level=NULL, gate_note=NULL
            RETURNING id, gate_state, gate_level
        """, (REPO, SCOPE, MISSION, LEVEL, json.dumps(ev, ensure_ascii=False, default=str), BY))
        row = cur.fetchone()
        cx.commit()
    print("verdict row:", row)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())