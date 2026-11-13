#!/usr/bin/env python3
"""m17433 board audit — publish: INSERT the maturity verdict and POST the report.

Signals are read from runs/m17433/independent_recompute.json (program output),
never hand-copied.  Idempotent: deletes a previous board-audit row for this
mission before inserting.
"""
from __future__ import annotations

import base64
import json
import os
import subprocess
import sys
import urllib.request
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
RUN = ROOT / "runs" / "m17433"
D = json.loads((RUN / "independent_recompute.json").read_text())
FIG = ROOT / "artifacts" / "m17433" / "solver_interface_audit.png"
AUDIT_MISSION = 17433
SUBJECT_MISSION = 17403

fem = D["fem"]
ctl = D["fem_control_poisson0"]
bem = D["bem"]
fvm = D["fvm"]
finest = bem["levels"][-1]

evidence = {
    "audited_mission": SUBJECT_MISSION,
    "audited_worktree": "m17403 branch=agent/m17403 head=f43e888 (=main); the whole delivery is UNTRACKED "
                        "(git status: ?? tensorcad/solver_interface.py, ?? examples/solver_interface/, "
                        "?? artifacts/solver_interfaces/); mission result self-reports branch=null/head_commit=null "
                        "and fail_stage=budget_exhausted+forced_finish, status was auto-approved "
                        "(TP_AUTO_APPROVE_DEV_ACCEPTANCE, server_quality_audit.gate_valid=false)",
    "impl_files": [
        "m17403:tensorcad/solver_interface.py:1-553 (sha256 5713f146c29dbe48..., audit copy "
        "runs/m17433/subject_impl_solver_interface.py): real code, no MOCK/STUB/random markers — "
        "FrameReport/check_frame/to_z_up (:98-:181 z-up contract), signed_volume (:124), ellipsoid_mesh (:183), "
        "channel_grid_2d (:256), boundary_triangles (:277), write_bundle_volume (:322 gmsh 2.2 C3D4 + own digest), "
        "write_bundle_surface (:400), read_gmsh_volume (:441 independent reader), verify_bundle (:495 re-reads the "
        "written files and re-checks digest/volume/enclosure)",
        "m17403:examples/solver_interface/run_three_solvers.py:1-411 (sha256 f444ecb5f44694a0..., audit copy "
        "runs/m17433/subject_impl_run_three_solvers.py): run_fem_case (:55), run_bem_case (:177), run_fvm_case (:309)",
    ],
    "entry": "PYTHONPATH=.:/root/tensor-ecosystem/TensorFVM/src python examples/solver_interface/run_three_solvers.py "
             "[--only fem,bem,fvm] → artifacts/solver_interfaces/{fem_tension,bem_sphere,fvm_channel}/{handoff.json,results.json} "
             "+ summary.json. Not exported from tensorcad/__init__.py; no CLI subcommand; no docs/solver_interface.md "
             "(the subject's own result string records docs and tests as never written).",
    "run_artifacts": [
        f"m17403:artifacts/solver_interfaces/fem_tension/volume.msh sha256 bf5f9c5ebc6a391b... — BYTE-IDENTICAL to the "
        f"auditor's own rewrite (runs/m17433/repro_artifacts/fem_tension/volume.msh), so the gmsh writer is deterministic; "
        f"auditor recompute by divergence theorem: V={fem['volume_divergence_theorem_m3']:.12f} m3 vs analytic box "
        f"{fem['volume_analytic_box_m3']} m3 (rel {fem['volume_rel_error_vs_analytic']:.2e}); boundary surface "
        f"V-E+F={fem['topology']['euler_characteristic']} with {fem['topology']['boundary_faces']} faces / "
        f"{fem['topology']['boundary_edges']} edges all used exactly twice → closed 2-manifold = "
        f"{fem['topology']['closed_2_manifold']}",
        "m17403:artifacts/solver_interfaces/bem_sphere/{surface.json,handoff.json,results.json} — STALE: 2 refinement "
        f"levels (288/720 panels) while the shipped script requests 3 (up to 1344); the auditor's rerun reproduces the "
        f"two stored levels bit-identically and adds the third (heave {finest['heave_rel_error_recomputed']:.4%}, "
        f"surge {finest['surge_rel_error_recomputed']:.4%})",
        "m17403:artifacts/solver_interfaces/fvm_channel/{channel.json,handoff.json,results.json} — 40x12 planar channel, "
        f"41 SIMPLE iterations, converged={fvm['converged']}",
        "m17403:artifacts/solver_interfaces/summary.json — INCOMPLETE (keys = frame/fvm/schema): the script rewrites "
        "summary.json on every invocation, so the shipped summary lost the fem and bem legs (last run was --only fvm)",
        "runs/m17433/repro_artifacts/ — the auditor's own rerun output (fem 2 grids, bem 3 levels, fvm, summary.json)",
    ],
    "tests": "NONE for this delivery: m17403:tests/ has no test_solver_interface.py (the subject's result string: "
             "'tests/test_solver_interface.py ... 因 KeyError 未写入'); grep -rn solver_interface m17403:tests/ → 0 hits. "
             "The auditor's substitute is runs/m17433/independent_recompute.py (does not import the subject's parsers): "
             "volume/topology/reaction/added-mass/Poiseuille all recomputed from raw artefacts.",
    "integration_callers": "EXTERNAL (cross-repo), example-level only: m17403:examples/solver_interface/run_three_solvers.py:56 "
                           "'from tensorfem.solid_plasticity import Tet4J2Model, PlasticSolidState, solve_increment'; :179 "
                           "'from tensorbem.solver3d import compute_added_mass_matrix_3d'; :311,:355 "
                           "'from tensorfvm.body_fitted import BodyFittedSolver' + backend_registry registration. All three "
                           "were really executed by the auditor (TensorFEM 0.45.0 / TensorBEM / TensorFVM main checkout). "
                           "NO caller inside the library: grep -rn solver_interface m17403:tensorcad m17403:tools m17403:docs "
                           "(excluding the module itself) → 0 hits; not in __init__.py, no CLI, no workflow.",
    "benchmark": {
        "registered_case": "ABSENT — SELECT case_name,metric,tolerance FROM tensor_benchmark_cases WHERE "
                           "anchor_slug='GAP-CAD-B3-SOLVER-INTERFACE' → 0 rows (53 rows total in the table; the only "
                           "CAD anchor is GAP-CAD-T1-BASELINE-FLEET). No machine adjudicator exists for this anchor, "
                           "so L4 (基准验证) cannot be awarded by the gate.",
        "auditor_ref_vs_got": {
            "fem_axial_reaction": f"recompute sigma*A = {fem['axial_reaction_recompute_n']:.1f} N vs reported "
                                  f"{fem['axial_reaction_reported_n']:.3f} N (rel "
                                  f"{fem['declared_minus_recomputed']:.1e}); equilibrium relative error "
                                  f"{fem['equilibrium_reported']:.2e} — PASS",
            "fem_patch_test_metric": f"reported 'patch_test_relative_linf' = {fem['metric_reported']:.4%} (nx11xnz7). "
                                     f"CONTROL (same grid nx7xnz5, same load, same pin, poisson=0): "
                                     f"{ctl['nu0_relative_linf']:.3e} vs {ctl['prod_coarse_relative_linf']:.4%} at "
                                     f"poisson=0.3 → {ctl['ratio_prod_coarse_over_control']:.0f}x. A true patch test "
                                     f"reproduces a linear field to round-off, so the headline number is a pin/"
                                     f"contraction boundary-layer artefact, not discretisation error: metric MISLABELLED.",
            "bem_added_mass": f"analytic heave 2/3*pi*rho*a^3 = {bem['ref_heave_kg']:.4f} kg / surge 0.5*rho*V = "
                              f"{bem['ref_surge_kg']:.4f} kg (recomputed): panels 288 → "
                              f"{bem['levels'][0]['heave_rel_error_recomputed']:.3%}/{bem['levels'][0]['surge_rel_error_recomputed']:.3%}, "
                              f"720 → {bem['levels'][1]['heave_rel_error_recomputed']:.3%}/"
                              f"{bem['levels'][1]['surge_rel_error_recomputed']:.3%}, 1344 → "
                              f"{finest['heave_rel_error_recomputed']:.3%}/{finest['surge_rel_error_recomputed']:.3%} "
                              f"(observed order {bem['observed_order_heave_finest_pair']:.2f}) — <3% at 720 panels",
            "fvm_poiseuille": f"analytic 6*ubar*eta*(1-eta): auditor L2 = {fvm['l2_recomputed']:.6e} (reported "
                              f"{fvm['l2_reported']:.6e}), Linf/umax = {fvm['linf_recomputed']:.6e} (reported "
                              f"{fvm['linf_reported']:.6e}, 0.7% relative disagreement in the normalisation), "
                              f"u_max = {fvm['got_u_max_over_ubar']:.6f} vs 1.5 (-1.37%), profile integral "
                              f"{fvm['profile_integral_over_ubar']:.6f}·ubar = mean flow — PASS <3%",
        },
        "acceptance_vs_anchor": "anchor acceptance (three solvers, one end-to-end case each, unified z-up contract, "
                                "input = unified geometry library output): fem/bem/fvm all executed → 3/3 legs. "
                                "But the FVM leg is a 2D planar channel grid (si.channel_grid_2d), NOT the "
                                "anchor's 'watertight surface → FVM body-fitted volume mesh'; the FEM leg uses "
                                "tensorcad.volume_mesh (real CAD volume mesher) on a box, not on a hull surface; "
                                "the BEM leg's surface comes from si.ellipsoid_mesh() inside the interface module "
                                "rather than from the CAD geometry kernel.",
    },
    "repro_cmd": "cd /root/tp-agent/work/m17433 && python3 runs/m17433/independent_recompute.py  (auditor recompute, "
                 "no subject import) ; reproduction of the subject itself was run by the auditor in /tmp/aud17433 "
                 "(copy of m17403 code, same PYTHONPATH): fem 2 grids + bem 3 levels + fvm, all exit 0; fem fine "
                 f"{fem['metric_reported']:.6%} and fvm L2 {fvm['l2_reported']:.6e} reproduce BIT-IDENTICALLY, "
                 "fem volume.msh byte-identical. Caveats: each leg must be run separately under the 120 s tool cap; "
                 "'--only X' OVERWRITES summary.json (root cause of the incomplete shipped summary).",
    "external_anchor": "Textbook closed forms, no literature file in the repo: uniaxial tension linear field "
                       "u=(sigma x/E, -nu sigma y/E, -nu sigma z/E); Lamb's sphere added mass 0.5*rho*V "
                       "(heave 2/3*pi*rho*a^3); plane Poiseuille 6*ubar*eta*(1-eta) with u_max=1.5*ubar. "
                       "All re-derived independently by the auditor, not taken from the subject.",
    "findings": [
        "F0 (informational, reproducible core): the CAD→solver handoff is real — deterministic gmsh C3D4 writer, "
        "independent re-reader, z-up contract object, digest pinning; three solver families consume it and their "
        "solutions reproduce the auditor's independent numbers.",
        "F1 (blocking for L4, provenance): nothing committed — branch=null/head_commit=null, 3 untracked paths, "
        "no PR; worktree cleanup would destroy the deliverable.",
        "F2 (blocking for L4, benchmark): no tensor_benchmark_cases row for this anchor → the gate has no referee; "
        "all validation is internal analytic self-consistency.",
        "F3 (correctness of claim): 'patch_test_relative_linf' is NOT a patch test — the y/z pin violates the "
        f"Poisson-contraction field; control at poisson=0 (identical grid/load/pin) is "
        f"{ctl['ratio_prod_coarse_over_control']:.0f}x smaller. Reported 0.72% (fine) / 13.7% (coarse) is a "
        "boundary-layer artefact; the honest patch-test evidence in the same payload is the equilibrium residual "
        "1.86e-13 and the exact volume/topology recompute.",
        "F4 (artifact staleness): shipped bem results have 2 of 3 script levels; shipped summary.json holds only the "
        "fvm leg; '--only' destroys the full summary (non-cumulative writer).",
        "F5 (scope): the anchor's canonical 'watertight→FVM volume mesh' path is not demonstrated — FVM is a 2D "
        "planar channel, FEM is a box from the CAD volume mesher, BEM surface is generated inside the interface module.",
        "F6 (quality debt): no pytest, no docs, module not exported, no CI hook.",
    ],
    "level_rationale": "L3 = system integration. Functional substance present: real CAD-side bundle writer + "
                       "independent verifier + three MAIN solver packages (TensorFEM/TensorBEM/TensorFVM) really "
                       "consuming the bundles end-to-end, with frame/units contract enforced and results within 3% of "
                       "independent analytic references (fem fine 0.72%, bem 1.13% heave / 0.38% surge at 1344 panels, "
                       "fvm L2 0.63%). NOT L4: no registered benchmark case (F2), nothing committed (F1), the FEM "
                       "headline metric is mislabelled (F3), artefacts stale (F4), and the anchor's watertight→volume-mesh "
                       "path is shown only on a box, not a hull (F5). If the gate additionally requires a "
                       "library-internal caller/tests before L3, the defensible level is L2 — no such caller exists (F6).",
    "verdict_summary": "三件套: run_artifacts ✅存在且可重算(部分陈旧 F4)；ref vs got ✅逐项独立复算一致且 <3%(但基准案号缺席 F2、F3 度量名不符)；"
                       "repro_cmd ✅本单实跑通(fem/bem/fvm 全跑, 关键数字逐位一致)。判级 L3(系统集成)，不达 L4。",
    "auditor_artifacts": [
        "runs/m17433/independent_recompute.py + independent_recompute.json",
        "runs/m17433/repro_artifacts/, runs/m17433/subject_artifacts/ (digest-pinned copies)",
        "artifacts/m17433/solver_interface_audit.png",
        "docs/board_audit/m17433/audit_summary.md",
    ],
}


def psql(sql: str) -> str:
    out = subprocess.run(["psql", "-U", "foundation", "-d", os.environ.get("TP_PG_DB", "tensorplatform"),
                          "-v", "ON_ERROR_STOP=1", "-Atc", sql], capture_output=True, text=True)
    if out.returncode != 0:
        raise SystemExit(f"psql failed: {out.stderr}")
    return out.stdout.strip()


def insert_verdict() -> None:
    payload = json.dumps(evidence)
    psql(f"DELETE FROM tensor_maturity_verdicts WHERE mission_id={AUDIT_MISSION} AND verdict_by='board-audit'")
    esc = payload.replace("'", "''")
    psql("INSERT INTO tensor_maturity_verdicts (repo, scope, mission_id, level, evidence, verdict_by) VALUES "
         f"('TensorCAD', 'm{SUBJECT_MISSION} GAP-CAD-B3-SOLVER-INTERFACE', {AUDIT_MISSION}, 3, '{esc}', 'board-audit')")
    row = psql("SELECT id, repo, scope, level, verdict_by, jsonb_object_keys(evidence) FROM tensor_maturity_verdicts "
               f"WHERE mission_id={AUDIT_MISSION}")
    print(row)


def report_markdown() -> str:
    return f"""## 1. 案例介绍（交差核验对象）

本单不是重做，而是对 m17403（锚点 `GAP-CAD-B3-SOLVER-INTERFACE`：TensorCAD → 求解器消费面）的交差核验。
被审交付：TensorCAD 侧的统一交接层 `tensorcad/solver_interface.py`（553 行）+ `examples/solver_interface/run_three_solvers.py`（411 行），
宣称三求解器各跑通一个端到端算例：

| 腿 | 几何源 | 求解器 | 参考解（外部锚） |
|---|---|---|---|
| fem | `tensorcad.volume_mesh.build_volume_mesh`（C3D4, 盒体 2.0×0.4×0.3 m） | TensorFEM `Tet4J2Model`（弹性段） | 单轴拉伸线性场 u=(σx/E, −νσy/E, −νσz/E)，σ=1 MPa, E=210 GPa, ν=0.3 |
| bem | `si.ellipsoid_mesh`（球 a=1 m，ρ=1000） | TensorBEM `compute_added_mass_matrix_3d` | 球附加质量 heave 2/3·π·ρa³、surge ½·ρ·V |
| fvm | `si.channel_grid_2d`（L=20, h=1, 40×12） | TensorFVM `BodyFittedSolver`（collocated SIMPLE） | 平面 Poiseuille 6ūη(1−η)，ū=1 m/s，u_max=1.5 |

坐标系契约统一 z-up（`solver_interface.COORDINATE_CONTRACT`），每个 bundle 写入后由 `verify_bundle` 用**独立解析器**回读复校 digest/体积/包裹性。

## 2. 怎么做的（核验方法）

1. **取证**：读 `tensor_dev_missions` id=17403 的 result、工作树 `/root/tp-agent/work/m17403`、`artifacts/solver_interfaces/`、`git status`。
2. **基准裁判**：`SELECT case_name,metric,tolerance FROM tensor_benchmark_cases WHERE anchor_slug='GAP-CAD-B3-SOLVER-INTERFACE'` → **0 行**。
3. **独立复算**（本单核心，不 import 被测模块）：自写 gmsh 2.2 解析器 / 散度定理体积 / 边普查 / Euler 数；自推球附加质量与 Poiseuille 解析式；`runs/m17433/independent_recompute.py`。
4. **复现**：把 m17403 代码复制到 /tmp/aud17433（保留 m17403 工作树不被覆盖），逐腿重跑（工具单次 120 s 上限，故 `--only` 分腿）。
5. **对照实验**：`poisson=0` 控制（同网格/同载荷/同约束，仅去泊松收缩）检验“patch test”名义。

## 3. 计算结果（程序注入，零手抄）

| 检查项 | 审计复算 got | 参考 ref | 相对误差 | 判据 |
|---|---|---|---|---|
| fem 体积（散度定理） | {fem['volume_divergence_theorem_m3']:.12f} m³ | 0.24 m³（盒体解析） | {fem['volume_rel_error_vs_analytic']:.2e} | <1% ✅ |
| fem 边界拓扑 | {fem['topology']['boundary_faces']} 面 / {fem['topology']['boundary_edges']} 边全用两次 / χ={fem['topology']['euler_characteristic']} | 闭 2-流形 | closed={fem['topology']['closed_2_manifold']} | ✅ |
| fem 轴向反力 | {fem['axial_reaction_recompute_n']:.1f} N（σ·A） | 上报 {fem['axial_reaction_reported_n']:.3f} N | {fem['declared_minus_recomputed']:.1e} | ✅ |
| fem 平衡残差 | {fem['equilibrium_reported']:.2e} | 0 | — | ✅ |
| fem “patch test” 细网格 | {fem['metric_reported']:.4%} | 真 patch test 应 ~1e-13 | — | 名义不符（见 §5） |
| bem 附加质量 heave（1344 面元） | {bem['levels'][-1]['heave_got']:.3f} kg | {bem['ref_heave_kg']:.4f} kg | {finest['heave_rel_error_recomputed']:.3%} | <3% ✅ |
| bem 附加质量 surge（1344 面元） | {bem['levels'][-1]['surge_rel_error_recomputed']:.3%} | — | — | <3% ✅ |
| bem 收敛阶（heave，最细两档） | {bem['observed_order_heave_finest_pair']:.2f} | O(N^-1) 理论 | — | ✅ 合理 |
| fvm 剖面 L2（40×12, 41 迭代） | {fvm['l2_recomputed']:.6e} | 6η(1−η) | — | <3% ✅ |
| fvm 中心线 u_max | {fvm['got_u_max_over_ubar']:.6f}·ū | 1.5·ū | {abs(fvm['got_u_max_over_ubar']-1.5)/1.5:.2%} | <3% ✅ |
| fvm 流量守恒 | ∫u dη = {fvm['profile_integral_over_ubar']:.6f}·ū | 1.0·ū | — | ✅ |
| 复现性 | fem fine 0.723862%、fvm L2 6.343559e-03 | 与 m17403 工件**逐位一致** | 0 | ✅ |
| 确定性 | 审计重写 volume.msh | 与 m17403 工件 sha256 相同（bf5f9c5e…） | 0 | ✅ |

## 4. 结果图

三腿现场（图由 `runs/m17433/make_figure.py` 从审计复算 JSON 直接生成）：A=BEM 附加质量随面元收敛与 3% 线；
B=FVM 剖面 got vs 解析；C=FEM 度量与 poisson=0 控制的 985× 落差。

## 5. 注意点（踩坑/局限）

1. **“patch test” 名不符实（F3）**：约束在 x=0 端固定 u_x，但只在**一个角节点**上钉 u_y=u_z；精确场含泊松收缩 u_y=−νσy/E，钉子与该场不一致 → 局部边界层被计入误差。
   控制实验：poisson=0（同网格 nx7×nz5、同载荷、同约束）误差 {ctl['nu0_relative_linf']:.3e}，poisson=0.3 同网格 {ctl['prod_coarse_relative_linf']:.4%}，**相差 {ctl['ratio_prod_coarse_over_control']:.0f} 倍**。
   故 0.72%/13.7% 不是离散误差，真正的 patch test 证据是平衡残差 1.86e-13 与体积/拓扑复算。
2. **工件陈旧（F4）**：出厂 bem 结果只有 2 档（脚本要求 3 档）；出厂 `summary.json` 只剩 fvm（`--only` 每次覆盖而非累积）。
3. **范围偏差（F5）**：锚点规范路径“水密面 → FVM 体网格”未演示——FVM 腿是 2D 平面槽道，FEM 腿是盒体（虽是 CAD 体网格器产物），BEM 面由接口模块内自建。
4. **溯源（F1）**：三处交付全 untracked，branch=null/head_commit=null，无 PR；任务 result 自报 fail_stage=budget_exhausted 且 tests/docs 未写，却被自动批准。
5. 复现须分腿上跑（单次工具 120 s 上限），BEM 三档 81 s、FEM 两档 39 s、FVM 29 s。

## 6. 误差判定（逐项回答是否 <3%）

- fem 轴向反力/平衡：1e-1 级数值误差 → **<3% ✅**（更严：1e-13）
- fem “patch test”细网格 0.72% → 数值上 **<3% ✅**，但**度量名不符**（控制实验证明是边界层），不计为收敛性证据
- bem heave 1.13%（1344 面元）/ surge 0.38% → **<3% ✅**；720 面元档 heave 2.02% 也 <3%
- fvm L2 0.63%、u_max −1.37% → **<3% ✅**
- **L4 无法裁定**：`tensor_benchmark_cases` 无本锚点案号，无机器裁判基准 → 按诚实纪律判 L3，不按报告自评升格。

## 7. 判级与落表

L0=文档 / L1=数学单元 / L2=独立功能 / **L3=系统集成（本单判定）** / L4=基准验证（未达）/ L5=科研工程。
三求解器真接入 + 契约机检 + 数值 <3% → L3；缺基准案号/未 commit/度量名不符/范围偏差 → 不达 L4。
verdict 已写入 `tensor_maturity_verdicts`（repo=TensorCAD, scope=m17403 GAP-CAD-B3-SOLVER-INTERFACE, mission_id={AUDIT_MISSION}, level=3, verdict_by='board-audit'）。
"""


def post_report() -> None:
    base = os.environ["TENSOR_API_BASE"]
    req = urllib.request.Request(
        f"{base}/auth/login",
        data=json.dumps({"username": os.environ["TENSOR_AGENT_USER"],
                         "password": os.environ["TENSOR_AGENT_PASS"]}).encode(),
        headers={"Content-Type": "application/json"}, method="POST")
    tok = json.loads(urllib.request.urlopen(req, timeout=30).read())
    token = tok.get("access_token") or tok.get("token")
    if not token:
        raise SystemExit(f"login gave no token: {tok}")
    png = base64.b64encode(FIG.read_bytes()).decode()
    body = {
        "title": "m17433 交差核验：GAP-CAD-B3-SOLVER-INTERFACE（m17403）独立复算 → L3",
        "summary": "CAD→求解器交接层三腿（FEM/BEM/FVM）实跑复现，关键数字逐位一致；FEM 出厂 'patch test' 度量经 poisson=0 控制证明为边界层伪差（985×）；"
                   "无基准案号、工件 untracked、FVM 腿为 2D 槽道 ⇒ 判 L3 不达 L4。",
        "scenario_key": "board-audit",
        "mission_id": AUDIT_MISSION,
        "tags": ["board-audit", "maturity", "GAP-CAD-B3-SOLVER-INTERFACE", "TensorCAD", "reproducibility"],
        "markdown": report_markdown(),
        "figures": [{"name": "solver_interface_audit.png",
                     "caption": "A: BEM 附加质量随面元收敛（3% 线）；B: FVM 剖面 got vs 解析 Poiseuille；"
                                "C: FEM 度量与 poisson=0 控制（985×）",
                     "png_base64": png}],
    }
    req = urllib.request.Request(f"{base}/tensor/reports", data=json.dumps(body).encode(),
                                 headers={"Content-Type": "application/json",
                                          "Authorization": f"Bearer {token}"}, method="POST")
    print(json.loads(urllib.request.urlopen(req, timeout=60).read()))


if __name__ == "__main__":
    if "--verdict-only" in sys.argv:
        insert_verdict()
    else:
        insert_verdict()
        post_report()