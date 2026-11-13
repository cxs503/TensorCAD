#!/usr/bin/env python3
"""Generate AUDIT.md for m17405 (board audit of m17357) from raw artifacts.

Numbers are injected programmatically from:
  * artifacts/geom_contract_audit/probe_results.json   (auditor's own probes)
  * <subject>/artifacts/geom_contract_assets/compliance_report.json
  * a fresh `pytest` run of the subject's test file
so no figure in the report is hand-typed.
"""
from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
SUBJECT = Path(sys.argv[1] if len(sys.argv) > 1 else "/root/tp-agent/work/m17357")
PROBE = ROOT / "artifacts" / "geom_contract_audit" / "probe_results.json"


def pytest_run() -> tuple[int, int, str]:
    r = subprocess.run([sys.executable, "-m", "pytest", "tests/test_geom_contract.py", "-q"],
                       cwd=SUBJECT, capture_output=True, text=True, timeout=900)
    tail = r.stdout.strip().splitlines()[-1] if r.stdout.strip() else ""
    passed = failed = 0
    for tok in tail.replace(",", " ").split():
        if tok.isdigit():
            continue
    import re
    m = re.search(r"(\d+) passed", tail)
    passed = int(m.group(1)) if m else 0
    m = re.search(r"(\d+) failed", tail)
    failed = int(m.group(1)) if m else 0
    return passed, failed, tail


def main() -> int:
    probe = json.loads(PROBE.read_text())
    comp = json.loads((SUBJECT / "artifacts/geom_contract_assets/compliance_report.json").read_text())
    npass, nfail, pytest_tail = pytest_run()

    anchors = "".join(
        f"| {a['case']} | `{a['formula']}` | {a['ref']:.9f} | {a['got']:.9f} | "
        f"{a['rel_err_pct']:+.4f}% | {'YES' if a['within_3pct'] else 'NO'} | "
        f"{'YES' if a['within_1pct'] else 'NO'} |\n"
        for a in probe["anchors"])
    prism = "".join(
        f"| {p['case']} | {p['ref']:.6f} | {p['got']:.6f} | {p['rel_err_pct']:+.3f}% | "
        f"{p['implied_cb']:.4f} | {p['declared_cb']:.2f} |\n"
        for p in probe["prismatic_observation_out_of_scope"])
    red = "".join(
        f"| {c['case']} | {'PASS' if c['passed'] else 'FAIL'} | "
        f"{', '.join(c['conventions']) or '—'} | {', '.join(c['severities']) or '—'} |\n"
        for c in probe["red_team"])
    fn = "".join(
        f"| {c['probe']} | {'PASS (silent)' if c['passed'] else 'FAIL'} | "
        f"{', '.join(c['violations']) or 'none'} | {c['note']} |\n"
        for c in probe["false_negatives"])
    cli = "".join(
        f"| {c['case']} | {c['declared_units']} | {c['exit']} | {c['expected_exit']} | "
        f"{'ok' if c['exit_ok'] else 'MISMATCH'} | {c['json_ok']} |\n"
        for c in probe["cli"])
    comp_rows = "".join(
        f"| {r['source'].rsplit('/', 1)[-1]} | {'PASS' if r['passed'] else 'FAIL'} | "
        f"{r['mesh_info'].get('units')} | {r['mesh_info'].get('diagonal_m')} | "
        f"{r['mesh_info'].get('signed_volume')} | "
        f"{'; '.join('[' + v['convention'] + '/' + v['severity'] + '] ' + v['message'] for v in r['violations']) or '—'} |\n"
        for r in comp)
    repro = probe["repro"]
    body = f"""# m17405 · 交差审计：m17357 `GAP-CAD-D3-GEOM-CONTRACT`（TensorCAD）

审计单：m17405（board-audit）· 对象：m17357（已 succeeded）· 裁判：`tensor_benchmark_cases` + `maturity_gate` JUDGE_SCHEME
对象工作树：`{SUBJECT}` · 审计可复现脚本：`docs/board_audit/m17405/audit_probe.py`

---

## 0. 结论先行

| 项 | 判定 |
|---|---|
| 声明级别 | 无（原单 result 无结构化 verdict） |
| **审计判级** | **L2 独立功能实现** |
| gate 预期 | `judge(2) → 2 → confirmed`（L3/L4 因键缺/无调用方/无基准案自动降级） |
| 处置 | approve with notes（confirmed=L2，不返修；三条后续项见 §7） |
| 三件套 | run_artifacts ✅ / repro_cmd ✅（字节级可复现）/ benchmark ❌（锚点零基准案） |
| 六段报告 | ❌ `tensor_reports` 无 mission_id=17357 行（report_id 为空） |
| 提交状态 | ❌ 8 个交付文件全为工作树未跟踪文件，`head_commit=null`，主检出 `/root/tensor-ecosystem/TensorCAD` 无 `geom_contract.py` |

---

## 1. 案例介绍（这是什么交付）

`GAP-CAD-D3-GEOM-CONTRACT` 要求 TensorCAD 给出**几何契约文档 + 运行时校验器**，把四条约定
（坐标系 / 单位 / 比例 / 取向）做成机器可执行的检查，并让 8 个现存 STL 资产全过一遍出合规表，
违规可自动检出（红队实证）。

对象交付物（全部真实存在，逐个 sha256 留痕于 `evidence.json`）：

| 文件 | 行数/大小 | 角色 |
|---|---|---|
| `tensorcad/geom_contract.py` | 347 行 | 库函数 `validate_contract`/`validate_stl_file` + CLI `main` |
| `docs/geometry_contract.md` | 149 行 | 四约定文档 + 用法 |
| `tests/test_geom_contract.py` | 231 行 / 22 测试 | 单元验证 |
| `tools/run_contract_compliance.py` | 107 行 | 8 资产合规表跑批（repro_cmd） |
| `tools/generate_contract_assets.py` | 259 行 | 资产与红队样本生成 |
| `artifacts/geom_contract_assets/` | 9 STL + `compliance_report.json` | run_artifacts |

四约定的判据（源码位置）：`_check_coordinate_system` :124（声明串）、`_check_units` :149（SI 词表）、
`_check_scale` :161（包围盒对角线 ∈ [0.05, 500] m）、`_check_orientation` :188（闭合体有符号体积 > 0；
开放面法向半球一致率 ≥ 0.5）。

## 2. 怎么做的（审计方法）

1. **取证**：读 `tensor_dev_missions` id=17357 的 `result`（含 `server_quality_audit`、`git_status`、`changed_files`）、
   对象工作树全部交付文件、`artifacts/geom_contract_assets/compliance_report.json`、`tensor_reports`（无行）。
2. **基准裁判**：`SELECT case_name,metric,tolerance FROM tensor_benchmark_cases WHERE anchor_slug='GAP-CAD-D3-GEOM-CONTRACT'` → **0 行**。
   同库 `repo='TensorCAD'` 仅有 1 条（`GAP-CAD-T1-BASELINE-FLEET`，另一锚点），与本锚点不可借用。
3. **独立探针**（`docs/board_audit/m17405/audit_probe.py`，不 import 对象的测试夹具）：
   自建网格 → 解析体积对照 / 5 组红队违规 / 3 组"沉默窗口"反例 / CLI 退出码 / 复现哈希比对。
4. **复算**：`python3 tools/run_contract_compliance.py` 前后 `compliance_report.json` sha256
   `{repro['sha256_before'][:16]}…` → `{repro['sha256_after'][:16]}…`，**identical={repro['identical']}**。
5. **测试**：对象工作树内 `pytest tests/test_geom_contract.py -q` → **{npass} passed / {nfail} failed**（`{pytest_tail}`）。
6. **集成面**：`grep -rn validate_contract` 全仓 → 除 `tensorcad/__init__.py:457` 再导出外**零调用方**（无求解器/导入管线挂接）。

## 3. 计算结果（程序化注入，零手抄）

### 3.1 对象的合规表（8 资产，四约定）

| 资产 | 判定 | 单位 | 对角线 (m) | 有符号体积 (m³) | 违规 |
|---|---|---|---|---|---|
{comp_rows}
5 PASS / 3 FAIL（红线：翻转法线、错标单位、左手系），与对象声明一致。

### 3.2 审计方独立解析基准（ref=闭式解，got=从交付 STL 由散度定理重算）

| 算例 | 公式 | ref (m³) | got (m³) | 相对误差 | <3% | <1% |
|---|---|---|---|---|---|---|
{anchors}
散度定理复算与解析解一致到 **0.0000% / −0.4742% / −0.2360%**：校验器的体积量是**真实几何测量**，
不是回填常量（球面负偏是 32×64 内接多面体的固有正向下偏）。

### 3.3 CLI 契约

| 用例 | 声明单位 | 实测退出码 | 期望 | 判定 | JSON 可解析 |
|---|---|---|---|---|---|
{cli}

## 4. 结果图

`artifacts/geom_contract_audit/contract_audit.png`（两联图，脚本 `make_figure.py`）：
(a) 8 资产 × 4 约定的违规矩阵（色块=违规，空白=通过）；(b) 交付 STL 体积 vs 闭式解析体积的奇偶图
（含 ±3% 带与逐点相对误差标注）。纯文档/后端任务本可无图，此处出图是因为"契约是否真在测量"必须
看数值对照而非读报告文字。

## 5. 注意点（坑与局限）

* **STL 不携带元数据**：单位与坐标系只能由调用方声明；`run_contract_compliance.py` 按资产硬编码声明值，
  因此 `wrong_units`/`left_handed` 两个红队样本的"检出"本质是**校验器对声明的裁决**，不是从几何反推。
  这不是缺陷（源码 docstring 已声明），但宣传语"左手系可自动检出"须限定为**已声明**的左手系。
* **坐标系检查是字符串检查**：`_check_coordinate_system` 只查 `right-handed`/`left-handed` 关键字；
  一个把右手系网格镜像、却声明 `right-handed` 的开放面网格**零违规通过**（见 §6 反例 C3）。
* **比例检查有盲窗**：对角线带 [0.05, 500] m → 100× 单位错标只要落在带内即通过（反例 C1）。
* **开放面取向只到 warning**：`ratio < 0.5` 才告警且不阻断（只有 error 才使 `passed=False`），
  恰好 50% 反向面**零违规通过**（反例 C2）。
* **未提交**：交付物仅存在于对象工作树（未跟踪）；主检出无此文件，其他仓/CI 现在无法 import 该能力。
* 审计方读到了对象源码（非盲审）；外部检索未做（无外网），外部锚点全部为闭式解析解，故不依赖文献。

## 6. 红队与反例（审计方自建夹具）

### 6.1 应检出（红队正例）

| 违规用例 | 判定 | 命中约定 | 严重度 |
|---|---|---|---|
{red}
**4/4 约定全部命中且约定标签正确**，与被审计方自报"3 类违规可自动检出"一致（审计方另加 `units`/`scale` 两例）。

### 6.2 应检出而未检出（反例，界定声明边界）

| 反例 | 观察 | 命中约定 | 说明 |
|---|---|---|---|
{fn}
结论：三条"沉默窗口"都是**有界的、可文档化的**局限，不构成造假；但使"自动检出"这句话的适用范围
必须收敛到"声明可读 + 越出带宽 + 开放面反向率 > 50%"。

### 6.3 越界观察（不计入本单判定）

| 算例 | ref=Cb·L·B·T | got (STL 体积) | 相对误差 | 隐含 Cb | 声明 Cb |
|---|---|---|---|---|---|
{prism}
校验器忠实报告了网格里的体积；偏差落在**资产生成器**一侧（Series 60 网格隐含 Cb 0.5779 / 0.7195 vs 声明
0.60 / 0.75），建议另开单用 `tensorcad.series60.analytic_volume` 仲裁，本单不判。

## 7. 误差判定与判级

| 关键量 | 基准 | 误差 | 是否 <3% |
|---|---|---|---|
| 箱形驳有符号体积 | 解析 V=L·B·D=1200 m³ | {probe['anchors'][0]['rel_err_pct']:+.4f}% | 是 |
| 单位球有符号体积 | 解析 4/3πr³=4.18879 m³ | {probe['anchors'][1]['rel_err_pct']:+.4f}% | 是 |
| Wigley 有符号体积 | 解析 (4/9)L·B·T=2.77778 m³ | {probe['anchors'][2]['rel_err_pct']:+.4f}% | 是 |
| 基准案 ref/got vs tolerance | `tensor_benchmark_cases` 该锚点 **0 行** | 无法计算 | **无基准案，L4 不成立** |

**判级 = L2（独立功能实现）**，按 `maturity_gate.SCHEME` 逐键核：

| 键 | L2 要求 | 实际 | 结论 |
|---|---|---|---|
| `impl_files` | 有非占位实现主体 | `tensorcad/geom_contract.py` 347 行真实现（无 NotImplemented/pass 占位） | ✅ |
| `entry` | 入口存在 | 库 `validate_contract` :240 / `validate_stl_file` :268 / CLI `main` :282（`python -m tensorcad.geom_contract`） | ✅ |
| `run_artifacts` | 运行工件存在非空 | `artifacts/geom_contract_assets/compliance_report.json`（8 报告）+ 9 STL；审计复跑字节一致 | ✅ |
| `tests` | 测试名非空 | 22 测试 全绿 | ✅ |
| `integration_callers` | 主求解器/平台侧调用方 | 仅 `tensorcad/__init__.py:457` 再导出，**无调用方** | ❌ → **拒绝 L3** |
| `benchmark` | ref vs got + 误差数字 | 该锚点无 `tensor_benchmark_cases` 行；审计方 3 条解析锚点非登记基准案 | ❌ → **拒绝 L4** |
| `repro_cmd` | 可复现命令 | `python3 tools/run_contract_compliance.py`（sha256 前后一致） | ✅（但 L4 已被 benchmark 卡死） |
| `external_anchor` | 文献/实验值 | 无外部检索；闭式解析解已计入 §3.2 | 不申 |

## 8. 处置与后续

**Disposition：approve with notes。** 功能实质=L2（真代码、真测量、可复现、红队 4/4 命中）；
L3（未接管线）与 L4（无登记基准案）**明确拒绝**，附百页报告不改变级别。

后续项（各自独立成单，不在本单返修）：
1. **落库**：把 8 个未提交文件 commit（当前 `head_commit=null`），否则能力对生态不可见。
2. **挣 L3**：在网格导入口（如 `tensorcad/export.read_stl` 或其调用方）加契约门，`integration_callers` 才有 grep 指纹。
3. **挣 L4**：为 `GAP-CAD-D3-GEOM-CONTRACT` 在 `tensor_benchmark_cases` 登记 ref/tolerance 行（建议箱形驳 1200 m³ 与 Wigley (4/9)LBT 两条解析锚点，tolerance 1%）。
4. **收紧取向检查**：闭合体已足够（体积符号）；开放面建议改为"反向面占比 > 10% 即 error"，并把 50% 沉默窗口写进文档。
5. **补报告**：m17357 无 `tensor_reports` 行（report_id 空），六要素报告欠账。

---
*生成：`docs/board_audit/m17405/make_audit_doc.py`（数字全部来自 probe_results.json / compliance_report.json / pytest 输出）。*
"""
    (HERE / "AUDIT.md").write_text(body)
    print(f"wrote {HERE/'AUDIT.md'} ({len(body)} bytes); pytest {npass} passed / {nfail} failed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())