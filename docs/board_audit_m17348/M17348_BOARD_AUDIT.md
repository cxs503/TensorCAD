# m17348 交差核验报告 — GAP-CAD-A3-CURVATURE-FAIRNESS（曲率场与光顺性诊断）

| 项 | 值 |
|---|---|
| 审计单 | m17400（board-audit，交差核验·非重做） |
| 被审单 | m17348（TensorCAD，status=succeeded，claimed_by=tp-impl-7，session 12573） |
| 锚点 | `GAP-CAD-A3-CURVATURE-FAIRNESS` |
| 被审工作树 | `/root/tp-agent/work/m17348`，分支 `agent/m17348`，HEAD `f43e888`（"Merge #398 (m17333) via push-merge"） |
| 被审报告 | `tensor_reports` id=16354（5066 字符，六段齐全，图 sha256 `d110827…45e3f4`） |
| 审计复现脚本 | `scripts/audit/m17348_board_audit.py`（一键重跑本报告全部数字） |
| 审计原始证据 | `artifacts/audit_m17348/{audit.json, verdict_evidence.json, sha256_manifest.json, delivered/, logs/, audit_figure.png, report_figure.png}` |
| 裁定 | **L4**（光顺能量/积分核）＋逐点曲率场 L2 级；交付态阻塞项 F1；落表 `tensor_maturity_verdicts`（verdict_by=board-audit） |

> 独立性限制：被审单与本单同为 agent 执行，审计者未参与 m17348 的实现；全部数字由本单脚本在**被审工作树的字节副本**上重跑得到，未采信报告文本。复核口径＝`tensor_benchmark_cases` 登记判据 + 解析闭式 + 本单重跑。

---

## 0. 一句话结论

被审单声称的**三件套真实存在且可复现**：10/10 交付文件落地、43 测试全绿、11/11 runner claim 绿、`results.json` 重跑逐叶对比仅 1 处差异（绝对路径字符串）、`curvature_crosscheck.json` 重跑 **sha256 逐字节相同**。**但**：① 交付物**不在任何 git ref**（分支上 0 字节），平台若按分支推送，本能力将整体丢失；② 验收头号条款「球主曲率 1/R <1e-6」测的是**既有解析场路径**，不是本单新建的 mesh 曲率场，新 mesh 场逐点精度实测 1.6% / 12.8% / 40%（中位），跨过 15% 精度线（被审报告 §5/§6 已如实披露，未粉饰）。

---

## 1. 三件套取证（逐项，指针到行/工件）

### 1.1 run_artifacts —— 真实存在，可重跑，可逐字节复核

| 工件 | bytes | sha256(前16) | 重跑结果 |
|---|---|---|---|
| `docs/fairness_energy/results.json` | 11562 | `d1b3cc77ffc1d4c5` | 重跑后逐叶对比 **1 处差异**（`/figure` 绝对路径），数值叶全部 ≤1e-12 |
| `docs/fairness_energy/report.md` | 3089 | `89dc41b601747565` | runner 自动重写，非手抄 |
| `docs/fairness_energy/plot.png` | 97370 | `d11082784201f83d` | 与报告 16354 内嵌图 **sha256 完全一致** |
| `artifacts/curvature_crosscheck.json` | 13581 | `88079e170c74e367` | 重跑后 **sha256 完全一致**（bit-identical） |
| `tensorcad/fairness_energy.py` | 24854 | `b5582bfd48192c5e` | 585 行：`MeshCurvatureField:100`、`principal_curvature_field:158`、`smoothing_energy:402`、`fairness_energy_report:554` |
| `tensorcad/wigley.py`（改） | 57803 | `7bde6eac8666b508` | `git diff` = **209 insertions, 0 deletions**（纯增量，未改动既有函数） |
| `tests/test_fairness_energy.py` | 16368 | `eb954a56dccbed58` | 20 pins；与原 `test_curvature.py` 合计 **43 passed**（重跑 43.4 s） |

字节副本已落 `artifacts/audit_m17348/delivered/`（含 `sha256_manifest.json`），即使被审工作树被回收，证据链不丢。

### 1.2 benchmark —— ref vs got（数值表，全部由本单脚本从 `results.json` 重算）

**A. 解析球主曲率（闭合体，1449 探针 = 7 半径 × 207）**

| 量 | ref | got | 误差 | 判据 |
|---|---|---|---|---|
| `max abs \|k−1/R\|` | 1/R，(R=0.1…10 m) | `principal_curvatures(field, probes)` | **2.248e-07 1/m**（相对 2.25e-08） | 任务榜 <1e-6 ✔ |

⚠ **该行测的是既有解析隐式场求值器**（`tensorcad/curvature.py`），不是本单新建的逐顶点 mesh 曲率场（见 §3 F2）。

**B. 球光顺能量（闭式 8π=25.132741、4π=12.566371）**

| 路径 | got | 相对误差 | 判据 |
|---|---|---|---|
| field 路径（level 32，28032 faces） | 25.146062 | **+0.0530%** | 自身 pin ≤0.1% ✔ |
| conditioned mesh 路径（level 32） | 26.573625 | **+5.7330%** | 自身 pin ≤6% ✔ |
| 恒等式 k1²+k2²=4H²−2K 残差 | — | 4.44e-16 | 内部一致性 ✔ |

**C. Wigley 光顺能量 vs 解析船侧（ref 21.537877，湿面积 14.879063 m²）**

| nx×nz | faces | smooth E | crease E | vs 解析 |
|---|---|---|---|---|
| 21×7 | 518 | 12.9265 | 157.48 | −39.98% |
| 41×13 | 1998 | 17.1565 | 252.65 | −20.34% |
| 81×25 | 7838 | 19.3525 | 439.34 | −10.15% |
| 121×37 | 17518 | 20.0849 | 622.66 | **−6.75%** |
| h→0 三点 Lagrange 外推 | — | **21.5503** | — | **+0.0577%** |

解析参考自身的两条独立校验：面积双求积（手写 Jacobian 复合 Gauss）4.44e-16；p-extension（order vs 2·order）9.90e-16 → 参考值不是同一条代码路径自证。

**D. 逐点 mesh 曲率场 vs 独立算子（Cohen-Steiner/Morvan normal cycle；`tools/crosscheck_curvature_thirdparty.py`）**

| case | faces | H vs 精确场（中位） | K vs 精确场（中位） | H vs normal-cycle（中位 / 最大） | K vs normal-cycle（中位 / 最大） |
|---|---|---|---|---|---|
| sphere | 15584 | 1.58% | 0.05% | 4.35% / 80.1% | 7.50% / 78.6% |
| ellipsoid | 15584 | **12.83%** | **22.43%** | 8.01% / 101.3% | **22.58% / 439.6%** |
| wigley_flank | 2242 | n/a | n/a | **40.00% / 141.2%** | **99.87% / 184.0%** |

→ **逐点值系统性跨过 15% 精度线**（除球面 H/K 的 median），积分量则收敛到 0.06%。被审报告 §5「逐点网格曲率仍不可靠…达不到 15% 判据」与 §6「逐点不满足 <15%，如实报 FAIL」**与实测一致**，无粉饰。

### 1.3 repro_cmd —— 可执行，已在本单复跑（非采信）

```
python examples/fairness_energy/run_example.py           # rc=0, 27.7 s, 11/11 claims green
python tools/crosscheck_curvature_thirdparty.py          # rc=0, 11.0 s, json sha 逐字节相同
python -m pytest tests/test_fairness_energy.py tests/test_curvature.py -q   # rc=0, 43.44 s, 43 passed
```

**但**：以上在**工作树字节副本**上成立，不在任何 git ref 上成立（§2 F1）。

---

## 2. 基准裁判（唯一裁判）

```
SELECT case_name, metric, tolerance FROM tensor_benchmark_cases
 WHERE anchor_slug='GAP-CAD-A3-CURVATURE-FAIRNESS';
--> 0 rows        （全表 53 行；repo='TensorCAD' 仅 GAP-CAD-T1-BASELINE-FLEET 1 行）
GET /api/v1/tensor/bench/l4-verdict/17348
--> {"status":"CLAIMED_NOT_PROVEN","rows":0,"failed_checks":[...,"within_tolerance"]}
    reason: no tensorbench_results row（全表仅 2 行）
```

**裁判缺席（F3）**：本锚点从未登记算例 → L4 的机器口径（平台按 ref/got/error 复算一行）不成立；本次 L4 依「任务榜 goal 文本给出的容差（1e-6 1/m；Wigley 稳定性）+ 解析闭式 + 独立求积」判定，并在证据里显式标注容差来源为**外部给定而非被审单自定**（1e-6 与「对拍 3 例」来自 owner 任务榜原文）。本单**不**代登记 `tensor_benchmark_cases` 行——由被审单自身数字生成裁判属自我循环，留待 owner/协调者按规范登记。

---

## 3. 发现（按严重度排序）

| # | 级别 | 发现 | 证据 |
|---|---|---|---|
| **F1** | **high（交付阻塞）** | 交付物**不在任何 git ref**：`git for-each-ref` 扫全部分支/远端（含 20+ refs）无一处含 `tensorcad/fairness_energy.py`；`git status` 8 项未提交/未跟踪；HEAD 停在他单合入提交 `f43e888`。若平台按 `agent/m17348` 推送，**本能力 0 字节入库**，m17348 的 succeeded 是无仓库载体的 succeeded。 | `audit.json:delivery_state`；sha256 清单已保全 10 文件，可一键恢复 |
| **F2** | medium | **头号验收线的被测对象错位**：`sphere \|k−1/R\|<1e-6` 由 `principal_curvatures(field, …)`（解析隐式场）测得，属既有能力；本单新建的 `principal_curvature_field(mesh)` 只有「conditioned 球 mean H ∈ (0.95,1.20)」这一 ±20% 松 pin（`tests/test_fairness_energy.py:288`），实测 +8%。报告 §5/§6 已披露逐点失败，但 §1 的条款满足度易被读成「新的曲率场达 1e-6」。 | `run_example.py:86 study_sphere_everywhere` 用 `principal_curvatures`；crosscheck 实测表（§1.2 D） |
| **F3** | medium | 裁判缺席：0/53 行；平台 `l4-verdict` 判 CLAIMED_NOT_PROVEN。L4 判定依据为解析闭式 + 独立求积，非登记算例。 | §2 |
| **F4** | low-medium | 「Wigley <3%」是 **h→0 外推值**（+0.058%）；单网格最细 121×37 为 **−6.75%**（>3%，<15%），21×7 更差到 −40%。报告表已列全四档，但「是否 <3%」的结论必须带「外推极限」限定，单次运行值不达 3%。 | `results.json:wigley.rows / richardson` |
| **F5** | low | runner 第 8 条 claim「crease energy grows under refinement」把 `value` 硬编码 `0.0`、`limit=457.5-…差值(465.175)`，**恒真不可伪**；发散方向只由 pytest（`coarse.crease_energy < fine.crease_energy`，`tests/test_fairness_energy.py:265`）保真。故「11/11 claims green」含 1 条空判据。 | `run_example.py` claim 表达式 |
| **F6** | low | 集成面弱：新模块未在 `tensorcad/__init__.py` 导出（`fairness_energy*` 全无），`grep` 除自身 example/tool/test 外 **0 调用者** → 任务榜 L3「接入主求解器/平台」不成立。 | `audit.json:integration` |
| **F7** | info（诚实项） | `meshmagick` 对拍子项 **WITHHELD** 取证齐全：`import` 失败、内网索引无包、沙箱禁出口三路尝试 + 日期写入 `curvature_crosscheck.json:provenance`；替代参考=精确场 + Cohen-Steiner/Morvan normal cycle。替代结果反过来暴露了逐点精度不足（F2），属诚实披露。 | `tools/crosscheck_curvature_thirdparty.py:47-60` |

---

## 4. 判级（L0–L5）

| 级 | 判据 | 结论 | 依据 |
|---|---|---|---|
| L0 | 文档成果 | — | 有真代码，非文档 |
| L1 | 数学/单元验证 | 通过 | 43 测试含闭式 8π/4π、Gauss-Bonnet 4π（rel 1e-9）、双求积一致 1e-12、缩放/重排/平移不变 |
| L2 | 独立功能实现（输入→计算→输出） | **通过** | `tensorcad/fairness_energy.py` 585 行真实现；三种能量 + field/mesh 双路径 + 折角桶分离；CLI 式 runner 产出工件 |
| L3 | 系统集成（接入主求解器/平台） | **不成立** | 无 `__init__` 导出、无跨仓/外部调用者（F6） |
| L4 | 基准验证（完整算例+误差分析+可复现） | **通过（积分核）** | 解析 21.537877 vs 外推 21.5503（+0.058%）；8π 达 +0.053%；四档网格误差表；三例对拍误差表；一键复跑逐叶/逐字节一致 |
| L5 | 科研工程成果 | 不成立 | 无新的科学结论或工程交付物；属诊断能力建设 |

**裁定 = L4**，限定表述：**光顺能量（积分量）达 L4**；**逐顶点曲率场为 L2 级**（有真实现与误差分析，但精度的 15% 线未达，且验收头号条款测的不是它——F2）；**L3 集成面未成立**（F6）；**裁判缺席**（F3）使 L4 只能由解析锚而非登记算例复核。

与 m17333（TensorCAD IGES，board-audit 判 L1：无裁判 + 测试未端到端 + 测试削弱）的差异在于：本单有闭式锚 + 独立算子 + 一键复现 + 双向诚实披露，故不降到 L1。

---

## 5. 未做/未证（not_run，不算 pass）

* **未登记 benchmark 算例**（有意不做，避免自我循环裁判）→ L4 的机器复核口径仍缺。
* **未跨仓验证调用**：未在 TensorPhysics / ship-performance-platform 中实调 `fairness_energy_report`（F6 的直接后果）。
* **未取 meshmagick 对拍**（沙箱/索引均不可得，属环境限制，非被审单遗漏）。
* **未重跑全仓测试套件**（仅 faithfulness 20 + curvature 23 = 43；全仓 262 项不在本单范围，被审单报 262/262 为历史基线）。
* **未修 F1**（本单为核验单；仅在 `artifacts/audit_m17348/delivered/` 做字节保全，恢复命令见 §6）。

---

## 6. F1 恢复路径（供协调者/平台执行）

```bash
# 被审工作树仍在：直接把审计保全的字节写回并提交
cp -a artifacts/audit_m17348/delivered/. <TensorCAD 工作树>/
cd <TensorCAD 工作树> && git add tensorcad/fairness_energy.py tensorcad/wigley.py \
  tests/test_fairness_energy.py examples/fairness_energy tools/crosscheck_curvature_thirdparty.py \
  docs/fairness_energy.md docs/fairness_energy artifacts/curvature_crosscheck.json
git commit -m "m17348: curvature field and fairness energy (GAP-CAD-A3)"
# 校验：sha256 应与 artifacts/audit_m17348/sha256_manifest.json 逐项一致
```

---

## 7. 判据 vs 平台精度标准（>15% 不可接受）

| 关键量 | 误差 | 判定 |
|---|---|---|
| 解析球主曲率 vs 1/R（field 路径） | 2.25e-08 相对 | ✔ 远优于 15% |
| 球 E_bending vs 8π（field 路径） | +0.0530% | ✔ |
| 球 E_bending vs 8π（conditioned mesh） | +5.733% | ✔ |
| Wigley E_bending（h→0 外推） | +0.0577% | ✔ |
| Wigley E_bending（最细网格单次） | −6.75% | ✔（<15%，但 >3%，见 F4） |
| 逐点 H / K（椭球、Wigley 侧面，中位） | 12.8–99.9% | ✘ **超出 15% 判据**（被审单已多次明示「逐点不可用」） |

**显式回答「是否 <3%」**：积分量 ✔（0.053% / 0.058%）；逐点量 ✘（12.8%–99.9%，honest FAIL，不作能力宣称）。

---

## 8. 处置（disposition）

**approve with findings**：

* **数字面**：批准。三件套真实、可复现、误差分析完整、诚实披露充分。
* **交付面**：**F1 为阻塞项**——能力真实但未入任何 ref，需按 §6 恢复并提交后方可「消账」；本单已完成字节保全，不依赖被审工作树存活。
* **宣示面**：能力边界必须带上 F2/F4/F5 的限定语（逐点曲率不可用于 <15% 判据；<3% 是外推极限；11 条 claim 含 1 条恒真）。
* 建议协调者：为 TensorCAD 的 A3 锚点登记 `tensor_benchmark_cases` 行（metric=sphere `|k−1/R|` 与 Wigley E_bending 外推，tolerance 取自任务榜），使后续同锚点不再「裁判缺席」。

---

## 附：图

`artifacts/audit_m17348/audit_figure.png`（三栏）：左＝Wigley smooth 能量 mesh→h=0 收敛与外推；中＝三例逐点 mesh 曲率 vs normal-cycle 的中位/最大误差（15% 线）;右＝smooth/crease 两桶分离（收敛 vs 1/h 发散）。`artifacts/audit_m17348/report_figure.png`＝被审报告 16354 内嵌图（与 `docs/fairness_energy/plot.png` sha256 相同）。