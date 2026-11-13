# 交差审计 m17397 — m17317 `GAP-CAD-T2-PARAM-SERIES`（TensorCAD Series 60 家族返修轮1）

审计人：board-audit（m17397）· 对象：mission **m17317**（m17252 返修轮1，锚点 GAP-CAD-T2-PARAM-SERIES）
· 审计工作树：`/root/tp-agent/work/m17397`（main，含全部已并 mission）
· 被审提交：`db0238e`（author `agent/m17317`，2026-11-13 04:17）

## 0. 结论

**级别 L2**（独立功能实现：真代码 输入→计算→输出；含真实算例、误差分析、可复现工件）。
**L3 未达成**：m17317 的提交 `db0238e` **不是 HEAD 的祖先**（`git merge-base --is-ancestor db0238e HEAD` → NO），
`pr_gate=null`（无 PR 合并记录）——CLI 入口 / `__init__` 导出 / delivery_page / run_artifacts **只存在于未合并提交**，
main 分支里没有这些文件（`tests/test_series60_cli.py`、`tools/series60_family/delivery_page.py`、`docs/m17317/` 均不存在）。
**L4 未达成**：`tensor_benchmark_cases` 对 `GAP-CAD-T2-PARAM-SERIES` 锚点 **0 行**（无注册判据/容差）；
3% 容差是 `tools/series60_family/audit.py:63` 的**代码内预置**，非台账预注册。

**关键缺陷（material）**：m17317 的「succeeded」是 `TP_AUTO_APPROVE_DEV_ACCEPTANCE` 自动批准，
**未核验 PR 合并**（`pr_gate=null`）。全部交付物只在未合并提交里，不在 main——任何下游消费者今天都调不到。
且 `fail_stage=budget_exhausted+forced_finish`（预算耗尽被迫收尾），部分检查未跑
（「未跑：HUL/ANA 旧 10 变异体、集群工件、与试验值对照」）。

## 1. 取证：三件套逐项核验

| 项 | 声明 | 实测 | 裁定 |
|---|---|---|---|
| run_artifacts | `docs/m17317/artifacts/*` | 在提交 `db0238e` 里：five_cb_gradient.json / transom_sweep.json / bulb_sweep.json / gate_power_series60.json / test_run.log / 3 张 png，均为真实数值（非占位） | ✅ 真（**但只在未合并提交**） |
| ref vs got | 五档 Cb × 五系数 vs `docs/data/series60_published_parents.json` | 逐行核对：全部 `within=True`（3% 门），最差 +0.34%（Cm @Cb=0.8）；**但性质是闭合**（Cp/Cm 硬编码在 `_PARENT_FORMS`，+0.38% 测的是公开表自身 Cb=Cp·Cm 互不自洽） | ✅ 数值真，**性质=闭合**（F1） |
| repro_cmd | `python3 tools/series60_family/audit.py` | 在提交 `db0238e` 里存在，写 `docs/m17317/artifacts/series60_family/*`；**不在 main** | ✅ 可执行（未合并） |

基准裁判：`SELECT ... FROM tensor_benchmark_cases WHERE anchor_slug='GAP-CAD-T2-PARAM-SERIES'` → **0 行**（全表 53 行均无此 slug）。
报告六段：tensor_reports id 16314，`## 1`–`## 6` 齐全（案例/方法/数字表/图/注意点/误差判定）。

## 2. 我自己跑出来的数（零手抄，程序化注入）

* `git show db0238e:docs/m17317/artifacts/series60_family/five_cb_gradient.json` → 五档 Cb 全部 `within=True`：
  Cb=0.6 块系数 -0.000355、Cb=0.7 -0.000355、Cb=0.8 -0.00044；midship 最差 +0.003449（Cb=0.8）；网格-闭式位移最差 -0.00044。
* `git show db0238e:docs/m17317/artifacts/gate_power_series60.json` → 14 个突变体全部 `status:'ok'`（全杀）。
* `git diff HEAD db0238e -- tensorcad/series60.py` → 余弦 keel 聚类（`_level_zeta`）+ VCB 符号修正，真代码。
* `git show db0238e:tests/test_series60_cli.py` → CLI 端到端（`--family series60 --gate --stl` → rc=0，gate passed + stl watertight）。
* `git merge-base --is-ancestor db0238e HEAD` → **NO**（未合并）。

## 3. 发现（按严重度）

### F1【material·声明问题】「succeeded」未核验 PR 合并
`pr_gate=null`，提交 `db0238e` 不是 HEAD 祖先。m17317 的「succeeded」是自动批准（`TP_AUTO_APPROVE_DEV_ACCEPTANCE`），
未核验合并。全部交付物（CLI 入口/`__init__` 导出/delivery_page/run_artifacts）只在未合并提交里，**main 里没有**。
任何按「已 succeeded」读 main 的门都会读到「能力不存在」。这是本单最重的发现——**不是造假**（提交与工件真实），
是**验收流程未核验合并**。

### F2【material·级别封顶】L3/L4 未达成
L3（系统集成）要求有 main 里的下游消费者——无（未合并，`git ls-tree HEAD` 无这些文件）。
L4（基准验证）要求 `tensor_benchmark_cases` 有注册判据——0 行。3% 容差是代码内预置，非台账。

### F3【中·声明问题】ref-vs-got 是「闭合」不是「独立」
`_PARENT_FORMS` 硬编码公开 Cp 阶梯（0.625…0.805），`Cm=Cb/Cp`、`Cp=Cb/Cm` 都是该阶梯的算术结果——
没有任何一个系数独立于公开表。+0.38%「对公开值」测的是公开表自身 Cb=Cp·Cm 互不自洽，不是实现精度。
（与 m17308 审计 F1 同源；m17317 报告文字诚实，但机器可读字段把 `prismatic/midship` 标为 `shape-derived (independent)`。）

### F4【低·完整性】forced_finish
`fail_stage=budget_exhausted+forced_finish`；「未跑：HUL/ANA 旧 10 变异体、集群工件、与试验值对照」——
部分验收项未执行。

## 4. 六段报告核验（tensor_reports id 16314，8933 字符）

① 案例（含几何/参数表/provenance）✅ ② 方法 ✅ ③ 数字表 ✅（与工件一致）④ 图 3 张 ✅
⑤ 注意点（球鼻是独立实体、τ 是宽度比、网格首阶亏损、不验证阻力/稳性）✅ ⑥ 误差判定（含「<3%」逐项回答）✅
**但**⑥ 表首行把「闭合」当「published comparison」，与 F3 同源——文字层面已注明，不算虚标。

## 5. 未做 / 不能做的

* 未复核公开母型数值的文献原始出处（沙箱禁外网，无 Todd 原表可对）——F4 的 WITHHELD 即此。
* 未验证带球鼻/方尾船型的水动力（阻力/稳性）——本单也未声称。
* 未跑全仓回归（m17317 声称 157 passed，我只核验了提交内工件与 CLI 测试代码）。
* 未检查 ship-performance-platform 侧是否消费该几何（跨仓，超出本单）。

## 6. 处置

`approve with notes`（非阻塞）：**判 L2**（功能实质=真代码+真工件+ref-vs-got 在容差内），
**不进 L3**（未合并入 main，无下游消费者）**也不进 L4**（无注册基准判据）。
**关键返修**：m17317 的「succeeded」必须补 PR 合并（把 `db0238e` 合入 main），否则能力只在未合并提交里，
任何按 main 读的门都判「不存在」。F3 的 `role` 字段应改为 `input ladder (closure check)`。

Verdict 行：`('TensorCAD', 'm17317 GAP-CAD-T2-PARAM-SERIES', 17397, 2, evidence, 'board-audit')`（id 568）。
