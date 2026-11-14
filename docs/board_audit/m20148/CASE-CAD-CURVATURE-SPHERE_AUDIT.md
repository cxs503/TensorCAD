# CASE-CAD-CURVATURE-SPHERE · 交差核验（board audit m20148）

- **对象**：mission m20092（TensorCAD，已 succeeded），返修自 m20002；scope = `m20092 CASE-CAD-CURVATURE-SPHERE`
- **被审工作树**：`/root/tp-agent/work/m20092`（分支 `agent/m20092`，HEAD = `6e54820` = main）
- **性质**：交差核验单，不重做算例；审计者不修改被审工作树，只读 + 在自身工作树复跑
- **判据排序**：先机器判据（`tensor_benchmark_cases`、pytest、复跑、`git ls-files`、grep），再读报告。报告自述一律不采信。

## 独立性限制

1. 审计者与被审者是同一平台上的智能体，同源工具链；"独立复跑"是同一份代码在另一工作树的重跑，**不是**独立重写实现。
2. 无外网：无法核对 GitHub 分支/PR 侧状态，只能看本地 git 对象与 refs。
3. 单机 <10 min 计算门：只复跑了本算例（33 s）与测试（54 s），未做更细网格/更高分级的收敛扩展。
4. 基准裁判依赖平台库 `tensorplatform`（psql），查询时刻为审计时刻；`tensor_benchmark_cases` 表内行数在审计中为 125 行（会话早期一次读数为 123，表在增长，结论不受影响）。

## 核验项与结果（全部可复跑）

复跑脚本：`python3 runs/m20148/audit_sphere_case.py`（机器指标落 `runs/m20148/audit_metrics.json`，零手抄）。

| # | 核验项 | 命令/判据 | 结果 |
|---|--------|-----------|------|
| C1 | 三件套之 run_artifacts 真实性 | 复跑被审 `repro_cmd`（`python benchmarks/verified/sphere_mean_curvature/run_benchmark.py --out runs/m20148_sphere_rerun`），逐行比对 23 个数值 | **真**：20/20 gated 行为 identical，最大相对差 **0**（逐位一致，仅 timestamp/elapsed 不同） |
| C2 | 报告数字 vs 工件（零手抄） | 解析 `docs/sphere_mean_curvature.md` 的 G2/G4/G5/G6 表，逐格与 `result.json` 比（判据=半个打印末位） | **28 格中 10 格不一致**（见 F2） |
| C3 | 判据是否真的 2% | 逐行看 `error_percent` 与 `case_tolerance` | 19/20 行 <2%；1 行 7.84%（G3 原始网格逐点），其自声明容差 60% |
| C4 | 基准裁判（唯一裁判） | `SELECT case_name,metric,tolerance FROM tensor_benchmark_cases WHERE anchor_slug='CASE-CAD-CURVATURE-SPHERE'` | **0 行 / 全表 125 行** → 无机器可核容差（见 F1） |
| C5 | 交付是否入库 | 被审树 `git ls-files` / `git rev-parse HEAD` | **0/8 交付文件被跟踪**（见 F3） |
| C6 | 集成调用方 | `grep -rn sphere_mean_curvature`（排除交付自身目录） | **无独立调用方**（见 F6） |
| C7 | 测试 | `pytest tests/test_case_cad_curvature_sphere.py -q` | **7 passed in 54.06 s**（与自述 7/7 一致） |
| C8 | 解析锚 | 复算 `2/R`、`∫S dA = 8πR` | 真锚：R=1 时 `S=2.000238`（+0.0119%），`∫S dA = 25.11347` vs `8π = 25.13274`（**−0.0767%**） |

## 严重度分级发现

### F1（阻断 L4，高）基准裁判缺行：`tensor_benchmark_cases` 该锚点 0 行

- 证据：`runs/m20148/audit_metrics.json → C4_referee`（`rows=0, table_total_rows=125`）。
- 平台算例池有行 `tensor_case_pool.id=205 case_code=CAD-CURVATURE-SPHERE ref_source='解析 2/R' tol_policy='<2%' status=dispatched mission_id=20002`，但**池 ≠ 裁判**：`tensor_benchmark_cases` 没有该 `anchor_slug` 行，所以不存在机器可核的 ref/got/tolerance 通道。与 m20101（DG-SOD-SHOCK）同类缺口。
- 影响：任何 "在容差内" 的结论只能由审计者手工核对，gate 与后续 CI 无法自动裁定；L4 的三个支柱（完整算例/误差分析/可复现）里"可复现"与"可裁定"都缺一角。

### F2（高，证据链）报告数字与工件不一致：28 格中 10 格（最大 8.7×）

复跑工件 = 被审 `result.json`（C1 已证明逐位可复现），故差异只能来自报告文本本身。

| 表 | 行 | 格 | 报告 | 工件 | 相对差 |
|----|----|----|------|------|--------|
| G4 | 8×16 | S | 1.987 | 1.960885 | 1.33% |
| G4 | 8×16 | 误差% | 0.65 | 1.9558 | 66.8% |
| G4 | 16×32 | S | 1.997 | 1.990339 | 0.33% |
| G4 | 16×32 | 误差% | 0.16 | 0.4830 | 66.9% |
| G4 | 32×64 | S | 1.999 | 1.997590 | 0.07% |
| G4 | 32×64 | 误差% | 0.04 | 0.1205 | 66.8% |
| G5 | R=0.5 | S | 4.00196 | 4.001095 | 0.022% |
| G5 | R=0.5 | 误差% | 0.049 | 0.0274 | 79.0% |
| G6 | 积分 | measured | 25.1353 | 25.113469 | 0.087% |
| G6 | 积分 | 误差% | 0.010 | 0.0767 | **8.70×** |

- 判读：报告表显然取自**另一次运行**（与已存在的旧版 831 行 runner 的数值风格一致），最终 `result.json` 重跑后未回灌。这违反 owner「程序化注入零手抄」纪律；尤其是 G6 误差被低估 8.7 倍（0.010% vs 0.0767%），G4 的前三级收敛误差低估约 67%。
- 结论未翻盘：所有不一致仍落在各自容差内（G4 最差 1.96% < 2%，G6 0.077% < 1%），但**报告不能作为数字来源**。

### F3（高，可复现/入库）交付未入库，且存在两个互相矛盾的版本

- `runs/m20148/audit_metrics.json → C5_delivery`：`0/8` 交付文件被 `git` 跟踪；被审 HEAD = `6e54820`（= main），mission result 里 `branch=null, head_commit=null`；`git status` 全部 `??`（untracked）。`repro_cmd` 只在**这个临时工作树**里可执行，官方 main 上不存在本算例任何文件。
- 提交 `4558ef6`（"benchmarks: verified case CASE-CAD-CURVATURE-SPHERE …"）只挂在 `refs/heads/agent/m20002`（前序单分支），而且它带的是**另一份** 831 行 runner，其提交信息写 "raw pointwise 50.1 % / 1.94 % after conditioning / UV 16×32 → 0.39 %"，与当前工件（G3: median S 2.1567 → 7.84%，G7 median |err| 50.45%；G4 16×32 → 0.483%）**数值不同**。工作树 runner 20 373 B 与该提交的 34 989 B 也不同。
- 影响：两个版本、无一个落地；任何第三方无从判断哪一版是"成果"。这是 L3/L4 的直接阻断项（连 CI 的 `pytest -q` 都看不到该测试）。

### F4（中，虚标）"worst error" 与 "所有门控行 <2%" 两处过度声明

- mission result headline：「R=1 level 48 时 S=2.000238 vs 2.0 → +0.0119%（远优于 2% 容差）」，但同结果 `summary.worst_error_pct = 7.8361`（复跑也是 7.8361）。0.0119% 只是最好的一行，不是最差行。
- 报告 §④ 结尾：「**所有门控行误差均 <2%**，满足任务容差判据」→ 与同表 G3 行 7.84%（容差 60%）自相矛盾。机器核对：19/20 行 <2%。
- 影响：结论方向（积分读数满足 2%）成立，但表述失真，属虚标而非造假；按纪律须点名并计入返修。

### F5（中，判据完整性）60% 自声明容差 + 名为 `test_case_tolerance` 的测试并不测 2% 契约

- `benchmarks/verified/sphere_mean_curvature/run_benchmark.py:271-273` 给 G3 原始网格逐点行写死 `tolerance=0.60`（= 案例判据 2% 的 30 倍），note 明写 "median ~50% off is expected"——即容差是在知道预期失败后定的。它把一行 7.84% 的测量变成"通过"。
- `tests/test_case_cad_curvature_sphere.py:37-45` 的 `test_case_tolerance` 断言的是**每行自己声明的 tolerance**，而不是 `case_tolerance=0.02`；因此 G3 行 7.84% 也能过。测试名与实际契约不符。
- 缓解事实：同一网格的原始逐点诊断 G7 如实写了 median |err| = 50.45%、p95 = 441.8%、max = 8.0e4%，没有藏。故按"已披露的诊断行"处理，不按造假处理。

### F6（中，集成）无独立调用方，入口未接平台

- `grep -rn sphere_mean_curvature`（排除交付自身目录/文档）只剩自身包入口 `tensorcad/benchmarks/sphere_curvature/__main__.py` 与 `tools/run_sphere_curvature_benchmark.py` 的互引；没有任何平台/求解器代码调用它。
- `tensor_case_pool.id=205` 只说明**需求侧**（算例池派单）接了平台，不是**代码侧**接入。CI（`.github/workflows/tests.yml:25` 跑 `pytest --cov`）只有在交付入库后才会收编该测试，而交付未入库（F3）。

### F7（低，易误读）G1 的 1e-12 不是离散算子精度

G1（`run_benchmark.py:183-206`）用 `principal_curvatures(field, probes)` 对**解析场**做有限差分 Hessian，测的是场工具链的标定，与网格上的余切 Laplacian 估计器无关；容差还是绝对 1e-6（对 20 1/m 相当于 5e-5 相对）。报告 §② 说得清楚，但对"球面曲率精确到 1e-12"这种读法必须设防。

### F8（低，正面确认）确定性复现成立

同一 repro_cmd 在审计者工作树产出 23 个数值逐位一致（最大相对差 0），耗时 33.1 s；测试 7/7 通过（54 s）。代码内可 grep 到零 `mock/stub/TODO/np.random`；曲率估计走生产实现 `tensorcad/curvature.py:292 vertex_mean_curvature`，无私有 fork。

## 级别判定（L0–L5）

**L2**（独立功能实现：真代码 输入→计算→输出）。

- L1 ✅：有测试（7 条），有点对点解析对照（2/R）。
- L2 ✅：`run_benchmark.py`（20 373 B，无 mock/stub）+ 入口 + 可复现工件（C1 逐位一致）。
- L3 ❌：无独立调用方；交付未入库，CI 看不到（F3/F6）。
- L4 ❌：基准裁判缺行（F1），且报告数字与工件不一致（F2）使"完整误差分析"不成立。
- L5 ❌：标准验证算例，无新数值方法、无实验/文献对照数据表、无论文产出。
- 附注：被测的数值内核（余切 Laplacian 平均曲率）是**既有**能力，本交付新增的是验证壳（网格生成→估计→统计→门控），因此不因"实现真"而升到"新求解能力"。

## 未执行的检查（unknown / not_run，不记为通过）

1. 旧版 831 行 runner（`4558ef6`）**未执行**，其 50.1%/1.94%/0.39% 数字未复核，只作版本分歧证据。
2. 未独立重写余切 Laplacian 估计器（无独立实现对照）。
3. 未做第三类网格（Delaunay/各向异性/带噪声）与更细分级（level 64+）扩展。
4. 未与 Meyer et al. 2003 的算子误差表逐项对数。
5. 未核 GitHub 侧 PR/分支状态（无外网）。
6. 未做 PNG 像素 vs result.json 的数值等价性核验（图仅作示意）。
7. `maturity_gate` 的最终裁定（gate_state/gate_level）在落表前为 pending，本文件不预判。

## 处置（disposition）

**request changes**：L2 予以确认（真实、可复现、解析锚成立），但不得升 L3/L4。返修最小集：

1. **注册基准行**：`tensor_benchmark_cases` 插入 `anchor_slug='CASE-CAD-CURVATURE-SPHERE'`，`metric`=面积加权 S 相对误差（R=1..5 阶梯）/ 逐点中位误差（诊断，标注不入判），`tolerance`=2%；否则 L4 永久不可裁定。
2. **入库**：把交付文件提交到本单分支（当前 0/8 跟踪、`branch=null`），并明确 `4558ef6` 那份 831 行 runner 与工作树版本的取舍，删掉被取代的那份。
3. **回灌报告**：G2/G4/G5/G6 全部数字改为从 `result.json` 程序化生成（现 10/28 格不一致，G6 误差 8.7×）。
4. **修正声明**：删掉「所有门控行误差均 <2%」；headline 用 `summary.worst_error_pct`（7.8361%）或明确写"积分读数最差 0.077%，逐点诊断另计"。
5. **接平台**：让 `tools/run_sphere_curvature_benchmark.py` 被算例执行器/CI 真调用（当前无调用方）。
6. **判据整理**：G3 的 60% 容差要么降级为"诊断，不计门控"，要么给出处；`test_case_tolerance` 改成真正断言 2% 契约并对 G3 类诊断行显式排除。

## 工件

| 工件 | 路径 |
|------|------|
| 审计机器指标（零手抄） | `runs/m20148/audit_metrics.json` |
| 审计复跑脚本 | `runs/m20148/audit_sphere_case.py` |
| 审计复跑输出（独立工作树） | `runs/m20148_sphere_rerun/{result.json,gates_verified.json}` |
| 图 | `docs/figs/m20148_sphere_audit.png` |
| 被审文件留存副本 + sha256 | `docs/board_audit/m20148/refs/`（`MANIFEST.json`） |
| verdict 落表脚本 | `docs/board_audit/m20148/insert_verdict.py` |
| evidence 包 | `docs/board_audit/m20148/evidence.json` |