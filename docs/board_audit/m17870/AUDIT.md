# 交差审计 · m17833 GAP-CAD-D2-GOLDEN-REGRESSION · TensorCAD 几何金标回归套件

审计单 m17870 · 对象 m17833（succeeded）· 被审提交 `f2222d1e516f`（base `9f12643dbc8c`，尚未并入 origin/main）·
判据 = `tensor_benchmark_cases` 4 行（id 65, 66, 67, 68）+ audit 侧独立复跑。

**结论：L4（基准验证）成立——`confirmed` 前提是打掉下面唯一一个 blocking 缺陷（repro_cmd 少一级 `..`，实测 rc=2）。**
本次审计不是在读报告，是把交付件在自己的树里重跑：`--check` rc=0、benchmark 30 行**逐位一致**、47 个测试全绿，
另外把 m17833 self-claim 的五个数字逐个对表（见 §⑤）。

## ① 案例介绍（物理背景 / 几何 / 参数）

被审对象不是一条船，而是 **4 条几何金标的回归闸门**：把 TensorCAD 的几何/静水力生产管线
（`wigley.generate_wigley`、`series60.generate_series60`、SDF 球/椭球 → marching tetrahedra、`dtmb5415.DTMB5415Space`）
钉在"参考值可被重新计算的"锚点上，几何一改就响。

| 案例 | 参考源（`reference_source`，运行时重算，非表内字面量） | 几何 sha256（前置 16） | 耗时 s | 断言量 | 判 |
|---|---|---|---|---|---|
| wigley | tensorcad.wigley.wigley_analytic_hydrostatics (closed forms) | 0ef07ca68299b621 | 1.8 | 8 | PASS |
| series60 | tensorcad.series60.analytic_* over the exponents fitted from Cb | f9da70a467f6d427 | 4.2 | 8 | PASS |
| sphere_ellipsoid | (4/3)πr³, 4πr² and (4/3)πabc over the exact primitive fields | a0f68b295bec37a3 | 16.1 | 9 | PASS |
| dtmb5415 | docs/data/dtmb5415/holdout_fixture.npz:per_config_rms_m[0] | 64fc9d7694be20f5 | 2.0 | 5 | PASS |


粗/细两级网格（e.g. Wigley/Series-60 coarse 81×25 → fine 161×49），容差 **由两级之差推导**：`tol = 2·|fine−coarse|/|ref|`，下限 1e-4；
另有一条 `1e-9` 基线漂移锁（`docs/golden/tensorcad_geometry_golden.json`，sha256 记在 `run_record.json`）。

## ② 怎么做的（审计方法 / 工具链）

1. **取证**：读 `/root/tp-agent/work/m17833` 的 result/changed_files，逐项落地检查 impl/entry/artifacts/tests/callers/benchmark/repro_cmd。
2. **独立复跑（不采信自报）**：
   - `PYTHONPATH=. python3 -m tensorcad.golden --check --json ...` → **rc=0, 35.21 s**, 4/4 case PASS, 0 error finding；
   - `PYTHONPATH=. python3 benchmarks/verified/golden_geometry/run_benchmark.py --out /tmp/audit17870/bench` → rc=0，**30 行 / 9 门 PASS**，24.1 s；
   - 与交付工件 `artifacts/golden_benchmark/result.json` 逐指标对比：**最大绝对差 0.0e+00（即逐位相同）**，几何指纹 3351b57966fc4a44 一致；
   - `pytest tests/test_golden_regression.py tests/test_golden_benchmark.py -q` → **47 passed, 217.27 s**（审计侧重跑，非引用其自报）。
3. **基准裁判对表**：用 `tensor_benchmark_cases`（anchor_slug='GAP-CAD-D2-GOLDEN-REGRESSION'）的 4 行 metric 串**机器解析**
   （`docs/board_audit/m17870/collect_evidence.py`），把每行的 `got` 与交付工件 `measured` 逐项比对：缺失 0 项，最大差 2.5e-09（DB 串只保留 10 位有效数字）。
4. **攻击面**：把交付件"可复现"的宣称当假设来打——**真去执行 README 承诺的 `repro_cmd.sh`**，于是抓到了唯一 blocking 缺陷。
5. **工具链**：python3.11.6 / loongarch64（与交付同构），matplotlib 出图，psql 只读取裁判行。

## ③ 计算结果（ref vs got；容差逐指标，非统一放宽）

裁判行（4 例 / 24 个可断言指标，全部 `status=passed`）：

| row | case | 指标数 | 该行最大 |rel err| | 该行最坏容差 | 余量 | 判 |
|---|---|---|---|---|---|---|
| 65 | CASE-CAD-D2-GOLDEN-WIGLEY | 8 | 1.953e-04 | 4.329e-03 | 22.17x | passed |
| 66 | CASE-CAD-D2-GOLDEN-SERIES60 | 6 | 8.591e-03 | 1.831e-02 | 2.13x | passed |
| 67 | CASE-CAD-D2-GOLDEN-SPHERE_ELLIPSOID | 9 | 8.860e-04 | 5.457e-03 | 6.16x | passed |
| 68 | CASE-CAD-D2-GOLDEN-DTMB5415 | 1 | 2.645e-08 | 1.000e-04 | 3780.72x | passed |

最差 8 行（全部 30 行见 `artifacts/golden_benchmark/benchmark_table.md` 与审计侧 `audit_reproduction.json`）：

| case | metric | reference | measured | rel err | 容差 | 判 |
|---|---|---|---|---|---|---|
| series60 | block_coefficient | 0.6 | 0.5948455588 | -8.591e-03 | 1.831e-02 | PASS |
| series60 | displacement_volume_m3 | 5.4 | 5.353610029 | -8.591e-03 | 1.831e-02 | PASS |
| series60 | maximum_section_area_m2 | 0.864 | 0.8566071822 | -8.557e-03 | 1.811e-02 | PASS |
| sphere_ellipsoid | sphere_volume_m3 | 4.188790205 | 4.185079128 | -8.860e-04 | 5.457e-03 | PASS |
| sphere_ellipsoid | ellipsoid_volume_m3 | 6.283185307 | 6.277618692 | -8.860e-04 | 5.457e-03 | PASS |
| sphere_ellipsoid | sphere_surface_area_m2 | 12.56637061 | 12.56063757 | -4.562e-04 | 2.816e-03 | PASS |
| wigley | block_coefficient | 0.4444444444 | 0.4443576416 | -1.953e-04 | 4.329e-03 | PASS |
| wigley | displacement_volume_m3 | 2.777777778 | 2.77723526 | -1.953e-04 | 4.329e-03 | PASS |


交付件自报的 9 道门（审计侧复跑**逐门复现**）：

| gate | verdict | detail |
|---|---|---|
| G1_golden_library_complete | PASS | cases=['dtmb5415', 'series60', 'sphere_ellipsoid', 'wigley'] |
| G2_wigley_reference_match | PASS | metrics=8 failing=[] errors=[] |
| G2_series60_reference_match | PASS | metrics=8 failing=[] errors=[] |
| G2_sphere_ellipsoid_reference_match | PASS | metrics=9 failing=[] errors=[] |
| G2_dtmb5415_reference_match | PASS | metrics=5 failing=[] errors=[] |
| G3_stored_references_reproduce | PASS | worst 0.000e+00 on wigley/displacement_volume_m3 over 24 references |
| G4_baseline_drift_lock | PASS | drift findings=0 |
| G5_platform_error_limit | PASS | worst series60/block_coefficient = -8.591e-03 |
| G6_geometry_fingerprint | PASS | wigley=0ef07ca68299; series60=f9da70a467f6; sphere_ellipsoid=a0f68b295bec; dtmb5415=64fc9d7694be |

## ④ 结果图（图 `audit_golden.png`）

左图 = 本次审计重跑的 30 行 `|measured−reference|/|reference|`（蓝=解析锚点，灰=恒等式），红短线=逐指标推导容差，
红点线=3% 精度线，黑虚线=平台 15% 拒收线；右图 = 裁判表 4 行的"该行最大误差 vs 该行最坏容差"。
图由 `collect_evidence.py` 直接从 result.json 与 DB 行生成（零手抄）。

## ⑤ 注意点 / 发现的缺陷（严重度排序）

1. **[blocking, 已给补丁] `benchmarks/verified/golden_geometry/repro_cmd.sh` 交付即不可执行**：`repo=$here/../..` 少一级 →
   解析成 `<repo>/benchmarks`，脚本立刻 `python3: can't open file '<repo>/benchmarks/benchmarks/...'`，**实测 rc=2**。
   README 与 result.json 都把"Reproduce with repro_cmd.sh"写成复现入口，所以这是包装层缺陷，不是数值缺陷：
   底层 `run_benchmark.py` 直接调用 rc=0 且数值逐位一致。修复 = 一行（`../..`→`../../..`），
   补丁 `docs/board_audit/m17870/repro_cmd_fix.patch`，已在交付树副本上 `git apply --check` 通过并实跑
   **rc=0 / 35.8 s / 30 行 / 与交付工件最大差 0.0e+00**。
2. **[medium] 4 个金标案例里只有 2 个有真外部锚点**：Wigley 闭式解、球/椭球 `(4/3)πr³ / 4πr²` 是解析真值；
   **Series-60 的参考量是同一模型族的 `analytic_*`**（8.59e-3 是离散自洽，不是文献表对拍）；
   **DTMB5415 的 holdout_rms 参考值就是该模型自己的预测**（吻合 2.6e-8，等于同义反复），该类只有 `Cb=V/(LBT)` 与镜像不变式有信息量。
   这不否定"回归闸门"的用途（它本来就是抓漂移），但 README 若把 4 例都表述为"参考对拍"就说过头了。
3. **[low] 工件 provenance 落后一个提交**：`result.json` 记 `git_commit 3c267f79`（中间提交），
   `tensor_benchmark_cases` 4 行记 `commit=9f12643`（base），最终提交是 `f2222d1`；数值经审计复跑证明未受影响，但溯源链指错提交。
4. **[low] 台账未收口**：`feature_ledger` 的 `GAP-CAD-D2-GOLDEN-REGRESSION` 仍是 `status='open'`；m17833 结果里 `ledger_gate=unanchored`。
5. **[unknown, 不记 pass] GitHub CI 未验证**：`.github/workflows/golden-geometry.yml` 走 `runs-on: ci-loong64` 自托管 runner，
   审计无 runner 权限，**只能确认文件接线正确（含 schema/断言步骤），不能确认它跑过**——按 unknown 记，不计入 L4 证据强度。
6. **[positive] 负面对照是真的**：`test_hand_edited_reference_is_rejected`、`test_missing_baseline_fails_closed`、
   `test_quick_mode_still_catches_a_drift`（1e-8 漂移在 pre-commit 快模式仍 FAIL）、`test_benchmark_rejects_a_perturbed_reference`——
   闸门不是装饰，47 项测试是断言不变式而不只是执行代码。

## ⑥ 误差判定（关键量 vs 基准，是否 <3%）

| 判据 | 阈值 | 实测 | 是否 <3% |
|---|---|---|---|
| 24 个解析锚点指标的最大相对误差（series60/block_coefficient） | 3% (平台精度线) | 0.859% | **是** |
| 30 行中超过 3% 的行数 | 0 | 0 | **是** |
| 30 行中超过 15%（平台拒收线）的行数 | 0 | 0 | **是** |
| 每行 |rel err| ≤ 该行推导容差 | 全部满足 | 0 违反（30/30） | **是** |
| G5 平台线（最差断言误差 ≤ 0.15） | 0.15 | 8.591e-03 | **是** |
| 复现一致性 | 逐位相同 | 最大差 0.0e+00 | **是** |
| GitHub CI 实跑 | 跑过 | 未验证 | **unknown（不计 pass）** |

诚实边界：最外一层真外部锚点只覆盖 Wigley + 球/椭球（24 个可断言指标里的 17 个），Series-60/5415 的"通过"
是模型内自洽/漂移锁；因此本判定给的是 **L4 = 工程基准回归验证**，不是"几何算法对外部实验数据的 3% 验证"。判级按功能实质，不按报告页数。
