# 交差核验：m18646 GAP-CAD-C3-LINES-FAIRING（型值表三向光顺）— 判 **L2**

- 核验单：m18778（board-audit，交差核验，非重做）
- 对象：mission **18646**（status=succeeded，repo=TensorCAD，finished 2026-11-13 15:58Z）
- 交付 commit：`4fa7b2a`「TensorCAD: three-directional offset-table fairing (GAP-CAD-C3-LINES-FAIRING)」，7 files / +2204
- 工作树：`/root/tp-agent/work/m18646`（分支 `agent/m18646`）
- 独立限制：核验在同一台 loongarch64 主机、同一 Python 3.11 上重跑；未做跨架构/跨版本复跑；
  未取得任何人工手光顺型值表；`tensor_benchmark_cases` 有该表但在本锚 0 行；核验者未修改被核对象任何文件。

## 1 取证：三件套逐项验收

| 三件套 | 声称/实测 | 判定 |
|---|---|---|
| 实现文件 | `tensorcad/lines_fairing.py` 1155 行：离散密切圆曲率 `signed_curvature`、翻转计数 `curvature_sign_flips`（2 % floor）、曲率变分能量 `line_bending_energy`、带状 Cholesky `_cholesky_banded_solve`、三阶差分算子 + λ 几何二分的经典光顺样条 `_energy_kernel`、MLS 二次 + 形状屏蔽的 `local` 核、等 y 值反解 `_section_crossing`/`_pass_buttocks`、驱动 `fair_offset_table` | **真代码**（无 mock/随机分支） |
| 入口 | `fair_offset_table(table, FairingSettings) -> FairingResult`；`PYTHONPATH=. python3 examples/lines_fairing/run_example.py` | 实测 exit 0，7.0 s |
| run_artifacts | `docs/lines_fairing/results.json`（461 行）、`docs/lines_fairing/report.md`、`docs/lines_fairing/plot.png`（99 441 B）、`docs/lines_fairing.md` | 4 件全部存在 |
| tests | `tests/test_lines_fairing.py` 34 项，实测 `34 passed in 6.60s`（独立重跑，非引用） | **PASS** |
| repro_cmd | `PYTHONPATH=. python3 examples/lines_fairing/run_example.py` → 重生成 `results.json` 与提交版 **逐字节等价（json 相等，True）**；打印 `{"acceptance":"PASS"}` | **可复现** |
| 基准裁判 | `SELECT case_name,metric,tolerance FROM tensor_benchmark_cases WHERE anchor_slug='GAP-CAD-C3-LINES-FAIRING'` → **0 行**（全表 110 行，TensorCAD 仅 6 行，全为 GAP-CAD-D2-GOLDEN-*） | **无判据行** |

## 2 基准判据（feature_ledger 验收）+ 实测

验收原文：`光顺前后型值偏差<给定容差 且 曲率符号翻转数单调下降（≥2 船数据）；与手工光顺参考对比`

| 条款 | 实测（results.json，独立复跑一致） | 判定 |
|---|---|---|
| ① 偏差 < 容差 | wigley 0.00892 ≤ tol 0.01；series60 0.01000 ≤ 0.01（浮点 1e-12 松弛）；kvlcc2 0.05000 ≤ 0.05；`fair_offset_table` 末段超带即 `raise FairingError` | **达标**（注意：三条全部贴带饱和，属"约束起作用"而非富余） |
| ② 翻转数单调降，≥2 船 | total flips：wigley 219→81→59；series60 450→224→182→166→146→139；kvlcc2 272→90→47→33→29 | **达标**（3 船） |
| ② 附注 | 单调性是**驱动层构造保证**（pass 只在不劣化 `(flips, energy)` 时提交，退化则半松弛重试/拒收）。分方向：kvlcc2 水线族 0→0→15→19→21（**变差**），站线族 272→74→…→7 补偿；文档显式披露，属诚实但非"每方向都变好" | 披露成立 |
| ③ 与手工光顺参考对比 | 无人工光顺型值表。替代：`hand_fairing_reference`（独立 MLS 二次平滑）+ 解析半宽面。wigley：与本驱动 max 0.0051 m / rms 0.0011 m（核验者独立复算，与 docstring 一致）；kvlcc2 上该参考无带约束，自身发散（rms 0.36 m、max 2.61 m），对比无意义 | **部分达标（代理参考）** |

精度（相对平台 >15 % 不可接受线）：以解析半宽面为参照，wigley 21×17 RMS 0.00216→0.00137 m（local）/0.00050（energy），相对半宽 5 m 为 **0.03 %** 量级；kvlcc2 33×25 RMS 0.01122→0.01060（−5.6 %），**max 误差 0.0200→0.0559 m 变差**（带约束的代价，文档已披露）。

## 3 判级（L0–L5）

- **L2 成立**：真实现（输入→计算→输出）+ 数学/单元验证（圆曲率 1/R 误差<1e-9、带状 Cholesky vs 稠密解残差<1e-12、已知拐点数 2k−1、解析面回复、确定性）+ 3 船完整算例 + 误差分析 + 一键复现（逐字节等价）。
- **L3 不成立**：全生态 `grep -r "lines_fairing\|fair_offset_table"` 于 TensorCAD/TensorPlatform/ship-performance-platform → **零调用方**（仅 `.git/worktrees/m18646/index` 命中）。无 CLI/registry/审计流水线接入，无平台 API 端点。
- **L4 不可裁**：`tensor_benchmark_cases` 本锚 0 行，无 ref/tolerance 判据行 → "基准裁判是唯一裁判"在本单**不存在裁判**；且 3 条船中 2 条的"解析真值"就是生成噪声表的同一参数面（kvlcc2 自述"NOT the official KRISO offset table"），属自洽/合成基准而非外部锚；手工光顺参考为替代件。
- **L5 不可达**（无外部实验/文献位形对拍、无跨机构复现）。

→ **级别 = L2**（"级别=功能实质"；报告篇幅不改变级别）。相较 m17622（GAP-CAD-A3-CURVATURE-FAIRNESS）判 L4 的依据是"主审计管线真调用 geometry_audit stage `smoothing_energy`"——本单无该集成，故不越级。

## 4 缺陷清单（按严重度，阻断与否）

- **F1（阻断落地，生态侧不可用）**：`git merge-base --is-ancestor 4fa7b2a origin/main` → **NO**；`git branch -a --contains 4fa7b2a` 仅 `agent/m18646`。交付只在工作树/本地分支，main 与 GitHub 均无。处置：随验收自动推送后复核，否则本条独立扣分。
- **F2（基准登记缺口，阻断 L4）**：`tensor_benchmark_cases` 本锚 0 行。建议种子行：`case_name=CASE-CAD-C3-FAIR-WIGLEY/series60`，metric=`max|faired-input| ≤ tol` 与 `total_flips 单调非增`，tolerance=`tol=1e-2/1e-2 m`，ref_source=`解析半宽面 + 本仓 run_example.py`。无此行则 L4 永远不可机器裁。
- **F3（报告欠账）**：`public.tensor_reports` 中 mission_id=18646 **无记录**（`report_id` 为空），即本单未发布平台论文级报告 → 守望线会派补报单并记欠账（本核验单不代偿，只记账）。
- **F4（数字夸大，非阻断）**：测试 docstring 内三处"实测"数字与代码可达结果不符，属过程流水式硬编码：
  - `test_energy_kernel_is_available_and_band_limited`：称 wigley energy `0.00216 → 0.00032`；实测 21×17 下 **0.00049（mp=4/bl=4）/ 0.00050（默认 5/5）**，夸大 ≈1.5×；
  - 同处称 series60 `energy 0.0048 / local 0.0011`；实测 21×17 为 **0.00630 / 0.00171**，25×21 为 0.00452 / 0.00149；
  - `test_fairing_moves_a_scattered_table_towards_the_analytic_surface`：称 kvlcc2 33×25 "同一运行改善到 0.00847"；实测 **0.00990（mp=4/bl=4）/ 0.01060（mp=5/bl=5，即 run_example 配置）**，改善率 5.6 %–11.8 % 而非 24.5 %。
  这三处不影响机器可裁条款①②，但违反"数字必须绑定可复现配置"的证据纪律，应在代码/docstring 中改为程序化注入。
- **F5（非阻断，结构限制）**：`buttock_levels=4` 时三船各只有 **1 条** iso-y 轮廓达到 `min_curve_points`（其余层级长度不足），故"三向光顺"的第三向覆盖率低；buttock-only 探针实测仍真移动型值（wigley 50/357、series60 201/525、kvlcc2 57/825 点），非空壳。分族预算软约束：buttock 层单调修复可把单点位移推到 0.0179 m > budget 0.0125 m（仍在总容差 0.05 m 内），文档"band_fraction 使末投影只是形式"的措辞略强。
- **F6（非阻断，自洽性）**：两条船(kvlcc2/series60)的"analytic"参照即生成面本身（`offset_forms` 自述 kvlcc2 非官方 offset），故降 RMS 只在 wigley 上具备外部锚意义。

## 5 检查项

已跑：34 项 pytest（独立重跑）；repro_cmd 全跑并与提交 results.json 比对（相同）；噪声-解析面 RMS/max 复算（6 组 kernel×hull）；kvlcc2 33×25 三组 (max_passes, buttock_levels) 复算；buttock-only 探针（3 船）；hand_fairing_reference 数值复算；docstring 数字复核；`tensor_benchmark_cases` 查询（110 行/本锚 0 行）；`tensor_reports` 查询（0 行）；git 祖先关系判定；全生态调用方检索。
未跑：`tests/` 全量回归（仅本锚 34 项）；loong64 集群复跑；其他船型/真实型值表（无数据）；人工手光顺对拍（无物）；推送到 GitHub 后的合并复核（无权限，按纪律不 push）。

## 6 处置

**request changes → gate 记 `reopened`**：
1. F1 推送/合并后复核 L3 前的可用性（本核验单不 push）；
2. F2 向 `tensor_benchmark_cases` 注册本锚 2 行判据（否则 L4 无裁判）；
3. F4 把三处 docstring 数字改为程序化注入（由 run_example 写出的 results.json 生成，禁止手抄）；
4. 若要上 L3/L4：把 `fair_offset_table` 接入 TensorCAD 的光顺审计管线（如 geometry_audit 的 `lines_fairing` stage）并留调用点证据；
5. F3 补发平台报告（本单发布的是核验报告，不是被核对象的交付报告）。