# 交差核验审计：m17352 · GAP-CAD-C1-5415-PARAM-REBUILD（TensorCAD）

审计人：m17434（board-audit）· 对象 = m17352 交付（**非重做**）· 结论：**L2 · reopened**
审计时点：2026-11-13 · 被审工作树 `/root/tp-agent/work/m17352`（交付文件全部为 untracked，未入任何 commit）

## 0. 独立性声明与手段
- 我不改被审代码；把 m17352 的 `tensorcad/ tools/ docs/` 复制到 `/tmp/aud17434` 后**从零重跑**两条复现命令，
  另写独立探针 `tools/audit_m17434/probe_gradient.py` 复核梯度、独立脚本复核距离/几何。
- 我**没有**重算 9001×27 数据库解析（只按 sha256 核源文件），也**没有**跑 9001 全配置留出（只复跑拟合脚本生成的 split）。

## 1. 三件套核验（逐项，带指针）
| 项 | 声明 | 核验 | 判 |
|---|---|---|---|
| run_artifacts | `runs/m17352_5415_param/{design_space.npz, design_space_holdout.csv, design_space_metrics.json, acceptance.json, sampled_configs.csv, gate_samples.csv, 3×png}` | 7/7 文件存在；复跑后 `acceptance.json` 584 个键**逐键相同**（仅 `wall_clock_s` 32.9→30.0），`design_space.npz` 数组 `np.array_equal` 全真、`design_space_metrics.json` 100 键仅 elapsed 不同 | REAL |
| benchmark ref vs got | 报告自称"配置0 逐点 PASS / 抽样 6/20 FAIL / 梯度 FAIL" | `tensor_benchmark_cases WHERE anchor_slug='GAP-CAD-C1-5415-PARAM-REBUILD'` → **0 行**；`tensor_benchmark_cases WHERE repo='TensorCAD'` → 仅 1 行（BASELINE-FLEET）；API `/tensor/bench/results?mission_id=17352` → `total:0`；`/tensor/bench/l4-verdict/17352` → `CLAIMED_NOT_PROVEN`（failed_checks 含 `ref_got_error_present/within_tolerance`） | **无基准条目 → L4 不成立** |
| repro_cmd | `PYTHONPATH=. python3 tools/dtmb5415_param/fit_design_space.py` + `run_acceptance.py`（32.9 s） | 两条均在 /tmp 副本跑通：拟合 35–88 s 输出逐字节一致，验收 30 s 输出逐值一致 | REAL |

三件套缺第二件（基准注册条目），其余两件经复跑为真。

## 2. 独立复核的数字（不是抄它的）
| 量 | m17352 声明 | 我复算 | 差 |
|---|---|---|---|
| STL→config0 网格点距 中位 | 0.0111 m | **0.01065 m**（去重前 494778 点，中位/梭 95 分位 0.0730） | 一致（其 82794 = 去重顶点） |
| 同 <分辨率占比 | 0.99902（门 >0.95） | **0.99922** | 一致，PASS 成立 |
| 反向 网格→STL 距离（它没做） | — | 中位 0.00874 m，p95 0.0709 m，max 0.462 m | 双向都 ~1 cm，几何确实重建得好 |
| 源数据 sha256 | 16b14ba5…bee2 | `sha256sum` 命中 | REAL |
| 参考 STL sha256 | d40b17b1…97fd | 命中，164926 面 | REAL |
| 留出集 RMS（shipped PCA-40 二次） | 2.45 %B（0.2516 m） | 复跑 metrics：0.25156 m / 2.4495 %B | 一致 |
| 对照 M1 线性 / M2 全二次 | 3.95 % / 2.36 % | 3.9528 % / 2.3636 % | 一致 |
| 27 变量体积梯度 | max rel err 1.0229（FAIL） | h=1e-3/1e-4/1e-5 **一律 1.023/1.037，且 g/fd 中位比 = 2.0001** | 数值可复现，但**诊断错了**（见 F1） |
| dV/dv 顶点雅可比 | 1.44e-8 | 矩阵链确为精确（与 FD 一致） | REAL |

## 3. 发现（按严重度）
**F1【阻断·可修】体积梯度是"整体系数 2"的 bug，不是结构断裂。**
`tensorcad/dtmb5415.py:562-576` 在散射前 `G *= 2.0`，而 `volume()`=2∫A dx 的镜像因子已由整船顶点散射承担；
实测 g/fd 比值 27/27 变量恒为 **2.0001**（跨 h 不变，`probe_gradient.py` 可复现）。删掉该行后 max rel err 由 **1.023 → 0.0115**（v=0），
v=全 0.2 时 0.0249：**门（<1%）仍差一点**——残余 ~1–2 % 才是它说的"散射合成"权重问题（我另测：把 x 项整体清零误差反涨到 2.3 %，说明 x 项方向对、权重仍欠准）。
即：honest-FAIL 的报告结论"未打通"方向对，但**根因判断与工作量估计错**（不是"预算所限未修完"，是漏了一个 ×2 后有 ~1 % 残差）。
补丁：`docs/audit_m17434/dtmb5415_gradient_factor2.patch`。

**F2【阻断·诚实性】水线 WATERLINE_Z=6.055 是魔法常数，且与它自己写的推导不成立。**
docstring 称"∫A dz = README 换算体积 = 0.549 t/998.5·(142/5.720)³ = 8412 m³"：该式实算 **8.412 m³（差 1000 倍，t→m³ 未换算）**；
我按同一判据反解水线得 **z = 6.137**（此处 ∇=8412），而常数 6.055 处 ∇=**8239.9（−2.05 %）**。
即"derived, not quoted"（引用的 `tools/dtmb5415_param/verify_geometry.py` **根本不存在**，docstring 悬空引用）不成立。

**F3【阻断·报告自相矛盾】** 同一物理量两处结论相反：§3.1 表把参考写成 `8.412 m³` / `Cb 0.00051`，算出 **+97853.77 % FAIL**；
§6 表同一行写"是（<3%）"而参考显示 `8 m³`。正确数（同 README 判据）= **−2.05 % / Cb −2.3 %，PASS**。
报告读者无法判定到底是 FAIL 还是 PASS；且这个"灾难性 FAIL"会污染聚合看板。

**F4【中】hydrostatics 里有一个硬编码量**：`tensorcad/dtmb5415.py:615 "keel_z_m": -3.279432` 是常数（取自参考 STL），
不是参数化模型输出（模型云 z_min = **−3.0245**，与 STL 差 **0.2550 m**）；换任意 v 都返回同一个舷龙骨高 → MOCK 字段。

**F5【中】验收门 1 的分辨率自指**：门限 = 该网格自身边长中位（0.455 m），故"99.9 % 点在门内"对粗网格几乎自动成立；
它真正有力的证据是中位 1.07 cm（≈2.3e-4·B），但这被门限噪声掩盖。另：config0 云/网格 z_min 比参考 STL 浅 0.255 m（龙骨/艉鳍区），
门限过粗不报警——几何先验证 Cb/L/B 用 Cb 时也受影响。

**F6【中】抽样门把结论写成了"验收 2 PASS/FAIL"，实质是模型形式误差**：20 抽样 6/20、均匀抽样 0/20（失败项 100 % 集中在
`half_breadth_non_negative`，即远离基准的 v 会让舷侧穿过中线）。报告已如实写 FAIL（无粉饰），但**验收 2 的门本身（在 [-1,1]^27 全域
要求网格质量）与"回归代理 2.45 %B"的能力不相容**，重做同一模型永远过不了 → 换代目标应是"收紧设计域 + 声明适用包络"。
**F7【轻】无单元测试**：全仓无 `tests/*5415*`；`tests` 槽为空，可复现性由我复跑验收脚本承担。
**F8【轻】无集成调用方、无 commit**：`grep -rn dtmb5415` 仅命中自身与 tools；工作树 4 个交付目录全为 `??`，未入分支 → 平台侧不可达。

## 4. 级别判定
- L0/L1 超越：有真代码、真数据、数值可复现。
- **L2 = 成立**：`tensorcad/dtmb5415.py`（727 行）真输入(27 维 v)→计算(点云/网格/水线静水力/体积/雅可比)→输出；
  配置 0 子案与独立外部参考（平台 `DTMB5415_Hull_exact.stl`，sha256 核对）双向 ~1 cm；整条流水线可复现到逐值一致。
- L3 **不成立**：零调用方（无 API、无求解器接入、无 tests 引用）。
- L4 **不成立**：`tensor_benchmark_cases` 0 行 + `/tensor/bench/results` 0 行 → 无 ref/got/tolerance 机器判据；
  且三条预置验收中 2.5 条未过（抽样 6/20、均匀 0/20、梯度 102 %）。
- L5 不成立。百页报告不改级别。

## 5. 裁定（disposition）
**request changes（reopened）** —— 不是走过场"已成功"。可在不改模型的前提下修好的三件事：
1. 删 `G *= 2.0`（F1）并把残余 1–2 % 的权重定位清楚；把 docstring"matches to round-off"改成实测值。
2. 用 z=6.137 或把「水线=README 体积反解」实现成 `verify_geometry.py` 真脚本（F2/F3），修 §3.1/§6 单位与结论。
3. 删掉 `keel_z_m` 硬编码，改由 `v` 计算（F4）；把 keel 0.255 m 差异写进局限。
随后（可选）登记 `tensor_benchmark_cases` 条目（config0 点距中位 ≤1.5 cm、留出 RMS ≤3 %B、∇ within 3 %）才有 L4 通道。

## 6. 未做的检查（unknown，不是 pass）
- 未验证 9001 全配置网格质量（只跑预置 20+20 抽样）；未验证 `.mat` 解析的 2250 点布局对**其它列**是否同样成立（只核了第 0 列不变量的代码路径）。
- 未做与文献 Cb=0.506 的第三方独立对标（我用的是数据集 README 表 8 412 m³ 这条链）。
- 未跑 `--patch-scatter` 之外的第二类根因修复实验（x 项权重再参数化）；残余 1–2 % 的精确来源仍 unknown。
