# m17622 交差核验：GAP-CAD-A3-CURVATURE-FAIRNESS（对象 m17465）

裁判口径：本单只核验 m17465 的交付是否真存在、数值是否在容差内、级别是否虚标。级别=功能实质，
百页报告不改变级别。核验时间 2026-11-13T07:20–07:30Z，审阅者 board-audit。

## 1 案例介绍（被审对象）

- 锚点验收条款（来自 m17465 goal 原文）：①解析球主曲率处处 1/R，误差 <1e-6；②Wigley 光顺能量数值稳定；
  ③与 meshmagick 对拍 3 例。
- 物理/几何对象：单位球族 R=0.1…10 m（7 档，闭式 k1=k2=1/R）；Wigley 解析舷侧（L=10 B=1 T=0.625，n=m=2，
  Gauss-Legendre 解析积分 E_bending=21.537877 / E_willmore=5.384469 / 湿面积 14.879063）；
  参数化船体（L=10 B=2 T=0.8，nx=41 nz=13，排水体积 7.446554 m³，Cb=0.4654）。
- 待核验量：三条光顺能量 E_bending=∫(k1²+k2²)dA、E_willmore=∫H²dA、E_gauss=∫K dA，及其离散不变性。

## 2 怎么做的（核验方法/流程）

1. **取证**：读 m17465 工作树 `/root/tp-agent/work/m17465`、tensor_reports id=16448、run artifacts。
2. **基准裁判**：`SELECT case_name,metric,tolerance … WHERE anchor_slug='GAP-CAD-A3-CURVATURE-FAIRNESS'`
   → **0 行**（全表 55 行无任何 CAD 锚）。故裁判改为「锚点验收条款 + 交付自带 gates_verified.json」，
   并逐行独立重算，不采信报告文字。
3. **独立复跑**（不覆盖被审工件，输出落 /tmp 与本单 runs/）：
   - `python benchmarks/verified/curvature_fairness/run_benchmark.py --out /tmp/audit17622/bench` → 46.7 s，exit 0。
   - `python -m tensorcad.fairness_cli sphere|hull`、`python -m tensorcad.cli fairness wigley` → exit 0。
   - `python tools/crosscheck_curvature_thirdparty.py --out /tmp/audit17622/crosscheck.json` → exit 0。
   - `python -m pytest tests/test_fairness_energy.py tests/test_fairness_entry.py -q` → **29 passed, 112.07 s**。
4. **工件完整性**：按 runs/m17465/manifest.json 重算 sha256 → **13/13 一致**（含字节数）。
5. **落地态核验**：`git ls-remote` + `git merge-base --is-ancestor` + GitHub API PR 列表。

## 3 计算结果（裁判逐项，程序化比对，零手抄）

| 条款/门 | 参考（ref） | 实测（got） | 误差 | 容差 | 判定 |
|---|---|---|---|---|---|
| G1 球主曲率 R=0.1…10（7 档） | 1/R（闭式 10→0.1 1/m） | 10.000000000000385 … 0.09999999999998727 | 审阅者另算 **max&#124;k−1/R&#124;=2.25e-7 1/m**(R=0.1) | 1e-6 绝对 | **PASS** |
| G2 球能量·场路径 level24 | 25.132741 (8π) / 12.566371 (4π) | 25.156293 / 12.578147 | +0.0937 % | 0.1 % | PASS |
| G3 球能量·网格路径 level32 | 同上 | 26.573625 / 12.763146 / 12.239480 | +5.733 / +1.566 / −2.601 % | 6 % | PASS |
| G4 Wigley 光滑桶 @121×37 | 21.537877 解析舷侧 | 20.084938（E_bending） | **−6.746 %** | 15 % | PASS |
| G4' same | 5.384469 | 5.036300（E_willmore） | −6.466 % | 15 % | PASS |
| G5 Richardson 外推 | 21.537877 | 21.550310 | **+0.0577 %** | 0.5 % | PASS |
| G5' 参考自检（阶数翻倍） | 0 | 9.90e-16 | 绝对 | 1e-12 | PASS |
| G6 不变性（缩放/单元序/1.6 km 平移） | 0 | 0 / 2.27e-13 / 1.14e-11 | 舍入级 | 1e-10 | PASS |
| G7（诊断）椭球 K 中位偏差 | 解析 K | 22.43 % | **> 15 % 线** | 0.15 | **FAIL（如实保留）** |
| G7（诊断）球 H / K | 解析 | 1.58 % / 0.05 % | — | 0.15 | PASS |

复跑一致性：20/20 门行与交付 result.json **逐位相同**（含全部误差百分数，0 处不符）。
验收条款兑现度：① **达**（且用更严的 max 判据仍 <1e-6）；② **达**（G4/G5/G6：收敛自下方 −39.98→−20.34→
−10.15→−6.75 %，外推 +0.058 % 落在解析参考上，不变性 ~1e-11）；③ **未达/诚实 WITHHELD**
（meshmagick 三次尝试均失败，provenance 逐条入库，替代为 trimesh 法向环+闭式解，且椭球 K 22.43 % 破 15 %）。

## 4 结果图

`docs/board_audit/m17622/refs/` 与 m17465 交付内 `benchmarks/verified/curvature_fairness/figures/curvature_fairness.png`
三面板：球场路径 |E_bending−8π| 随 level 收敛；Wigley 光滑桶收敛 vs 折边桶发散（对数纵轴+解析参考线）；
三例外部位对拍中位偏差柱状图（红虚线=15 % 标准线）。本单不重绘（审计不改动被审工件）。

## 5 注意点（核验发现的缺陷/摩擦，按严重度）

- **F1 落地缺口（阻断生态收益）**：远端 main=5eae547 **不含**交付；交付在 `agent/m17465`=e030c60
  （PR **#422 审计时仍 open**）。`git merge-base --is-ancestor e030c60 origin/main` → NO。
  且推送门禁滤除 16 个文件，`runs/m17465/*`（repro_cmd.sh、sha256 manifest、图、CLI 记录）**只在本地工作树**，
  分支上无。→ 代码实质 L4，但「未合并 = 生态侧不可用」，合并前有效级别应视为 L0。
- **F2 基准登记缺口（平台级）**：`tensor_benchmark_cases` 该锚 0 行 → G3/G4/G6 的容差由作者自定，
  唯一来自条款的容差是 1e-6。建议把 G1/G5 行种子进表，使「基准判据是唯一裁判」真正可机器裁。
- **F3 条款③未兑现**：meshmagick 对拍 WITHHELD，替代件可信（Cohen-Steiner & Morvan 法向环 + 闭式解）。
- **F4 读数风险**：`all_gates_pass=true` 与 `diagnostics.G7.verdict=fail` 并存，只看旗标的消费者会误读「条款全达」。
  建议把「门」与「条款」分开两个字段。
- **F5 已澄清的加强结论**：G1 门用探针均值，审阅者另算 max 偏差仍 2.25e-7 < 1e-6，条款①不因均值口径而虚达标。
  报告中「≤2.2e-7 1/m」的写法经复算成立。

## 6 误差判定（逐项回答「是否<3%」）

| 关键量 | 误差 | <3 %？ |
|---|---|---|
| 球主曲率（max 口径） | 2.25e-7 1/m（相对 ~2e-8） | ✅ 是 |
| 球能量·场路径 | +0.094 % | ✅ 是 |
| Wigley Richardson 外推 E_bending | +0.058 % | ✅ 是 |
| 不变性残差 | ≤1.1e-11 绝对 | ✅ 是 |
| 球能量·网格路径（G3） | +5.73 % | ❌ 否（≤声明容差 6 %，网格离散固有） |
| Wigley 光滑能量直接网格（G4） | −6.75 % | ❌ 否（分辨率所限，靠外推消除） |
| 椭球 Gaussian 曲率（G7 诊断） | 22.43 % | ❌ 否（**破 15 % 标准**，如实 FAIL） |

结论：**可被参考值判定的量全部 <3 %；两处 5–7 % 与一处 22.4 % 如实保留、未粉饰。**

## 判级与处置

- **级别 L4（基准验证）**：impl_files 非占位 + 入口可跑 + run_artifacts 真存在且 sha256 自洽 +
  主审计管线真调用（geometry_audit stage `smoothing_energy`）+ ref vs got 数值对照 + 一键复现可执行，
  六项俱全并已独立复跑。**不判 L5**：外部位形（meshmagick）WITHHELD，且椭球 K 破 15 % 标准。
- 处置：**approve with notes**；随附 F1（合并落地，PR #422）、F2（基准入表）两条返修建议。
- 未跑检查如实列出：全量 262 测试、其他船型族 CLI、集群 loong64 架构复跑、meshmagick 对拍。