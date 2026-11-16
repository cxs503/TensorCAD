# GAP-CAD-D2 几何金标回归门禁（TensorCAD 落点）

> 本文由 `tools/golden/gap_d2_gate_docs.py` 从 `artifacts/gap-cad-d2-golden/*.json` 与两个
> 门禁文件本身程序化生成，正文数字零手抄（`--check` 可在 CI 里查陈旧）。
> 上游：TensorPapers 报告 **#20224**，套件 `papers/gap-cad-d2-golden-regression/`。
> 锚点：`GAP-CAD-D2-GOLDEN-REGRESSION`（任务 m22501，TensorCAD 仓侧收口）。

## 1. 这道门是什么，为什么仓里现在有两道

TensorCAD 的几何是生态的地基：`hull` / `wigley` / `series60` / `dtmb5415` 产形状，
`hydrostatics` + `mesh` + `geometry_audit` 把形状变成 Cb、∇、VCB、BM_T、湿表面积。
这些量漂 1 % 不会有任何测试报警，只会让下游阻力/耐波算例悄悄错。

仓里原有 `.github/workflows/golden-geometry.yml`（`tensorcad.golden`）守的是**仓内自己
committed 的基线**：四个金标 case、`FROZEN` 档判"是否漂移"。它的弱点写在它自己的设计里——
基线由本仓 `--bless` 生成，所以一次"改了代码顺手重 bless"的提交可以让它全绿。

本单装进来的 `geometry-golden.yml` 是**第三方的判据**：判据代码与金标数值都不住在这个仓里
（TensorPapers 的 dual-derived suite），它对该仓只暴露 `TENSORCAD_PATH`。每个金标值在进入
判定前必须由两条结构不同的推导给出同一结果（EXACT 档 1e-9 / QUAD 档 1e-6），否则直接抛
`GoldenError`，**不存在手敲回退值**。因此这道门不可能被"bless 成与自己一致"。

两道门并存不是冗余：`golden-geometry.yml` 便宜（~5 s，粗网格）、每次 `tensorcad/**` 改动都跑；
`geometry-golden.yml` 贵（~30 s，两档网格 + 检测力探针）、只在这一改动**能移动几何数字**时跑。

## 2. 触发路径（程序化取自 `geometry-golden.yml` 与 `.githooks/pre-commit-geometry-golden.sh`）

| 路径 | push 触发 | PR 触发 | pre-commit 钩子 |
|---|---|---|---|
| `tensorcad/hull.py` | yes | yes | yes |
| `tensorcad/wigley.py` | yes | yes | yes |
| `tensorcad/series60.py` | yes | yes | yes |
| `tensorcad/dtmb5415.py` | yes | yes | yes |
| `tensorcad/hydrostatics.py` | yes | yes | yes |
| `tensorcad/hydrostatics_analytic.py` | yes | yes | yes |
| `tensorcad/mesh.py` | yes | yes | yes |
| `tensorcad/topology.py` | yes | yes | yes |
| `tensorcad/exchange.py` | yes | yes | yes |
| `tensorcad/geometry_audit.py` | yes | yes | yes |
| `tensorcad/voxel.py` | yes | yes | yes |
| `.github/workflows/geometry-golden.yml` | yes | yes | no (CI only) |

**不触发的改动**：只改 `docs/**`、`README.md`、`CHANGELOG.md`、`tests/**` 等文档/测试文件时，
push/PR 路径过滤器不命中，pre-commit 钩子也在路径判断上直接退出（实测 0.02 s，见 §6）。
这条性质由 `tests/test_geometry_golden_gate.py::test_documentation_only_edits_are_not_watched`
与 `::test_hook_lets_a_documentation_only_commit_through_without_running_the_suite` 双保险钉住。

一个诚实的边界：**PR 路径过滤器本身无法在本仓 CI 里被自测**（一次 docs-only 的 PR 不触发
工作流，也就没有工作流能报告"我没触发"）。因此"该不触发时不触发"是靠单元测试对
`on: push/pull_request` 路径表做的静态断言，而不是靠 CI 跑绿——两种证据强度不同，这里写明。

## 3. 容差表（本次实测）

各 case 的整体判据：

| case | 族 | 档位 | 容差（相对） | max 实测误差 | 判定项 | WITHHELD |
|---|---|---|---|---|---|---|
| `wigley_coarse` | Wigley | coarse | 1.00e-04, 3.66e-04, 0.001055, 0.001172, 0.001172, 0.003160, 0.004329, 0.004329 | 1.95e-04 | 8 | 0 |
| `wigley_fine` | Wigley | fine | 1.00e-04, 3.66e-04, 0.001055, 0.001172, 0.001172, 0.003160, 0.004329, 0.004329 | 1.95e-04 | 8 | 0 |
| `series60_coarse` | Series-60 | coarse/invariant | 1.00e-09, 1.00e-04, 2.06e-04, 2.06e-04, 0.018110, 0.018312, 0.018312 | 0.008591 | 8 | 0 |
| `series60_fine` | Series-60 | fine/invariant | 1.00e-09, 1.00e-04, 2.06e-04, 2.06e-04, 0.018110, 0.018312, 0.018312 | 0.008591 | 8 | 0 |
| `sphere_ellipsoid_coarse` | SDF primitives | coarse | 1.00e-04, 0.002776, 0.002816, 0.005457, 0.005457 | 8.86e-04 | 10 | 0 |
| `sphere_ellipsoid_fine` | SDF primitives | fine | 1.00e-04, 0.002776, 0.002816, 0.005457, 0.005457 | 8.86e-04 | 10 | 0 |
| `dtmb5415_coarse` | DTMB 5415 | coarse/declared/invariant | 0, 1.00e-09, 1.00e-04, 0.030000 | 0.015693 | 11 | 0 |
| `dtmb5415_fine` | DTMB 5415 | declared/fine/invariant | 0, 1.00e-09, 1.00e-04, 0.030000 | 0.015693 | 11 | 0 |

逐量视图（同一量在粗/细网格上的容差不同，故给区间）：

| 量 | 符号 | 档位 | 容差区间（相对） | 全 case 最坏误差 | case 数 |
|---|---|---|---|---|---|
| `block_coefficient` | Cb | coarse/declared/fine | 0.004329 – 0.030000 | 0.015693 | 6 |
| `displacement_volume` | ∇ | coarse/declared/fine | 0.004329 – 0.030000 | 0.008591 | 6 |
| `ellipsoid_extent_x` | ellipsoid_extent_x | coarse/fine | 1.00e-04 – 1.00e-04 | 0 | 2 |
| `ellipsoid_extent_y` | ellipsoid_extent_y | coarse/fine | 1.00e-04 – 1.00e-04 | 0 | 2 |
| `ellipsoid_extent_z` | ellipsoid_extent_z | coarse/fine | 1.00e-04 – 1.00e-04 | 0 | 2 |
| `ellipsoid_surface_area` | ellipsoid_surface_area | coarse/fine | 0.002776 – 0.002776 | 4.50e-04 | 2 |
| `ellipsoid_volume` | ellipsoid_volume | coarse/fine | 0.005457 – 0.005457 | 8.86e-04 | 2 |
| `holdout_rms` | holdout_rms | coarse/fine | 1.00e-04 – 1.00e-04 | 2.65e-08 | 2 |
| `identity_cb_equals_cp_times_cm` | identity_cb_equals_cp_times_cm | invariant | 1.00e-09 – 1.00e-09 | 0 | 2 |
| `identity_cb_equals_volume_over_lbt` | identity_cb_equals_volume_over_lbt | invariant | 1.00e-09 – 1.00e-09 | 0 | 2 |
| `identity_mirror_doubles_y_extent` | identity_mirror_doubles_y_extent | invariant | 1.00e-09 – 1.00e-09 | 0 | 2 |
| `identity_mirror_preserves_x_extent` | identity_mirror_preserves_x_extent | invariant | 1.00e-09 – 1.00e-09 | 0 | 2 |
| `identity_mirror_preserves_z_extent` | identity_mirror_preserves_z_extent | invariant | 1.00e-09 – 1.00e-09 | 0 | 2 |
| `identity_open_shell_has_one_boundary_component` | identity_open_shell_has_one_boundary_component | invariant | 0 – 0 | 0 | 2 |
| `identity_open_shell_is_not_a_solid` | identity_open_shell_is_not_a_solid | invariant | 0 – 0 | 0 | 2 |
| `identity_volume_equals_cb_lbt` | identity_volume_equals_cb_lbt | invariant | 1.00e-09 – 1.00e-09 | 8.88e-16 | 2 |
| `length_stations` | length_stations | declared | 0.030000 – 0.030000 | 0.008908 | 2 |
| `longitudinal_centre_buoyancy` | LCB | coarse/fine | 1.00e-04 – 1.00e-04 | 1.78e-16 | 4 |
| `maximum_section_area` | A_mx | coarse/fine | 0.003160 – 0.018110 | 0.008557 | 4 |
| `prismatic_coefficient` | Cp | coarse/fine | 2.06e-04 – 0.001172 | 3.91e-05 | 4 |
| `sphere_extent_x` | sphere_extent_x | coarse/fine | 1.00e-04 – 1.00e-04 | 0 | 2 |
| `sphere_extent_y` | sphere_extent_y | coarse/fine | 1.00e-04 – 1.00e-04 | 0 | 2 |
| `sphere_extent_z` | sphere_extent_z | coarse/fine | 1.00e-04 – 1.00e-04 | 0 | 2 |
| `sphere_surface_area` | sphere_surface_area | coarse/fine | 0.002816 – 0.002816 | 4.56e-04 | 2 |
| `sphere_volume` | sphere_volume | coarse/fine | 0.005457 – 0.005457 | 8.86e-04 | 2 |
| `vertical_centre_buoyancy` | VCB | coarse/fine | 0.001055 – 0.001055 | 5.21e-05 | 2 |
| `waterline_beam` | waterline_beam | declared | 0.030000 – 0.030000 | 0.011848 | 2 |
| `waterplane_area` | A_wp | coarse/fine | 2.06e-04 – 0.001172 | 3.91e-05 | 4 |
| `wetted_surface_area` | S_wet | coarse/fine | 3.66e-04 – 3.66e-04 | 1.80e-05 | 2 |

本次判词分布：PASS=74。

`FROZEN` 档出现在 `dtmb5415_exact` 上：三角面数/顶点数/边界边数是**整数相等**（容差 0），
包围盒三向尺度容差 1e-9——它判的是"出厂工件还是不是那一个"，不是"精度有多好"，所以这一行
的 max 误差是 0 也不构成任何精度声明。

**容差是怎么来的（不是拍的）**：EXACT 档容差 = 该 case 网格离散误差的实测上界留 2 倍以上余量；
`series60_cb060_default`（?x? 默认网格）给到 4.5e-2，是因为该网格对 Cb 的欠填实测 1.7747%（`series60_convergence.json` 首档），
把"粗网格的已知缺陷"显式写成容差，而不是悄悄调大或把 case 删掉。细则见 TensorPapers
报告 #20224 与 `outputs/series60_convergence.json`。

## 4. WITHHELD：不能判的量就不判

_（本次没有 WITHHELD。）_

`WITHHELD` 永不产出判词、也永不被"补齐成 PASS"：5415 出厂 STL 是开口面（边界边存在），
无独立锚，所以 ∇/Cb 明写拒判；Series 60 的 Cm/Cp/VCB/Cwp/湿表面积没有表值转录，同样拒判。
只有 Cb/∇ 以**定义式**入判（积分必须回到请求的丰满度）。

## 5. 检测力（门有没有牙齿）

控制臂：先证明这道门在被测库"原样"时是全绿的——否则后面的"抓到变异"毫无意义。

| case | load-bearing 控制臂 |
|---|---|
| `dtmb5415_coarse` | PASS |
| `dtmb5415_fine` | PASS |
| `series60_coarse` | PASS |
| `series60_fine` | PASS |
| `sphere_ellipsoid_coarse` | PASS |
| `sphere_ellipsoid_fine` | PASS |
| `wigley_coarse` | PASS |
| `wigley_fine` | PASS |

共 8/8 控制臂 PASS。

探针（把被判量挪到"总误差 = k × 容差"）：

| 探针组 | 抓/应抓 | 命中率 | 说明 |
|---|---|---|---|
| 强变异族（k>1，必须 FAIL） | 24/24 | 1.0000 | **验收口径的"检测力 100 %"** |
| 几何闸探针（L/B/T 错船型） | 12/12 | 1.0000 | 闸不过即降级 WITHHELD |
| WITHHELD 拒判探针 | 3/3 | 1.0000 | 丢证据的量必须被拒判 |
| 弱变异族（k≤0.5，**不要求**抓到） | 6/12 | 0.5000 | 容差内的自然离散，抓到才算过度敏感 |
| 误报 | 0 例（`false_alarms=0`） | — | 控制臂 + k≤0.5 臂无一误判 |

`unexpected_misses=0`、`misses=[]`。
注意 `quantity_detection_rate` 只有 0.5000：那一栏含 k=0.5 的
**容差内**探针，抓不到才是正常（抓到反而说明容差虚紧）。验收要求的 100 % 指
`strong_detection_rate` 与 `gate_detection_rate`，两者都是 1.0000。

### 5.1 退出码不够：一个"测不出来的金标 case"必须算红

装门时做了一次**真回归**对照，不只做合成变异：把 TensorCAD 复制一份、在源码里注入一处改动，
再用同一个钩子去判它（复现：`tensorcad/wigley.py` 的 `half_breadth` 把 `0.5*beam` 改成
`0.515*beam`；`tensorcad/hydrostatics.py` 把水线面投影面积乘 1.02）。

| 注入的回归 | `suite.py` 退出码 | 未 PASS 的 case | 门怎么判 |
|---|---|---|---|
| Wigley half-breadth +3% | 0 | `wigley_coarse`=SKIP, `wigley_fine`=SKIP | `suite.py` 退出 0，被报告断言判红（钩子退出 1：`[golden] BLOCKED: a golden case was not measured`） |
| Series-60 Cb +5% | 1 | `series60_coarse`=FAIL, `series60_fine`=FAIL | `suite.py` 退出 1（FAIL 12 项）→ 第一道即红 |
| Sphere radius +10% | 1 | `sphere_ellipsoid_coarse`=FAIL, `sphere_ellipsoid_fine`=FAIL | `suite.py` 退出 1（FAIL 10 项）→ 第一道即红 |

第一行是本单的一个**实质发现**：几何模块改坏后 TensorCAD 自己的 `generate_wigley` 会拒绝出网格
（水线偏移偏离解析面 3.0e-2 m），套件把该 case 记成 `SKIP` / `UNAVAILABLE`——这符合它"缺测
永不 PASS"的诚实口径——**但 `suite.py` 仍然退出 0**。也就是说，只看退出码的门会在"两个金标
case 已经不再被测量"的情况下放行。本仓两个载体因此都不读退出码，而读报告：
`tools/golden/gap_d2_assert_report.py` 逐 case 断言 `verdict` 以 `PASS` 开头且 case 数 = 8，
再断言检测力的四个率；钩子与 CI 共用这一个脚本（`tests/test_geometry_golden_gate.py` 里钉住）。
上游 `suite.py` 的退出码语义未被改动——本单的写权限只在 TensorCAD 侧，加固做在这里。

## 6. 首跑证据（本仓 checkout，非 TensorPapers 侧）

| 脚本 | 退出码 | 墙钟 | 环境 |
|---|---|---|---|
| `suite.py` | 0 | 45.2 s | `TENSORCAD_PATH=/root/tp-agent/work/m22681` |
| `mutants.py` | 0 | 120.5 s | `TENSORCAD_PATH=/root/tp-agent/work/m22681` |

- 被检 checkout：`/root/tp-agent/work/m22681`，HEAD `45e4c3ade0dd`，Python 3.11.6
- 墙钟随机器负载波动（同一条命令在 §6 的钩子探针里实测 35.00 s）；可复现的是退出码与
  `run_meta.json` 里记的 sha256，不是秒数
- 两次运行原始输出：`artifacts/gap-cad-d2-golden/suite.log`、`mutants.log`
- 快照（供 §3/§5 表格注入）：

| 快照（`artifacts/gap-cad-d2-golden/`） | sha256 前 16 位 |
|---|---|
| `golden_suite.json` | `242489e38fcb4265…` |
| `hook_probe.json` | `b6eeae8c4d82a8a8…` |
| `mutation_report.json` | `1a022f0ae6f3c6ba…` |
| `regression_control.json` | `c0800aba1de54cb4…` |
| `series60_convergence.json` | `b5c2a6c14195482e…` |

pre-commit 链步实测（`git diff --cached` 语义用 `GOLDEN_STAGED_FILE` 复现，
`tools/golden/gap_d2_gate_probe.py` 生成 `artifacts/gap-cad-d2-golden/hook_probe.json`）：

| 暂存内容 | 钩子行为 | 退出码 | 墙钟 | 末行 |
|---|---|---|---|---|
| 仅暂存 `docs/*.md`、`README.md`（走仓内钩子链 `.githooks/pre-commit`） | 路径过滤器两道都不命中，直接返回，套件不启动 | 0 | 0.02 s | `(no geometry change detected)` |
| 同上，但只调 dual-derived 钩子 | 路径不命中，直接返回（0.02 s 量级：钩子本身的开销） | 0 | 0.02 s | `(no geometry change detected)` |
| `tensorcad/hull.py`（只调 dual-derived 钩子） | 跑 suite.py + mutants.py，全绿才放行 | 0 | 35.00 s | `[golden] ok` |
| `tensorcad/hull.py`（走仓内钩子链） | 先跑仓内 `tensorcad.golden --quick`，再跑 dual-derived suite + mutants | 0 | 40.00 s | `[golden] ok` |
| `tensorcad/hull.py` + 一份被注入 +3% Wigley 半宽的 TensorCAD 副本 | **`suite.py` 退出 0，仍被报告断言判红**（金标 case 变成 UNAVAILABLE） | 1 | 35.00 s | `[golden] BLOCKED: a golden case was not measured` |
| `tensorcad/hull.py` + `GOLDEN_DIR` 指向不存在的目录 | **明确拒绝**（不静默降级成放行） | 1 | 0.01 s | `refusing to guess` |
| `tensorcad/hull.py` + `GOLDEN_DIR` 未设、旁无 TensorPapers | 警告 + 放行（本地不可执行 ≠ 本地失败），CI 侧仍 fail-closed | 0 | 0.01 s | `this commit is not locally gated` |

## 7. 怎么用

```bash
# CI：push/PR 命中几何路径即自动跑（.github/workflows/geometry-golden.yml）
# 本地：钩子由仓内既有 hooksPath 安装器接管，无需另设 .git/hooks
tools/golden/install_git_hooks.sh          # core.hooksPath=.githooks（本仓既有机制）

# 套件位置解析顺序：$GOLDEN_DIR → <repo>/../TensorPapers/papers/gap-cad-d2-golden-regression
#                     → $HOME/TensorPapers/papers/gap-cad-d2-golden-regression
GOLDEN_DIR=/path/to/TensorPapers/papers/gap-cad-d2-golden-regression git commit ...

# 手动跑（与 CI 完全同一条命令）
cd $GOLDEN_DIR && TENSORCAD_PATH=$(git -C /path/to/TensorCAD rev-parse --show-toplevel) \
    PYTHONPATH=$(git -C /path/to/TensorCAD rev-parse --show-toplevel) python3 suite.py
```

`GOLDEN_DIR` 显式给出但目录不存在 = 硬错（不静默降级）；未给出且找不到 = 本地放行 + 警告，
因为"这个 checkout 旁边没克隆 TensorPapers"不等于"代码回归了"。CI 侧没有这个宽容：
拿不到 TensorPapers 或套件失败都会让 job 红。

## 8. 局限与改动清单（诚实）

- **只证明几何/静水力**：这道门对 Cb/∇/VCB/BM_T/S_wet 一类量负责，**不**证明阻力、耐波、
  操纵、结构适用性。几何对不等于船好。
- **控制臂 8/8 不等于全谱覆盖**：8 个 case 覆盖球/椭球/Wigley/Series60/5415，Series60 的
  高阶级系数、5415 的 ∇ 仍是拒判（无独立锚）。
- **粗网格容差是承认缺陷，不是校准**：`series60_cb060_default` 3.7 % 的 Cb 欠填被写进容差，
  细网格 case（321×81）才是精度声明。任何人想"把粗网格也收进 2 %"，得改套件而不是改容差。
- **runner 标签**：drop-in 原文写 `runs-on: ubuntu-latest`，本仓其余门禁用自托管 `ci-loong64`，
  这里改成 `ci-loong64`，否则该 job 在本生态里不会被派发。这一处与其他两处改动
  （TensorPapers checkout 的 token 回退、验收断言步）都写在 workflow 头部注释里。
- **上游退出码语义有口子（已在本仓兜住）**：`suite.py` 在 case 变 `UNAVAILABLE` 时仍退出 0。
  建议 TensorPapers 侧把 `UNAVAILABLE` 计入非零退出；在它改之前，TensorCAD 侧靠报告断言兜住
  （§5.1 的第一行就是这条口子的实测证据）。
- **钩子不是强制手段**：`git commit --no-verify`、或没装 hooksPath 的 checkout 都能绕过；
  权威判据是 CI。
- **`.git/hooks/pre-commit` 未直接写入**：本仓用 `core.hooksPath=.githooks`，钩子以
  husky 式链步（`.githooks/pre-commit` 末尾调用 `.githooks/pre-commit-geometry-golden.sh`）
  落地，这样改动是**可 review 的仓内文件**，而不是别人机器上的 .git 目录。

## 9. 相关文件

| 文件 | 作用 |
|---|---|
| `.github/workflows/geometry-golden.yml` | CI 门（权威，fail-closed），含验收断言步 |
| `.githooks/pre-commit-geometry-golden.sh` | 本地 husky 式链步（路径过滤 → suite + mutants） |
| `.githooks/pre-commit` | 仓内第一道（`tensorcad.golden --quick`）并链第二道 |
| `tests/test_geometry_golden_gate.py` | 两臂路径表一致性、docs-only 不触发、钩子行为 |
| `tools/golden/gap_d2_assert_report.py` | 判定层：读两份报告断言验收数字（钩子与 CI 共用） |
| `tools/golden/gap_d2_gate_probe.py` | 钩子七种决策模式的实测探针（写 `hook_probe.json`） |
| `tools/golden/gap_d2_gate_docs.py` | 本文生成器（`--check` 查陈旧） |
| `artifacts/gap-cad-d2-golden/regression_control.json` | 真回归对照（注入改动 → 套件与门怎么判） |
| `artifacts/gap-cad-d2-golden/` | 首跑快照：suite/mutants JSON、日志、run_meta |
