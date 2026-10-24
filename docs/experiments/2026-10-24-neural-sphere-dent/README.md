# 神经几何球体局部编辑：真实运行记录（2026-10-24）

本目录是一次**在本机真实执行**的神经几何实验记录，不是对代码或文档的复核。所有数字都来自
`examples/neural_geometry/` 中的脚本，并附有独立于脚本自身度量的复核结果、权重文件、
原始 stdout、环境清单与复现命令。

## 0. 先纠正一个状态前提

任务描述里的时间线已经过时：

| 事件 | 实际状态 |
|---|---|
| PR #8（revisioned online edit） | **已于 2026-10-24 03:46 UTC 合并**，head `93da961`；其 `tests` CI 最初两次失败（19 s），随后在分支上修复并转绿 |
| 提交 `7d4fdd0`（非有限指标写成 JSON `null`） | **已随 PR #9 于 04:00 UTC 合并进 main**，不是待验证状态 |
| PR #9（回滚报告修复） | 已合并 |
| PR #10（启用神经几何 CI） | 已合并；其后的 main 上神经几何任务**失败过一次** |
| PR #11（补 NumPy 依赖、凹陷采样贴近表面） | 已合并（04:50），main 上的神经几何任务随后**成功** |

所以"修复并确认 PR #8 的 CI"这一步已经无法按原样执行（PR #8 早已合并，而且它合并时
`neural-geometry-experiments` workflow 还不存在，是 PR #10 才加的）。`research/neural-geometry-online-edit`
分支也已被删除。**当前 main 上的代码已经包含 `7d4fdd0` 的内容**，本次运行正是对它的
真实检验（见第 4 节的回滚用例）。PR #8 也谈不上"暂不合并"，需要的话只能 revert。

## 1. 环境

| 项目 | 值 |
|---|---|
| Python | CPython 3.11.2 |
| PyTorch | 2.2.2+cu121（CPU 运行，`torch.cuda.is_available() == False`） |
| safetensors / numpy | 0.8.0 / 2.4.6 |
| 硬件 | 2 vCPU（Xeon @2.60 GHz）、3 GB RAM、无 GPU |
| OS | Linux 6.1.158 x86_64 |

说明：PyPI 当前最新 `torch`（2.14.1）会拉取整套 CUDA toolkit，因此固定到 2.2.2；
该 wheel 仍需 NVIDIA 运行时依赖才能 import，所以完整安装了它的 CUDA 依赖但在 CPU 上运行。
NumPy 2.x 与 torch 2.2.2 会打印一条无害的 `Failed to initialize NumPy` 警告。

## 2. 运行清单

| 运行 | 命令要点 | 退出码 | 事务结果 | 墙钟 |
|---|---|---|---|---|
| 冒烟 | `--base-steps 20 --edit-steps 20` | 1 | `rolled_back` | 7.8 s |
| 基线球 | `train_sphere.py --steps 1800` | 0 | revision 1 发布 | 21.3 s |
| 正式编辑 | `--base-steps 1200 --edit-steps 1000` | 0 | `committed`（revision 2） | 43.3 s |
| 确定性复跑 | 同上，独立进程 | 0 | `committed`，权重哈希逐字节相同 | 43 s |
| 回滚 A | `--depth 0.25 --sigma 0.15` | 1 | `rolled_back` | 46 s |
| 回滚 B | `--depth 0.45 --sigma 0.12` | 1 | `rolled_back` | 46 s |
| 回滚 C | 回滚 B 写入已存在 rev2 的路径 | 1 | `rolled_back`，磁盘 rev2 未被改动 | 46 s |

冒烟测试按预期回滚：它的存在只是证明程序能跑、指标可测量、状态能正确报告，
不代表算法通过，也没有被当作通过。

## 3. 结果

### 3.1 基线球（unseen radius 0.537，1800 步，seed 503）

| 指标 | 实测 | 阈值 | 判定 |
|---|---|---|---|
| `field_mae` | 0.003707 | ≤ 0.008 | PASS |
| `field_rmse` | 0.005329 | — | — |
| 表面径向误差 mean | 0.002184 | ≤ 0.010 | PASS |
| 表面径向误差 P95 | 0.005108 | — | — |
| 表面径向误差 max | 0.010619 | — | — |
| 射线零交点覆盖率 | 1.000 | ≥ 0.95 | PASS |

训练 18.75 s。权重 SHA-256：`3eb7b8bf15161d46772b9eea7852b79e839f58e79711e5fae2bd1b9dec0c7c5d`。

### 3.2 正式编辑（1200 基线步 + 1000 编辑步，seed 503）

| 指标 | 实测 | 验收阈值 | 判定 |
|---|---|---|---|
| `ray_coverage` | 1.000 | ≥ 0.98 | PASS |
| 凹陷目标深度 | 0.060 m | — | — |
| 凹陷区平均径向位移 | 0.04612 m | — | 区域内高斯平均 |
| 凹陷区平均表面误差 | **0.002321 m** | ≤ 0.025 | PASS |
| 保护区平均位移 | **0.001280 m** | ≤ 0.012 | PASS |
| 保护区 P95 位移 | 0.003312 m | — | — |
| 凹陷 / 保护区射线数 | 60 / 3776 | — | — |
| 编辑阶段耗时 | 23.79 s | — | — |

权重 SHA-256：`db1df71e2f681a958b9382904554c21135fedb5af692abf07932e50b57253cd7`。
独立进程复跑得到**完全相同的哈希和全部指标**（确定性成立）。

### 3.3 独立复核（不依赖脚本自身的度量代码）

用种子重放父模型（1200 步），哈希 `62764ef9…`，与运行时的父模型一致；
再用 4096 个固定种子方向、0.20–0.80 m 致密径向扫描 + 二分求零交点：

| 指标 | 实测 |
|---|---|
| 极点（0°）凹陷峰值深度 | **0.05762 m**（目标的 96.0%，误差 0.00076 m） |
| 极点候选表面半径 / 解析目标 | 0.44076 / 0.44000 m |
| 凹陷区（59 条射线）平均 / 最大位移 | 0.04649 / 0.05673 m，全部向内 |
| 保护区（3765 条射线）平均 / 中位 / P95 / P99 / 最大 | 0.001272 / 0.000995 / 0.003368 / 0.005262 / **0.010850** m |
| 保护区最大位移位置 | 175.7°（背侧），向内 |
| 保护区向外漂移比例 | 33.5% |

同时测到一件重要的事：**基线模型自身也不是完美球**。在 R = 0.5 处它的径向形状误差为
mean 0.00409 / P95 0.00973 / **max 0.01715 m**（约 3.4% R）。这比编辑造成的保护区位移
（mean 0.00127）大一个量级——也就是说，候选相对"真实球面"的绝对误差主要来自场表示能力，
而不是这次编辑。评价局部编辑必须用"相对父模型"的位移，脚本正是这么做的，这一点是对的。

### 3.4 回滚用例

| 用例 | 现象 | JSON | 产物 |
|---|---|---|---|
| 冒烟 20/20 | 保护区平均位移 0.0301 > 0.012 | 合法 | 未写目录 |
| depth 0.25 | 覆盖率 0.9922，但 `dent_ray_count = 0`，凹陷指标为 `null` | 合法 | 未写目录 |
| depth 0.45 | 覆盖率 0.9902，`dent_ray_count = 0`，指标 `null` | 合法 | 未写目录 |
| 写入已有 rev2 路径 | 磁盘权重哈希前后均为 `db1df71e…`，未被覆盖 | 合法 | 目录保持 rev2 |

回滚用例 A/B 的目标表面移出了脚本的 `[r-0.14, r+0.08]` 径向搜索窗（极点半径 0.25/0.45
小于 0.36），因此**一条凹陷射线都没找到**，`inf` 被 `json_safe_metrics` 转成 `null`，
`allow_nan=False` 序列化成功、退出码 1、没有发布任何东西。这正是 `7d4fdd0` 要修的路径，
现在有了真实证据：修好了。

## 4. 对照"真正完成"的三个条件

| 条件 | 结论 | 依据 |
|---|---|---|
| 几何有效 | **满足** | 峰值凹陷 0.0576 / 0.06，凹陷区平均误差 0.00232 ≤ 0.025，60 条凹陷射线全部向内 |
| 局部编辑有效 | **有条件满足** | 保护区 mean 0.00127 ≤ 0.012、P95 0.00337；但实测 max 0.01085，已达阈值的 90%，且验收只看 mean |
| 事务有效 | **部分满足** | 成功才发布、失败不写盘、复跑逐字节一致；但 revision 是硬编码的，发布不是原子操作，父版本从未落盘 |

## 5. 这次运行暴露的问题

### P1 — revision 是硬编码的，事务只支持 rev1 → rev2

`examples/neural_geometry/edit_sphere_dent.py:96,103,112` 把 `revision: 2` / `parent_revision: 1`
写死；回滚结果 `:238` 也固定输出 `"active_revision": 1`。

证据：回滚用例 C 写入一个**已存在 revision 2** 的路径，脚本仍报告
`active_revision: 1`，而磁盘上 `load_bundle()` 读到的是 revision 2。也就是说
"失败后活动版本仍是 1"这句话在真实文件系统上不成立——它只是常量。

### P1 — 发布不是原子的

`build_bundle()` 直接 `mkdir` 最终目录并逐个写 `weights/manifest/state/training/validation`。
没有临时目录 + `os.replace`，也没有"最后写 manifest"的提交点。本次没有崩溃所以没事，
但一次磁盘写失败或进程中断会留下半成品 revision，而 `load_bundle` 只校验哈希与路径，
可能读到"权重已换、manifest 未换"的状态。这距离"事务"还差一步。

### P2 — 凹陷深度指标的语义会误导

`dent_depth_measured_mean_m = 0.0461` 是**高斯凹陷区域内所有射线的平均位移**，而
`dent_depth_target_m = 0.06` 是**峰值深度**。两者并列会让人以为只做到了目标的 77%，
而独立复核的极点峰值是 0.0576（96%）。建议改名为
`dent_depth_region_mean_m`，并新增 `dent_depth_peak_m`。

### P2 — 容差键名与实际门禁不一致

`validation.json` 声明 `dent_depth_error_m: 0.025`，但代码用它约束的是
`dent_surface_error_mean_m`（`:222-226`）。键名应改为 `dent_surface_error_mean_m`，
否则报告读者会以为门禁卡的是深度而不是表面误差。

### P2 — 保护区只卡均值

3765 条射线的 max 位移 0.01085 m 已接近 0.012 m 的均值阈值；P99 为 0.00526。
少量射线可以漂到均值的 8.5 倍而不触发失败。建议门禁加入 P99 与 max，或至少把它们写进
`validation.json`。

### P3 — 失败原因不可诊断

`"reason": "validation thresholds not met"` 对"没有一条有效凹陷射线"和"误差超限"是同一句话。
应输出 `failed_checks: ["dent_rays_missing"]` / `["dent_surface_error"]`，并保留
`dent_ray_count`（已有）与门限值。

### P3 — CI 从不执行回滚分支

`neural-geometry-experiments.yml` 只跑两个应当成功的实验，因此 `7d4fdd0` 修的那条路径
在 CI 里从未被走到。建议增加一个 `continue-on-error` 的失败注入 job，断言：
退出码 1、stdout 是合法 JSON、没有产生新 revision。

### P3 — 父版本没有落盘

编辑实验的父模型只存在于内存里，输出目录只有 revision 2。本次通过"用种子重放父模型"
独立验证了父哈希（`62764ef9…`），但这依赖于训练过程逐位可复现；更稳的做法是基线训练
把 revision 1 写到磁盘，编辑时用 `--parent-bundle` 读取并在回滚后用哈希对比。

## 6. CI 状态

| 运行 | 工作流 | 结论 | 用时 |
|---|---|---|---|
| [40063938007](https://github.com/cxs503/TensorCAD/actions/runs/40063938007) | neural-geometry-experiments (main) | success | 2m16s |
| [40063937945](https://github.com/cxs503/TensorCAD/actions/runs/40063937945) | tests (main) | success | 35s |

更早的 main 运行 `40063511445`（PR #10 合并后）神经几何任务失败，由 PR #11 修复依赖与
采样后转绿。PR #8 分支上的 `tests` 最初也失败过两次（`40059592624`、`40059601575`），
在 `research/neural-geometry-online-edit` 上修好后转为成功；但当时还没有
`neural-geometry-experiments` 这个 workflow——该文件直到 04:01 才由 `925e110` 加入。
CI 的通过只能说明两个实验达到了各自脚本内的阈值；它不保存指标，也不执行
回滚用例——本次补上了这两件事。

## 7. 建议的下一步（最小改动，按顺序）

1. `build_bundle` 接受 `revision` / `parent_revision` 参数，从现有 manifest 读取并 +1；
   回滚时从磁盘读回 `active_revision`，读不到就输出 `null`，不要写常量。
2. 原子发布：写 `output.tmp-<uuid>/`，`fsync` 后 `os.replace`；manifest 最后写。
3. 指标改名 + 新增 `dent_depth_peak_m`、`protected_p99_displacement_m`、
   `protected_max_displacement_m`，容差键与指标同名。
4. 失败结果增加 `failed_checks` 列表。
5. CI 增加失败注入 job（`--depth 0.45`），断言退出码、JSON 合法、目录未创建。
6. 基线训练落盘 revision 1，编辑用 `--parent-bundle` 读入，回滚后比较哈希。
7. 以上稳定后，再把 `apply_edit → train → validate → commit/rollback` 抽象成通用接口，
   然后才扩到椭球/圆柱/船体。

## 8. 产物

```
README.md                       本报告
metrics.json                    全部运行、指标、哈希、独立复核的机器可读汇总
environment.json                Python / torch / 硬件
commands.md                     逐条复现命令（含 torch 安装注意事项）
logs/*.stdout.json              每次运行的原始 stdout（含回滚用例）
logs/dent-1200-1000.stderr.log  正式运行的 stderr
verification/verify_profile.py  角向剖面独立复核（重放父模型 + 二分求根）
verification/protected_stats.py 保护区位移分布独立复核
verification/*.json             上述脚本的输出
bundles/neural-sphere.ngo/          revision 1：manifest/state/training/validation/weights
bundles/sphere-dent-seed-503.ngo/   revision 2：同上
```

权重均 < 220 KB，`load_bundle()` 可校验 SHA-256；这些文件是本次真实计算的证据，
不是手工编辑的摘要。
