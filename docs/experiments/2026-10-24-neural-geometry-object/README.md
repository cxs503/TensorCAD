# 神经几何对象：真实运行记录（0.3.0）

本目录是本机真实执行 `examples/neural_geometry/run_object_pipeline.py` 的结果。
它验证的不是"脚本能跑完"，而是四个可判定的性质：字段可保存/重载、编辑是事务、
派生网格可用、失败编辑不发布任何东西。

```bash
python examples/neural_geometry/run_object_pipeline.py --output results/object
# 1m33s on 2 vCPU, CPU only
```

环境：CPython 3.11.2，torch 2.2.2（CPU），safetensors 0.8.0，numpy 2.4.6，无 GPU。

## 1. 拟合（revision 1）

| 指标 | 实测 | 门限 |
|---|---|---|
| field MAE | 0.005092 m | ≤ 0.008 |
| field RMSE | 0.007775 m | — |
| 表面径向误差 mean | 0.004086 m | ≤ 0.010 |
| 表面径向误差 P95 / max | 0.009899 / 0.015768 m | — |
| 射线覆盖率 | 1.000 | ≥ 0.95 |
| 训练耗时 | 12.41 s（1200 步） | — |
| 权重 SHA-256 | `62764ef9…3181` | 与独立重放一致 |

## 2. 编辑事务（revision 2）

| 指标 | 实测 | 门限 |
|---|---|---|
| 凹陷峰值（目标 0.060 m） | **0.052690 m** | 误差 ≤ 0.025 → 实测 0.007310 |
| 凹陷区平均表面误差 | 0.002825 m | ≤ 0.020 |
| 保护区平均位移 | 0.001209 m（有符号 −0.000773） | ≤ 0.012 |
| 保护区 P95 / P99 / max | 0.003037 / 0.004882 / **0.009653 m** | ≤ 0.020 / ≤ 0.030 |
| 射线覆盖率 | 1.000（4096 随机 + 64 定点锥内） | ≥ 0.95 |
| 编辑耗时 | 24.87 s（1000 步） | — |

凹陷峰值达到目标的 87.8%；保护区最大漂移 9.7 mm，是均值的 8 倍，仍低于 30 mm 门限——
这正是把 P99/max 也写进门禁的原因。

## 3. 派生网格（96³ 采样）

```
vertices 97332   triangles 194660
watertight true   edge-manifold true   euler characteristic 2   signed volume > 0
volume 0.515400 m^3   area 3.175588 m^2
obj 10006315 bytes   stl 9733084 bytes
```

球体本身解析体积 0.523599 m³；凹陷后 0.5154 m³，符合"向内压 6 cm"的预期。
两者都以哈希记录在 `bundle/manifest.json` 的 `derived[]` 里，因此任何一份导出的
STL 都能追溯回它的 revision 与分辨率。为控制仓库体积，本目录不包含 20 MB 的网格文件，
按上面的命令可逐位重放。

## 4. 失败注入

用 `--reject-depth 0.25 --sigma 0.15`（凹陷比搜索窗还深）触发真实失败：

```
status rolled_back
failed checks: dent_peak_error, dent_region_error_mean
active revision 2（未变）
weights / manifest 字节不变
history 追加 rev3 rolled_back，并保留 metrics 与 failed_checks
```

`history.json` 里的 `inf`/`nan` 已写成 `null`，报告本身始终是标准 JSON。

## 5. 这次运行同时固定的行为

* `fit` 不会覆盖已有对象（否则会静默丢掉一条血缘）。
* `expected_parent` 与磁盘版本不符时提交失败，两个写入者不会互相覆盖。
* 失败候选从不进入活动目录；成功候选先写 staging，校验哈希后再做目录级 rename。
* 父版本权重哈希记录在 `training.parent_weights_sha256`，可与 `history[0]` 对账。
* `tensorcad-neural verify` 能检出被篡改的权重（哈希不符）并返回 1，而不是抛异常。

## 6. 已知限制（实测）

* 网格质量：最小三角形面积约 3×10⁻¹⁶ m²，属于 marching tetrahedra 在对称轴附近的
  sliver 三角形。它们不影响水密性与体积（占总面积可忽略），但对要求高质量面元的
  求解器应加一步重网格化。
* 场精度 ~3% R：球面最大径向误差 15.8 mm，大于本次编辑造成的保护区漂移（mean 1.2 mm）。
  提升精度要靠更大的网络或更长的训练，不是格式问题。
* 提取耗时：96³（97 万采样点）约 8 s；48³ 约 1 s。CPU 上 160³ 需数十秒。

## 7. 产物

```
summary.json          全部指标、历史与检查结果
run.log               完整 stdout（含每一条 PASS）
bundle/manifest.json  revision 2 的完整元数据（含 derived 清单）
bundle/state.json     可编辑控制项（radius=0.5, dent depth=0.06, sigma=0.22）
bundle/training.json  这一版的训练记录与父权重哈希
bundle/validation.json 逐项检查与阈值
bundle/history.json   三次尝试：fit / dent_edit committed / dent_edit rolled_back
bundle/weights.safetensors  218916 字节，SHA-256 由 manifest 校验
```
