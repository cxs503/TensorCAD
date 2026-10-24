# 神经几何对象格式（NGOM 0.2 / tensorcad 0.3）

本文说明 TensorCAD 如何把一个几何体存成"可编辑的连续场 + 元数据"，以及为什么它比
STL/OBJ/IGES 更适合作为**权威**几何表示。文档描述的是已经实现并跑通的格式，不是设计草案。

## 1. 问题：网格文件不是几何模型

STL/OBJ 是三角面片快照，IGES/STEP 是 B-rep 快照。作为**交换**格式它们没问题，作为
**保存**格式它们丢的东西恰好是工程上最需要的：

| 丢失的信息 | 后果 |
|---|---|
| 生成参数 | 想改一个半径只能重新建模，无法回放"这个形状是怎么来的" |
| 单位与坐标系 | 读进来只能假设，米制/英尺制搞错不会报错 |
| 拓扑来源 | 面片密度是一次性决定，改分辨率等于重做 |
| 版本与操作历史 | 无法回答"这版和上版差在哪、为什么被接受" |
| 有效性证据 | 是否水密、是否自交，只能靠下游再去查 |
| 编辑语义 | "在这个位置压一个 6 cm 的凹坑"无法表达，只能重新导出一份新网格 |

神经几何对象把几何存成 `s = F_θ(x, c)`：`θ` 是网络权重，`c` 是可编辑状态，
零等值面就是表面。**分辨率不再是存储属性**——要 96³ 还是 512³ 的网格，是导出时的一次计算，
而不是保存时的一次决定。

## 2. 对象包结构

```
part.ngo/
  manifest.json          # 入口：架构、场语义、状态、血缘、派生件清单
  weights.safetensors    # 只含张量，不反序列化任何可执行对象
  state.json             # 可编辑控制项与当前取值（field_state）
  training.json          # 这一版是怎么训练出来的
  validation.json        # 实测指标、阈值、逐项判定
  history.json           # 每次提交/回滚尝试，最旧在前
  derived/               # 可选：surface.obj / surface.stl 等派生产物
```

`manifest.json` 的完整 JSON Schema 在 `schemas/neural-geometry-manifest-v0.2.schema.json`。
读取端只依赖标准库：`tensorcad.neural_geometry` 校验元数据、哈希与路径，从不执行包内代码。

## 3. manifest 0.2 的字段

```jsonc
{
  "format": "tensorcad.neural-geometry",
  "format_version": "0.2.0",
  "object_id": "urn:tensorcad:geometry:sphere-r0.5",
  "revision": 2,

  "lineage": {                       // 这一版从哪来
    "parent_revision": 1,
    "root_revision": 1,
    "operation": "dent_edit",
    "history_file": "history.json"
  },

  "model": {
    "family": "conditional_implicit_field",
    "architecture_id": "tensorcad.mlp_fourier.v1",
    "architecture_config": {
      "coordinate_dim": 3, "state_dim": 1, "output_dim": 1,
      "hidden_width": 128, "hidden_layers": 4,
      "coordinate_encoding": "fourier", "fourier_bands": 6
    },
    "weights": { "path": "weights.safetensors", "sha256": "…", "bytes": 218916 }
  },

  "field": {                         // 场是什么，怎么读
    "quantity": "truncated_signed_distance",
    "inside_sign": "negative",
    "surface_isovalue": 0.0,
    "length_unit": "m",
    "truncation_distance": 0.15,
    "continuity": "continuous",
    "resolution": "implicit; bounded only by the query point density"
  },

  "coordinates": {                   // 单位与坐标系不是注释，是数据
    "dimension": 3,
    "frame": "object_local",
    "units": "m",
    "bounds": { "min": [-0.58, -0.58, -0.58], "max": [0.58, 0.58, 0.58] }
  },

  "state": {                         // 可编辑控制项
    "field_state": { "radius": 0.5 },
    "controls": {
      "radius": { "value": 0.5, "unit": "m", "role": "field_state", "edit_mode": "conditioned" },
      "dent":   { "value": 0.06, "unit": "m", "role": "online_edit",
                  "edit_mode": "online_fit", "sigma_radians": 0.22, "direction": [1, 0, 0] }
    }
  },

  "training":  { "method": "online_fit_with_replay", "parent_weights_sha256": "…", "steps": 1000, "seed": 503 },
  "validation": { "status": "passed", "metrics": { … }, "tolerances": { … }, "failed_checks": [] },
  "provenance": { "producer": "TensorCAD", "created_utc": "2026-10-24T06:12:00+00:00", "command": [ … ] },

  "derived": [                       // 派生件：删掉也不损失任何信息
    { "kind": "surface_mesh", "format": "obj", "path": "derived/surface.obj",
      "sha256": "…", "bytes": 7231234,
      "metrics": { "watertight": true, "euler_characteristic": 2, "triangle_count": 194660 } }
  ],

  "state_file": "state.json", "training_file": "training.json", "validation_file": "validation.json"
}
```

设计约束（都有测试或代码强制）：

* **单位、坐标系、截断距离、符号约定**必须显式，不允许默认值推断。
* **未知字段保留而不丢弃**：`additionalProperties: true`，读端不会因为版本更新而丢元数据。
* **路径必须是包内相对路径**，`../` 与绝对路径一律拒绝。
* **权重必须有 SHA-256（可选字节数）**，不匹配即报错。
* **架构白名单**：只有声明过的架构可加载，防止"几何文件里藏一段 Python"。
* **0.1 清单仍可解析**：新增字段全部可选，旧包不会因为升级而打不开。

## 4. 事务模型

一次编辑不是"覆盖写文件"，而是一次事务：

```
apply_edit(request)
   ├── 读取磁盘上的活动版本（父版本 = revision N，父权重哈希）
   ├── 训练候选（replay 保护区域）
   ├── 度量（峰值深度、凹陷区误差、保护区 mean/P95/P99/max、射线覆盖率）
   ├── 逐项门禁 → 任一失败
   │      ├── 候选写入临时目录（从不出现为活动版本）
   │      ├── 删除临时目录
   │      ├── history.json 追加 {status: rolled_back, failed_checks, metrics}
   │      └── 活动版本字节不变
   └── 全部通过
          ├── 在 staging 目录内写好权重与全部元数据（含继承的历史）
          ├── 校验 staging 内的权重哈希与清单一致
          ├── expected_parent 与磁盘 revision 不一致 → 拒绝（乐观并发）
          └── 目录级 rename 原子发布；history.json 追加 {status: committed}
```

不变量（`tests/test_neural_store.py` 覆盖）：

1. 未完成的编辑不可见——候选只在隐藏 staging 目录里；
2. 校验失败时活动版本**逐字节不变**；
3. 每次尝试（成功与失败）都进 `history.json`，失败项带 `failed_checks` 与实测指标；
4. `expected_parent` 不等于磁盘版本时提交失败，两个人不会互相覆盖；
5. 非有限测量值（`inf`/`nan`）序列化为 `null`，失败报告本身永远可解析。

## 5. 派生件：STL/OBJ 的正确位置

`tensorcad-neural extract` 从活动版本派生网格：

```bash
tensorcad-neural extract --object results/sphere.ngo --resolution 96 --format obj,stl
```

实测（96³ 采样、marching tetrahedra）：

```
vertices 97332   triangles 194660
watertight true  edge-manifold true  euler characteristic 2
signed volume 0.515400 m^3   surface area 3.175588 m^2
```

派生件记录在 `manifest.derived` 里，带哈希与网格质量指标，因此任何人都能确认
"这个 STL 是哪个 revision、哪个分辨率导出来的"。删掉它们不损失信息，重新导出即可。

提取实现是 `tensorcad.neural.extract`：Kuhn 六四面体细分 + 逐面局部定向，
不使用 marching cubes 的 case 表，因此不依赖容差焊接也不会开裂；
在 48³ 下球体体积误差 −0.13%，96³ 下 −0.03%。

## 6. 与 STL/IGES 的关系

| | 神经几何对象 | STL/OBJ | STEP/IGES (B-rep) |
|---|---|---|---|
| 权威表示 | 场 + 状态 + 元数据 | 三角面片 | NURBS 面 |
| 分辨率 | 查询/导出时决定 | 固定 | 解析 |
| 可编辑参数 | 是（状态、编辑事务） | 否 | 是（特征树，但需 CAD 内核） |
| 版本与验收证据 | 内建（history + validation） | 无 | 取决于 PDM |
| 求解器直接可用 | 需先导出网格 | 是 | 需先离散 |
| 依赖 | 推理时需 torch，元数据层仅标准库 | 无 | 需 CAD 内核 |

因此本项目的定位不是"立刻取代 STEP"，而是：**把神经场作为可编辑、可追溯的权威表示，
把 STL/OBJ 作为它的导出**。B-rep 的精确曲面（倒角、螺纹、孔轴配合）目前仍只有 CAD 内核能表达，
`tensorcad.adapters.ocp` 是后续接缝（issue #4）。

## 7. 当前边界

* 架构只有一个：`tensorcad.mlp_fourier.v1`（坐标 Fourier 特征 + 状态向量 + MLP）。
* 已实现的形状族是球族（半径可条件化，可做局部凹陷）。椭球/圆柱/船体尚未接入，
  但格式与事务层与形状无关，扩展只需新增 `fits_*` 与相应 oracle。
* 场精度受表示能力限制：0.5 m 半径的球，径向误差 mean ≈ 0.004 m、max ≈ 0.017 m（≈3.4% R），
  比同一模型的编辑漂移大一个量级——瓶颈在网络规模和训练预算，不在格式。
* 推论时延：CPU 上 96³ 提取（约 97 万个采样点）≈ 8 s（本项目环境，2 vCPU）。
* 不是精确几何：对需要严格 G1/G2 连续或公差链的零件，仍应导出到 B-rep。
