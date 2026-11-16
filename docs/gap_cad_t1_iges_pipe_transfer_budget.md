# GAP-CAD-T1-IGES-PIPE：trim-on-surface 验收预算（两档预算+诚实标注）

**单号** m22709 ｜ **仓** TensorCAD ｜ **分支** agent/m22709 ｜ **证据目录** `artifacts/gap_t1_transfer_budget/`
**改动面** `tensorcad/geomquality/trim_mesh.py`（`tessellate_trimmed_model` / `import_iges` / `_gates` / `Tessellation`）+ `tests/test_iges_transfer_budget.py`

---

## 1. 案例介绍（物理背景 / 几何 / 参数）

金标样本是被审过两轮的第三方工业导出：**KP505 桨（模型缩尺 D=250 mm）**，
`benchmarks/anchors/sources/KP505_modelScale_250mm.igs`，1 042 550 B，sha256 `0bb54ae7…462f`，
单位 MM，D 段含 144/143/142（trimmed B-rep）——即锚点语义要求的裁剪面路线，不是自写 patch-set。

| 参数 | 值 |
|---|---|
| 面数（文件） | 8（4× B-spline 叶片面 + 毂回转面 + 2× 平面端面 + 1× 平面） |
| trim 表示 | 5 面走 142 参数空间像（parametric），3 面走 3-D 曲线反演（inverted） |
| 模型 extent（最大坐标绝对值） | 114.84 mm |
| 弦高预算 chord | 1e-3·extent = 0.11484 mm |
| 旧验收（exact） | 1e-6·extent = 1.1484e-4 mm |
| 新验收（transfer） | 1e-4·extent = 1.1484e-2 mm |

## 2. 怎么做的（方法 / 流程 / 关键设置）

1. **复现缺陷**：旧代码对 KP505 只出 7/8 面——`refused: {1: "inversion residual 2.213e-04 > tolerance 1.148e-04"}`。
   逐点探针（`invert_point` 沿环走 + 精确切向导数栅格回退）显示毂面 trim 环上第 55–60 号点在
   u≈6.39 处出现**平滑鼓包**（4e-5 → 2.2e-4 → 4e-5），不是求解器卡死（卡死会抖动、且密集栅格搜不到更优点），
   而是导出文件里那条 3-D 曲线与其回转面之间 ~0.22 µm 的转移偏差（3-D 曲线是曲面交线的样条近似）。
2. **定位判据错配**：网格边界顶点本来就由反演得到的 `(u, v)` 发点（`refine_parametric` →
   `S(u,v)`），也就是说**残差根本不进网格**：被拒的 1/8 模型、~210 条边界边，是拿一个网格精度
   1000 倍严的尺子量出来的。0.22 µm / 250 mm = 1.5e-6 相对量。
3. **改为两档预算 + 标注**（最小改动，无几何搬运）：
   * `residual_tolerance`（exact，1e-6·extent）保留为"干净反演"的尺子；
   * 新增 `transfer_tolerance`（1e-4·extent，即 chord 预算的 1/10）为**准入**预算；
   * 残差落在两档之间的面 → 照常三角化，逐面记 `residual_class="transfer_budget"`，
     `Tessellation.faces_in_transfer_band` / 报告 `import.faces_in_transfer_band` 指名道姓；
     G2 的 `limit` 改为 transfer，并附 `exact_limit` 与 `faces_in_transfer_band` 计数，
     **绝不把"放宽准入"写成"精确反演"**；
   * 超出 transfer 仍拒（错曲面/错变换/错配对都是 1e0 量级，比预算高 4 个数量级）；
   * 顺带修一处**空洞通过**：`max_inversion_residual` 现在统计含被拒面的全部实测残差，
     旧行为下"全部面被拒"会报 0.0 让 G2 假绿。
4. **工具链**：纯 stdlib 管线（`python -m tensorcad.geomquality.import_cli`）+ pytest；
   重测试（含 3 个 IGES 集成测试文件）递交集群 `cj4322`（8 核），本机只跑 <120 s 子集。

## 3. 计算结果（程序化注入，来自 `artifacts/gap_t1_transfer_budget/*.json`）

| 量 | 旧门（7/8） | 新门（8/8） |
|---|---|---|
| faces tessellated / 文件 | **7/8** | **8/8** |
| faces refused | 1（毂面，残差 2.213e-4 mm） | 0 |
| 最大 trim 反演残差 | 2.983e-5（被拒面不计入，旧口径失真） | 2.213e-4（含 transfer 档标注） |
| 三角数（repair 后） | 69 468 | 77 429 |
| free edges（repair 前→后） | 468 → 159 | 678 → **80** |
| non-manifold edges（后） | 38 | 462 |
| duplicate triangles（后） | 11 | 11 |
| components / closed | 28/2 | 126/2 |
| signed volume（模型单位³） | 5.9334e-5 | 8.3061e-5 |
| 表面积 | 0.24595 | 0.24925 |
| 门 | G1 FAIL, G2 PASS, G3/G4/G6/G7 FAIL, G5 PASS | G1 **PASS**, G2 **PASS**, G3/G4/G6/G7 FAIL, G5 PASS |

图：`artifacts/gap_t1_transfer_budget/figure.png`
（左=毂面残差 vs exact/transfer/chord/粗错配四档预算，log 轴；右=旧门/新门拓扑量对比，symlog 轴）。

## 4. 结果判定

* **达成**：锚点验收里"CLI+库函数"与"真实 IGES 全过门"的第一道硬门 **G1 8/8 面读入**由 7/8 修到 8/8，
  且 G2（trim 在曲面上）以新口径 PASS 并自带证据；这项修复不搬动任何几何（顶点仍由 `(u,v)` 生成，
  合成用例里网格顶点到曲面偏差实测 **0.0**）。
* **未达成（honest FAIL，不粉饰）**：G3 水密（80 条自由边）、G4 流形（462 条非流形边）、
  G6 单壳（126 分量）、G7（11 个重复三角）仍红。**并新增一条负面发现**：修好 8/8 后 repair 链
  （weld→hole-fill）制造的非流形边从 38 涨到 462（repair 前是 0），即 repair 在缝合 T 型接缝时
  把"自由边"换成了"非流形边"——这一步需要真正的拓扑缝合（`stitch_t_junctions` 只在显式传
  `seam_tolerance` 时启用，实测 1e-3 时 free edges 121 / non-manifold 509 / 分量 163，仍未闭合）。
  这是下一步的选题，不在本单范围。

## 5. 注意点（坑与局限）

* **放宽判据必须留痕**：`faces_in_transfer_band` 是这条改动的兑现物；只看 `G2 PASS` 会漏掉
  "这一面是带账通过的"。
* **残差与网格无关**：边界顶点从 `(u,v)` 发点，所以容差放宽不改变网格位置；这一点由合成用例
  `_on_surface_deviation ≤ 1e-12` 显式证明，也是本改动不构成"以放过换绿灯"的依据。
* 预算默认值绑定 extent（模型的相对量），不是绝对 mm；换船型/换缩尺自动随之缩放。
* 桨是 patch 装配体，第三方导出普遍带 T 型接缝，`tests/test_gap_cad_t1_iges_pipe.py` 的模块
  docstring 早已写明"全部闭门不要求"；本单没有为了绿灯去动 `weld_tolerance` 或关掉 gate。
* 本机 load15 ≈ 72/128，重测试递交集群；本机只跑 <120 s 用例，未起新的重算进程。

## 6. 误差判定（逐项，显式回答"是否 <3%"）

| 关键量 | 基准 | 实测 | 误差 | <3% ? |
|---|---|---|---|---|
| 面读入率（KP505，8 面文件） | 8/8 | 8/8 | 0 % | **是** |
| 毂面 trim 残差 | transfer 预算 1.1484e-2 mm | 2.213e-4 mm | 内缩 52 倍（1.9e-2 预算占比） | **是** |
| 残差 / 模型 extent | 0（理想） | 1.93e-6 相对 | — | **是**（远优于 3%） |
| 网格顶点到曲面偏差（合成用例） | 0 | 0.0（≤1e-15） | 0 % | **是** |
| 仍被拒的粗错配（合成 1e-2 相对） | 应拒 | 已拒（报 "transfer budget"） | 判据双侧均被测试钉住 | **是** |
| 水密 / 流形 / 单壳 | 0 自由边 / 0 非流形 / 1 壳 | 80 / 462 / 126 | 功能未达成 | **否（honest FAIL）** |
| 桨体积/面积外部基准 | 无公开体积基准（KP505 只公开 D 与敞水性能） | — | — | **不适用（无基准 ≠ 通过）** |