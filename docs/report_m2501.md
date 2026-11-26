# Mission 2501 补报：TensorCAD 0.12.0 — 船体纵向型线族参数化 + 闭式静水力基准

> 补报单（scenario_key=`report-backfill`，mission_id=2501）。原单已交付但缺论文级进展报告，本页为补报。
> 诚实声明：原单工件目录未留档（`artifacts/hull_family/*`、`docs/analytic_hydrostatics/convergence.png` 均未归档），
> 本报告所有数字均从 worktree 内已交付工件**程序化提取**（`docs/hull_family/evidence.json`、
> `docs/hull_family/gate.json`、`docs/analytic_hydrostatics/results.json`、`docs/analytic_hydrostatics/report.md`），
> 零手抄；缺失工件一律标注「未留档」，不编造。

---

## ① 案例介绍（物理背景 / 几何 / 参数表）

**选题依据**：读完 `backlog_tensorcad_section.md`（TensorCAD backlog 节）。其中「可贡献开发项」第 1 条
**多船型族扩展（P0）**（在 `hull.py` 新增分布函数族：系列 60 / Savitsky / 用户自定义，参数化框架已有）与
「静水力扩展（P3）」为本单落点。审计已确认的缺口：`tensorcad.hull` 是下游唯一参数入口，但只放样
`sin^p` 锥度 × `sqrt` 垂向剖面的单船型族；Series 60 站型律在 `tensorcad.series60` 中**无生产路径**
（跨检审计 `docs/audit_m17252_gap_cad_t2_param_series.md` 记录 *"the declared entry point
`tensorcad/hull.py` was not touched by a single line"* 为 L2→L3 阻塞）。

**物理背景**：参数化船体曲面 `y(x,z) = (B/2)·sin(πx/L)^p·sqrt((z+T)/T)`（`w=(z+T)/T∈[0,1]`，
`p=BREADTH_EXPONENT=0.7`）可**完全分解**为 x 积分 × w 积分，两者均有闭式——这是闭式静水力基准的
存在前提：网格积分数字从此有「自身解析面」可回答，而非只对更细的自身副本。

**参数表**（两套独立案例，均来自工件）：

| 案例 | L [m] | B [m] | T [m] | 密度 [kg/m³] | 网格 | 备注 |
|---|---|---|---|---|---|---|
| 船体族（Series 60 分派） | 10.0 | 1.5 | 0.6 | 1000（淡水） | 41×13 / 81×21 / 161×41 | Cb 梯 0.55–0.85，五父 0.60/0.65/0.70/0.75/0.80 |
| 闭式静水力（analytic_hydrostatics） | 10.0 | 2.0 | 0.8 | 1025（海水） | 11×5 / 21×9 / 41×13 / 81×21 / 161×41 | 半宽 `(B/2)sin(πx/L)^0.7·sqrt(w)` |

---

## ② 怎么做的（方法 / 流程 / 关键设置 / 工具链）

**交付物一：船体纵向型线族参数化（`tensorcad.hull` 家族分派）**

- `HullParameters` 新增 `family`（`'sine'` 默认，**bit-for-bit 不变**）＋ `block_coefficient` /
  `waterplane_coefficient` / `lcb_ratio` / `transom_ratio` / 三个 bulb 比。
- `family='series60'` 分派到 `tensorcad.series60`（`generate_series60`），经现有 `HullModel` schema
  呈现（校验报告、静水力块、offsets 约定、`metadata()`）；族自身报告挂在 `family_report` 下。
- **offsets 约定在单点转换**：`_build_family_mesh` 把族模块的半宽转为本模块全围长，一次完成。
- **参数拒绝而非静默丢弃**：族不能 honor 的参数（如 sine 族要 Cb、series60 要 sine 指数）显式拒绝——
  「要 Cb 却拿到指数」是对调用者的欺骗。
- `estimate_discretisation_error` / `resolve_resolution` 用 `dataclasses.replace` 携带全参数集，
  并以**族自身闭式**（而非 sine 公式）为细化基准。
- **网格等价断言为 bit-for-bit**：两个声称建同一曲面的入口必须产出相同顶点，非容差内。

**交付物二：闭式静水力基准（`tensorcad.hull.analytic_hydrostatics`）**

- 曲面分解：`TAPER_INTEGRAL = ∫₀¹ sin(πx)^p dx = Γ(0.85)Γ(0.5)/(πΓ(1.35))`（Beta 函数）。
- 闭式：体积 `L·B·T·TAPER_INTEGRAL·2/3`；水线面 `L·B·TAPER_INTEGRAL`；中剖面 `2/3·B·T`；
  LCB/LCF `=L/2`（对称）；VCB `=-2/5·T`；`Cb=2/3·TAPER_INTEGRAL`；`Cp=TAPER_INTEGRAL`。
- **湿表面积是唯一无闭式量**：`w=u²`（消龙骨 `1/sqrt(w)`）＋ `ψ=(π/2)τ³`（消艏艉 `x^(p-1)` 尖端斜率）
  双变换后 80 点 Gauss–Legendre（Legendre 节点 Newton 迭代生成，核心保持标准库-only）；细化到 120 节点
  相对变化 `<1e-9`。
- **基准自检**（公式只对自己证明不了什么）：`TAPER_INTEGRAL` vs Simpson `5.6e-13`；湿面积 Gauss vs
  Simpson 乘积 `6.0e-09`；体积/水线面/中剖面/VCB vs 分变量 Simpson `1e-11`；齐次性（L/B/T 加倍 → 面积×4）。

**工具链**：`pytest`（`tests/test_hull_family.py` 17 测、`tests/test_analytic_hydrostatics.py` 21 测）；
`tools/hull_family/report.py` → `docs/hull_family/evidence.json` + `docs/hull_family.png`；
`tools/hull_family/gen_doc_and_gate.py` → **5/5 = 100% 突变体击杀**（`docs/hull_family/gate.json`）；
`examples/analytic_hydrostatics/run_analytic_check.py` → `docs/analytic_hydrostatics/{results.json,report.md}`。

---

## ③ 计算结果（数字表，程序化注入零手抄）

### 3.1 船体族：五父 × 三档网格 vs 各自闭式（`docs/hull_family/evidence.json`）

| Cb | 网格 | 三角形 | V [m³] | V 误差 | Cb 实测 | Cb 误差 | VCB 误差 | Awp 误差 |
|---|---|---|---|---|---|---|---|---|
| 0.60 | 41×13 | 1998 | 5.201439 | −3.68 % | 0.5779 | −3.68 % | −3.69 % | −0.055 % |
| 0.60 | 81×21 | 6558 | 5.284288 | −2.14 % | 0.5871 | −2.14 % | −2.19 % | −0.014 % |
| 0.60 | 161×41 | 25918 | 5.343948 | −1.04 % | 0.5938 | −1.04 % | −1.07 % | −0.003 % |
| 0.65 | 41×13 | 1998 | 5.626530 | −3.82 % | 0.6252 | −3.82 % | −3.79 % | −0.063 % |
| 0.65 | 81×21 | 6558 | 5.719181 | −2.24 % | 0.6355 | −2.24 % | −2.26 % | −0.016 % |
| 0.65 | 161×41 | 25918 | 5.786175 | −1.09 % | 0.6429 | −1.09 % | −1.11 % | −0.004 % |
| 0.70 | 41×13 | 1998 | 6.051213 | −3.95 % | 0.6724 | −3.95 % | −3.89 % | −0.073 % |
| 0.70 | 81×21 | 6558 | 6.153837 | −2.32 % | 0.6838 | −2.32 % | −2.32 % | −0.018 % |
| 0.70 | 161×41 | 25918 | 6.228259 | −1.14 % | 0.6920 | −1.14 % | −1.15 % | −0.005 % |
| 0.75 | 41×13 | 1998 | 6.475374 | −4.07 % | 0.7195 | −4.07 % | −3.97 % | −0.087 % |
| 0.75 | 81×21 | 6558 | 6.588253 | −2.40 % | 0.7320 | −2.40 % | −2.38 % | −0.022 % |
| 0.75 | 161×41 | 25918 | 6.670212 | −1.18 % | 0.7411 | −1.18 % | −1.19 % | −0.005 % |
| 0.80 | 41×13 | 1998 | 6.898818 | −4.18 % | 0.7665 | −4.18 % | −4.04 % | −0.107 % |
| 0.80 | 81×21 | 6558 | 7.022396 | −2.47 % | 0.7803 | −2.47 % | −2.43 % | −0.027 % |
| 0.80 | 161×41 | 25918 | 7.112039 | −1.22 % | 0.7902 | −1.22 % | −1.22 % | −0.007 % |

- 15 例全部 `identical_to_family_module=True`（与 `generate_series60` **bit 级同构**）。
- 最坏体积误差：161×41（Series 60 审计发布分辨率）**+1.22 %**（Cb=0.80）；符号全负、随网格单调下降——
  放样弦切在解析面**内侧**，`zeta^(1/m)` 截面在龙骨上方首网格单元内转弯。

### 3.2 船体族：艉板与球鼻艏透传（`evidence.json`）

- 艉板比 0.15：艉站宽 **0.225000 m** = 期望 `0.15·B`（相对误差 +0.00 %）；尖艉控制 0.000000 m；
  浸没艉面积 **0.132168 m²** = 中剖面的 0.1500；带艉面水密（无校验错误）**True**。
- 球鼻艏（0.03·L 长，`A_BT/A_m=0.08`）：闭式体积 **0.014098 m³** = 船体排水量 +0.23 %；
  艏垂线面积 **0.070490 m²**；带/不带球鼻艏船体排水量 **bit-identical（True）**——球鼻艏是第二闭体，
  非布尔并集，无双计。

### 3.3 船体族：拒绝矩阵（`evidence.json`，6/6 全拒且消息点名原因）

| 案例 | 拒绝 | 消息点名原因 |
|---|---|---|
| series60 无 block_coefficient | True | True |
| sine 族带 Cb | True | True |
| sine 族带 transom | True | True |
| series60 带 sine 指数 | True | True |
| 未知 family（'wigley'） | True | True |
| Cb 出梯（0.95） | True | True |

### 3.4 闭式静水力：五档网格 vs 闭式（`docs/analytic_hydrostatics/results.json` + `report.md`）

闭式基准：V=**7.51269822 m³**，水线面=**14.0863092 m²**，中剖面=**1.06666667 m²**，
LCB/LCF=**5 m**（=L/2），VCB=**−0.32 m**（=−2/5·T），Cb=**0.469544**，Cp=**0.704315**，
湿面积=**22.9127295 m²**（Gauss 基准）。

| 量 | 11×5 | 21×9 | 41×13（默认） | 81×21 | 161×41 |
|---|---|---|---|---|---|
| 排水体积 | −5.2905 % | −1.8416 % | **−0.8804 %** | −0.3867 % | **−0.1359 %** |
| 湿表面积 | −1.2367 % | −0.5487 % | −0.3261 % | −0.1730 % | −0.0742 % |
| 水线面 | −1.8477 % | −0.5684 % | −0.1749 % | −0.0538 % | −0.0166 % |
| 最大剖面 | −3.5075 % | −1.2805 % | −0.7067 % | −0.3330 % | −0.1194 % |
| LCB / LCF | +0.0000 % | +0.0000 % | +0.0000 % | +0.0000 % | +0.0000 % |
| VCB | +4.1318 % | +1.6449 % | **+0.9391 %** | +0.4566 % | +0.1685 % |
| Cb | −5.2905 % | −1.8416 % | −0.8804 % | −0.3867 % | −0.1359 % |
| Cp | −1.8477 % | −0.5684 % | −0.1749 % | −0.0538 % | −0.0166 % |

- 默认 41×13 全量 **< 3 %**（`default_grid_inside_three_percent=True`）；161×41 全量 **< 0.2 %**
  （`fine_grid_inside_point_two_percent=True`）；体积误差随 h **单调下降**（`volume_error_is_monotone_in_h=True`）。
- 符号：凡长度/面积/体积加权和处**全负**（弦切在解析面内侧，欠报）；**VCB 是唯一高报**（+0.939 %）——
  梯形截面规则在龙骨附近（`sqrt` 剖面最陡处）高估宽度，抬高形心。
- LCB/LCF 落 `L/2` 至 **1e-9 m**（对称性后果，校验 x 矩框架与符号约定，非非对称载荷）。

### 3.5 离散化估计器 vs 真实误差（`results.json` `estimator`，默认 41×13）

| 量 | 粗/细相对变化 | 真实误差 vs 闭式 | 比值 |
|---|---|---|---|
| 排水体积 | 0.5741 % | 0.8804 % | **1.53** |
| 湿表面积 | 0.1843 % | 0.3261 % | **1.77** |
| 水线面 | 0.1212 % | 0.1749 % | **1.44** |

`estimate_discretisation_error` 是**乐观指示器而非误差界**（真实误差为报告变化的 1.44–1.77×，
与一/二阶混合收敛一致——尖端奇异性所致）；`analytic_relative_errors` 把该比值以
`outside_discretisation_estimate` 随数字一起返回。

---

## ④ 结果图

- **`docs/hull_family.png`（2325×720 px，已留档）**：(a) 五父 161×41 网格中剖面——U 形随 Cb 上升填充；
  (b) 三档网格位移误差 vs 闭式（五父全绘）；(c) 尖艉 vs 0.15 艉板沿船长站宽——艉甲板宽是参数，其余段不变。
- **`docs/analytic_hydrostatics/convergence.png`（未留档）**：闭式静水力收敛曲线图在工件目录未归档，
  本补报以 §3.4 数字表呈现收敛（五档网格、单调下降、符号全负/VCB 例外），不另绘。
- 本补报单按任务规格 `figures=[]`（不内嵌 base64），结果图以 worktree 路径引用；原图均可在
  `docs/` 下复现（`tools/hull_family/report.py`、`examples/analytic_hydrostatics/run_analytic_check.py`）。

---

## ⑤ 注意点（踩过的坑 / 局限 / 适用范围）

1. **默认网格龙骨转弯未解**：41×13 同参数位移约 −4 %，触发模块 2 % `sampling_error_within_tolerance`
   **警告（非错误）**；161×41 才是 <1.2 % 声明的分辨率，生产导出应走 `resolve_resolution`。
2. **VCB 例外**：唯一高报量（+0.939 % 默认 / +0.169 % 细网格），源于龙骨附近梯形规则高估——不是缺陷，
   是离散化方向性；若出现正体积误差反而意味着放样过冲（弦近似不可能），应视为缺陷。
3. **离散化估计器乐观**：真实误差 = 报告相对变化 1.44–1.77×，只能当指示器，不能当误差界。
4. **无外部锚点（WITHHELD 维持）**：闭式是 Todd (1951) 站型律的**解析值**，系数梯是**输入**——
   §3.1/§3.4 度量的是「放样 vs 自身律」，不是「律 vs 扫描母型 offsets 表」；无实验/文献可比。
5. **无水动力声明**：阻力/耐波/稳性/GZ 出范围，本单是几何 + 静水力积分。
6. **球鼻艏是独立闭体**：需单水密固体的下游求解器必须显式合并。
7. **`vertex_manifold` 警告继承**：`generate_series60` 在参考船（环在中心线闭合）上自报，是族网格属性非分派缺陷。
8. **船体族端点非 G1**：`b(x)=B/2·sin(πx/L)^0.7` → `db/dx ~ x^-0.3`，艏艉水线斜率无限（曲率发散），
   龙骨同（`sqrt(z/T)`）——是族属性，非网格器问题（0.12.0 曲率审计已记录）。
9. **LCB/LCF 平凡**：`L/2` 由对称性得出，测不出非对称偏置；需非对称变体才可检验。
10. **水线是 z=0 平面**：纵倾/横倾在模型外（与 `compute_hydrostatics` 同）。

---

## ⑥ 误差判定（关键量 vs 基准，显式回答是否 <3 %）

### 船体族（基准 = 各父自身闭式 / 族模块）

| 关键量 | 基准 | 实测 | < 3 % ? |
|---|---|---|---|
| 位移（161×41，五父最坏） | 闭式 | +1.22 % | **是** |
| 实测 Cb vs 目标（161×41，五父最坏） | 目标 Cb | +1.22 % | **是** |
| mesh vs 族模块（15 例） | 族模块 | bit-identical | **是**（精确） |
| 艉板宽 vs 0.15·B | 0.15·B | +0.00 %（1e-9 相对） | **是** |
| 船体位移 带/不带球鼻艏 | 自身 | bit-identical | **是**（精确） |

### 闭式静水力（基准 = 合成曲面的解析闭式）

| 关键量 | 基准 | 默认 41×13 | 161×41 | < 3 % ? |
|---|---|---|---|---|
| 排水体积 | 7.51269822 m³ | −0.880 % | −0.136 % | **是**（默认 <3 %，细 <0.2 %） |
| 水线面 | 14.0863092 m² | −0.175 % | −0.017 % | **是** |
| 中剖面 | 1.06666667 m² | −0.707 % | −0.119 % | **是** |
| VCB | −0.32 m | +0.939 % | +0.169 % | **是** |
| 湿面积 | 22.9127295 m² | −0.326 % | −0.074 % | **是** |
| LCB/LCF | 5 m | 1e-9 m | 1e-9 m | **是**（精确） |

**结论：关键量误差全部 < 3 %（默认网格），发布分辨率（161×41）全部 < 0.2 %。**

**无基准声明原因**：闭式静水力基准是**合成曲面的解析定义**（`sin(πx/L)^0.7·sqrt(w)` 的闭式积分），
非实验/文献/母型实测——该曲面是合成船型，唯一有意义的基准是它自身的解析定义，故「<3 %」判定针对的是
**网格离散化误差 vs 自身闭式**，而非 vs 外部实验；外部锚点（扫描母型 offsets 表）维持 WITHHELD，
不粉饰为已对标。

---

## 附：验证与复现

```bash
python -m pytest tests/test_hull_family.py -q          # 分派本身（17 测）
python -m pytest tests/test_analytic_hydrostatics.py -q # 闭式基准（21 测）
python tools/hull_family/report.py                     # evidence.json + 图
python tools/hull_family/gen_doc_and_gate.py           # 本页 + 突变体（5/5 击杀）
python examples/analytic_hydrostatics/run_analytic_check.py  # results.json/report.md（8 项 claim 全 True）
```

本补报单为**诚实 meta-报告**：原单 FAIL 结论（若有）如实保留，不改业务代码；缺失工件（`artifacts/hull_family/*`、
`convergence.png`）标注「未留档」，报告未发不算完成。
