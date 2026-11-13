## m17405 · 交差审计：m17357 `GAP-CAD-D3-GEOM-CONTRACT`（TensorCAD）

**判级 L2（独立功能实现）· approve with notes · verdict id=543 · gate 预演 judge(2)→2 confirmed**

### ① 案例介绍
`GAP-CAD-D3-GEOM-CONTRACT` 要求 TensorCAD 交付**几何契约文档 + 运行时校验器**：把坐标系/单位/比例/取向
四条约定做成机器可执行检查，对 8 个现存 STL 资产跑出合规表，并实证违规可自动检出。物理背景=船舶 CAD 网格
进入阻力/耐波/结构求解器前的**互换性契约**（单位错标 1000×、法线反转、左手系都会让下游阻力/水动力结果失效
而不报错）。对象工作树 `/root/tp-agent/work/m17357`，交付 6 个文件（校验器 347 行、文档 149 行、测试 231 行、
两个工具脚本、9 STL + 合规表），参数表=箱形驳 30×8×5 m、单位球 r=1 m、Wigley L=10 B=1 T=0.625 m、
Series 60 L=20 B=3 T=1.5 m（Cb 0.60/0.75）。

### ② 怎么做的
1. 取证 `tensor_dev_missions` id=17357（`head_commit=null`、8 路径未提交、`server_quality_audit.gate_valid=false`）。
2. 基准裁判 `tensor_benchmark_cases WHERE anchor_slug='GAP-CAD-D3-GEOM-CONTRACT'` → **0 行**（无 ref/tolerance）。
3. 自建独立探针 `docs/board_audit/m17405/audit_probe.py`（不 import 对象夹具）：解析体积对照 / 5 组自建违规 /
   3 组沉默窗口反例 / CLI 退出码 / 复现哈希。
4. 复现 `python3 tools/run_contract_compliance.py` → 前后 sha256 `ee18be55cb3dbdc1…` **identical=True**。
5. 测试 `pytest tests/test_geom_contract.py -q` → **22 passed / 0 failed**。
6. 集成面 grep `validate_contract` → 仅 `tensorcad/__init__.py:457` 再导出，**零调用方**。

### ③ 计算结果
对象 8 资产合规表（程序化读取对象工件）：

| 资产 | 判定 | 有符号体积 m³ | 违规 |
|---|---|---|---|
| box_barge.stl | PASS | 1200.0 | — |
| inverted_wigley.stl | FAIL | -2.771222144 | [coordinate_system/warning] coordinate frame is unspecified; right-handed is assumed; [orientation/error] closed mesh has non-positive signed volume (inward normals or inverted) |
| left_handed.stl | FAIL | 2.771222144 | [coordinate_system/error] left-handed coordinate system is not supported |
| series60_cb60.stl | PASS | 52.014385507 | — |
| series60_cb75.stl | PASS | 64.753744466 | — |
| sphere_r1.stl | PASS | 4.168925864 | — |
| wigley.stl | PASS | 2.771222144 | [coordinate_system/warning] coordinate frame is unspecified; right-handed is assumed |
| wrong_units.stl | FAIL | 2.771222144 | [coordinate_system/warning] coordinate frame is unspecified; right-handed is assumed; [scale/error] bounding-box diagonal 0.0100693 m is below minimum 0.05 m |

审计方独立体积基准（ref=闭式解析解，got=从交付 STL 用散度定理重算）：

| 算例 | 公式 | ref | got | 相对误差 | <3% |
|---|---|---|---|---|---|
| box_barge closed box | `V = L*B*D = 30*8*5` | 1200.000000 | 1200.000000 | **+0.0000%** | 是 |
| sphere r=1 faceted (cube-to-sphere, 32x64) | `V = 4/3*pi*r^3` | 4.188790 | 4.168926 | **-0.4742%** | 是 |
| Wigley parabolic hull L=10 B=1 T=0.625 | `V = (4/9)*L*B*T` | 2.777778 | 2.771222 | **-0.2360%** | 是 |

红队 5 例（审计方自建网格，不共用对象夹具）：

| 违规用例 | 命中约定 | passed |
|---|---|---|
| inward normals (closed tetra, reversed winding) | orientation | False |
| left-handed declared frame | coordinate_system | False |
| unrecognised unit token | units | False |
| below-minimum scale (1 mm tetra declared m) | scale | False |
| above-maximum scale (1000 m tetra) | scale | False |

### ④ 结果图
![contract audit](attachment:contract_audit.png)
左：8 资产 × 4 约定违规矩阵（红=error，黄=warning，空白=通过），5 PASS / 3 FAIL 与对象自报一致；
右：交付 STL 体积 vs 闭式解析体积奇偶图（±3% 带 + 逐点误差标注），三点全部落在 ±0.5% 内——
证明校验器度量的是真实几何，不是回填常量。

### ⑤ 注意点
* **未提交**：8 个交付文件全为对象工作树未跟踪文件，主检出 `/root/tensor-ecosystem/TensorCAD` 无
  `geom_contract.py`，其他仓/CI 现在无法 import 该能力（已按字节快照进 `docs/board_audit/m17405/snapshot/`）。
* **无六要素报告**：`tensor_reports` 无 mission_id=17357 行（report_id 空），本报告代为补齐并留痕欠账。
* **坐标系检查是声明串检查**：源码 docstring 已声明；镜像网格 + 声明 right-handed 的开放面网格零违规通过。
* **比例检查有盲窗**：100× 单位错标只要对角线仍落在 [0.05, 500] m 内即通过。
* **开放面取向只到 warning**：恰好 50% 反向面零违规通过（阈值为 ratio < 0.5）。
* 适用范围：STL 不携带单位/坐标系元数据，结论只对"调用方声明正确"的场景成立。

### ⑥ 误差判定（是否 <3%）
| 关键量 | 基准 | 误差 | <3% |
|---|---|---|---|
| 箱形驳体积 | 解析 L·B·D=1200 m³ | +0.0000% | 是 |
| 单位球体积 | 解析 4/3πr³=4.18879 m³ | -0.4742% | 是 |
| Wigley 体积 | 解析 (4/9)L·B·T=2.77778 m³ | -0.2360% | 是 |
| 登记基准案 ref/got vs tolerance | `tensor_benchmark_cases` 该锚点 0 行 | 无法计算 | **无基准案 → L4 不成立（honest）** |

越界观察（不计入本单判定，校验器忠实报告网格体积，偏差在资产生成器侧，建议另开单用
`tensorcad.series60.analytic_volume` 仲裁）：

| 算例 | ref=Cb·L·B·T | got | 误差 | 隐含 Cb | 声明 Cb |
|---|---|---|---|---|---|
| series60_cb60 vs Cb*L*B*T (L=20 B=3 T=1.5) | 54.000 | 52.014 | -3.677% | 0.5779 | 0.60 |
| series60_cb75 vs Cb*L*B*T (L=20 B=3 T=1.5) | 67.500 | 64.754 | -4.069% | 0.7195 | 0.75 |

### 判级与处置
L2 逐键：`impl_files` ✅（347 行实现，无占位）· `entry` ✅（库 :240/:268 + CLI :282）·
`run_artifacts` ✅（8 报告 + 9 STL，字节级可复现）· `tests` ✅（22 绿且含正反例不变式）·
`integration_callers` ❌ 无调用方 → 拒 L3 · `benchmark` ❌ 锚点 0 行 → 拒 L4 · `external_anchor` 不申。
反例 3 条已登记（`evidence.checks_not_run` 同单）。后续：commit 落库 → 接网格导入口挣 L3 →
登记解析锚点基准案（箱驳 1200 m³ / Wigley (4/9)LBT，tolerance 1%）挣 L4 → 收紧开放面取向阈值。

**Honest FAIL 项**：无登记基准案、无六要素报告、交付未提交——三项如实记录，不粉饰。
