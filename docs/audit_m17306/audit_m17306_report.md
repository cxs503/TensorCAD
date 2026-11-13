# 交差核验：m17250 · GAP-CAD-T1-IGES-PIPE（IGES/STEP → 水密网格进口管道）

**审计单** m17306 · **被审单** m17250（`succeeded`）· **仓** TensorCAD · **判定**：**L1**（声明 L1，机器门复核口径见末节）
**审计口径**：只核验不重做。所有数字来自 `artifacts/audit_m17306/checks.json`（`tools/audit_m17306_crosscheck.py`）与 `gold_quick.json`（`tools/audit_m17306_publish.py measure`），零手抄。
**独立性限制**：审计者与被审者同为平台智能体，本核验只依赖可复跑命令 + 机器判据（`tensor_benchmark_cases`、pytest、grep、`git ls-tree`），不采纳报告自述。

## 1. 案例（被审对象是什么）

m17250 声称补上 GAP-CAD-T1-IGES-PIPE：TensorCAD 原先只读「无裁剪 NURBS patch set」（`tensorcad/geomquality/io.py:625 read_iges`，同文件 `:652` 明写 `IGES trimmed surface (type 144) is out of scope`），因此工业界唯一被接受的输入（裁剪 B-rep）进不来。交付声明 = 新增 3 个模块（record/entity 图 → trimmed faces → 反演 → 三角化 → 修复 → 7 道门）+ CLI + 测试 + 证据运行，输入集 4 个文件：KP505 螺旋桨（第三方 CATIA V5-6R2012 导出，1,042,550 B，唯一真·trimmed 文件）、KVLCC2-class 与 Series 60（本仓自己的 `write_iges` 写出的无裁剪船体）、`wigley_iii.iges`（存在于 inputs 但未进 provenance、未进 manifest、未参与运行）。

被审问题的实质只有一句：**裁剪面能不能变成水密网格。**

## 2. 方法（怎么核的）

1. 取证：`git ls-tree -r HEAD` 与 `git log --all -- <path>` 核交付件是否在任何 ref 里；逐件 sha256/字节数。
2. 基准裁判：`SELECT case_name, metric, tolerance FROM tensor_benchmark_cases WHERE anchor_slug='GAP-CAD-T1-IGES-PIPE'`。
3. 复跑：把 m17250 整树复制到 `/tmp/audit_m17250`（不动原树），跑它自己写的 `repro` 命令与测试，比对工件数值是否逐位一致。
4. 保真度量化（报告没有的量）：把导入网格体积与 `tensorcad.geomquality.baselines` 的解析船体体积（series60 L10/B1.0/T0.625/Cb0.60；KVLCC2 L320/B58/T20.8/Cb0.8098）对齐，得出往返误差与反算 Cb。
5. 集成面：在生态仓 `/root/tensor-ecosystem/TensorCAD` grep 新模块的调用方。
6. 文档自述 vs 实测逐条对照（下表 3.3）。

## 3. 计算结果

### 3.1 三件套是否存在（run_artifacts / benchmark / repro_cmd）

| 核验项 | 实测 | 结论 |
|---|---|---|
| 交付件存在性 | 16/16 路径存在（源码 3 + 测试 1 + 文档 1 + runs 6 + artifacts 3 + 输入 1），`iges_trim.py` 46,871 B / `trim_mesh.py` 43,814 B / `import_cli.py` 6,524 B | 存在 ✅ |
| 是否入库 | `in_main_HEAD=False` 与 `committed_anywhere=None` **16/16**；`git branch --list agent/m17250` → 与 main 同为 `404be46`，`git log --all -- tensorcad/geomquality/iges_trim.py` 为空 | **未提交 ❌**（全部交付件是工作树里的 untracked 文件，任何一次清理即蒸发；`succeeded` 与「可交付」不等价） |
| run_artifacts | `kvlcc2_class.stl` 451,911 B、`series60.stl` 448,035 B 非空；`KP505_modelScale_250mm.stl` **49 B（只有 `solid`/`endsolid` 头，0 三角面）** | 部分 ✅（trimmed 分支无网格工件） |
| benchmark | `tensor_benchmark_cases` 共 46 行，`anchor_slug='GAP-CAD-T1-IGES-PIPE'` **0 行**（全表仅 1 行与 CAD 沾边，属 TensorFVM） | **无基准 ❌** → L4 不可能 |
| repro_cmd | 见 3.2，5 条命令全部跑通、数值逐位一致 | ✅（KP505 单文件在负载 load≈190 下需 78–115 s，整条 3 文件命令超过本机 120 s 单命令预算，属吞吐不是正确性） |

### 3.2 复跑（同机、同代码、同输入）

| 文件 | 原工件 | 复跑 | 一致 |
|---|---|---|---|
| series60.iges | 三角面 2228，体积 3.693120820714208 m³ | 2228 / 3.693120820714208 | ✅ 逐位一致 |
| kvlcc2_class.iges | 三角面 2244，体积 307508.1764356243 m³ | 2244 / 307508.1764356243 | ✅ 逐位一致 |
| KP505_modelScale_250mm.igs | route=trimmed，tessellated **0**/4，refused 4，三角面 0 | 同：0/4 refused，0 三角面，逐面残差 1.001e2–1.456e2 vs 容差 1.518e-4 | ✅ 拒绝行为可复现 |
| pytest `tests/test_iges_trim_import.py` | 文档称 10 passed / 1 xfail / 1 skipped | 10 passed, 1 xfailed, 1 skipped（skip 原因：`tests/test_iges_trim_import.py:203` 依赖 `runs/gap_t1_iges/evidence.json`） | ✅ 一致 |

### 3.3 文档自述 vs 实测（差异即审计发现）

| 文档位置 | 自述 | 实测 | 判定 |
|---|---|---|---|
| `docs/gap_cad_t1_iges_pipe.md:71` | 「4 faces tessellated, 4 refused」 | `KP505...report.json` `faces_tessellated=0`、`faces_refused=4` | **虚标**（把 0 写成 4） |
| `docs/gap_cad_t1_iges_pipe.md:57` | 表格 KP505 行「4 read, 4 refused」 | `faces_in_file=4`，被读成面并落成网格的是 0 | 措辞掩盖了 0 输出 |
| `docs/gap_cad_t1_iges_pipe.md:7,51` | 数字来自 `runs/gap_t1_iges/evidence.json` | 该文件**不存在**（`run_evidence.py:109` 会写它，工作树无），因此 :82-84 的金标准数字无工件背书 | 证据链断 ❌ |
| `docs/gap_cad_t1_iges_pipe.md:82` | `DTMB5415_Hull_exact.stl` 164,926 triangles | 164,926 faces ✅（494,778 顶点；文档 §1 写的「828,794 vertices」与实测不符） | 数字部分成立、顶点数错 |
| `docs/gap_cad_t1_iges_pipe.md:83-84` | gold 体积 30,368 m³，vs 水密 `DTMB5415_Hull.stl` 19,568.7 m³ | 30,368.142 m³ / 19,568.704 m³ ✅；但 trimesh 判定**两个文件都不是 watertight**（`is_watertight=False`，winding 一致） | 体积成立；「水密」定名不成立 |
| 全文 | 无误差判定节、无 Cb/L/B 校核 | 往返体积误差从未给出（见 3.4） | 六要素 ⑥ 缺 ❌ |

### 3.4 保真度（报告没有、但这是它能被量化的唯一质量轴）

| 船体 | 解析体积（baselines） | 导入网格体积 | 偏差 | 基线 Cb | 网格反算 Cb |
|---|---|---|---|---|---|
| Series 60 | 3.7500000055 m³ | 3.6931208207 m³ | **−1.517 %** | 0.6000 | 0.5909 |
| KVLCC2-class | 312,621.67 m³ | 307,508.18 m³ | **−1.636 %** | 0.8098 | 0.7966 |

两例同号、同量级偏低，说明不是随机噪声：是 `samples u=24, v=12` 的弦高误差在三面接缝与艏艉收敛段系统性丢失体积。**在 3% 内可接受，但没有任何 G1–G7 门检查它**（G5 只查「体积 > 0」），报告也从未声明，属「门通过 ≠ 几何正确」。

## 4. 图

`artifacts/audit_m17306/figure.png`：左=两例 IGES 往返体积相对解析船体的偏差（−1.52 %、−1.64 %，虚线 ±3 %）；右=由水密网格反算的 Cb 对基线 Cb（0.5909 vs 0.6；0.7966 vs 0.8098）。图由脚本从 `checks.json` 直接绘制，无手工数字。

## 5. 注意点（坑与边界）

* **「succeeded」不代表入库**：m17250 的三件套 16/16 都是 untracked，branch 与 main 同 commit。核验前必须 `git log --all -- <path>`，否则结论会建立在随时消失的文件上。
* **判据排序**：先查 `tensor_benchmark_cases` 有没有行，再读报告。本例基准行 0，报告写得再好也只能到 L2。
* `run_evidence.py` 的金标准审计（164,926 面的纯 Python `audit_mesh`）在本机 load≈190 下超过 15 min 未返回，本单判 **not_run**，改用 trimesh 的 O(n) 面数/体积/水密判定做交叉验证（3.3）；慢路径不作为判级依据。
* 绿的 G1–G7 全部来自**无裁剪**路线，而 `io.read_iges` 早已能读 entity 128 并 sew（`io.py:625`），所以这两条「PASS」不是新功能达成，是旧能力套新审计壳。
* 测试面覆盖：`tests/test_iges_trim_import.py` 对 tokenizer/反演/耳切/审计/修复/焊接逐项有闭式断言（面反演 1e-12、方块体积 8.0 m³ 1e-12、Euler=2），但**没有任何一条测试把 entity 144/143 文件走成网格**——trimmed 分支的输入→输出在整仓从未被自动化验证过。

## 6. 误差判定（逐项，显式回答「是否 < 3%」）

| 关键量 | 基准 | 实测 | 误差 | < 3% ? |
|---|---|---|---|---|
| Series 60 网格体积 | 3.7500000055 m³（`baselines` 解析） | 3.6931208207 m³ | −1.517 % | **是**（但无门、未声明） |
| KVLCC2-class 网格体积 | 312,621.67 m³（同上） | 307,508.18 m³ | −1.636 % | **是**（同上） |
| Series 60 Cb | 0.6000 | 0.5909 | −1.52 % | 是 |
| KVLCC2-class Cb | 0.8098 | 0.7966 | −1.64 % | 是 |
| gold_exact 体积 | 文档 30,368 m³ | 30,368.142 m³ | +0.0005 % | 是（文档数字成立） |
| plain 体积 | 文档 19,568.7 m³ | 19,568.704 m³ | +0.00002 % | 是 |
| **KP505 trimmed 面成功率** | 需 4/4 面 → 水密网格（锚点要求） | **0/4，0 三角面** | 功能未达成（不是精度问题） | **否** |
| 锚点基准误差 | `tensor_benchmark_cases` 无该锚点行 | 无 ref/got/tolerance | 无可比基准 | **不适用（无基准，非「通过」）** |

诚实结论：**精度轴无从失分（因为没有基准行），功能轴明确失分（4/4 面被拒）。**

## 7. 判级与理由（L0–L5）

**L1。** 依据 owner 判据「级别 = 功能实质」：

* 不是 L2：锚点功能 = 裁剪 IGES/STEP → 水密网格。唯一第三方 trimmed 文件 KP505 **4/4 面被拒、0 三角面**，且整仓**没有一条测试把 144/143 走成网格**，即 trimmed 分支从未有过一次「输入→网格输出」的成功执行。L2 要求的「输入→计算→输出」在锚点所指功能上不成立。
* 不是 L0：新码 2,354 行真实存在并真实执行（parser 读出 entity census `{100:8,102:12,108:2,110:12,120:2,124:8,126:53,128:4,142:8,144:8}`，逐面反演残差、修复链、7 门全是真算），组件数学有闭式验证。
* 不是 L3：`grep -rn 'trim_mesh|iges_trim|import_cli' /root/tensor-ecosystem/TensorCAD/tensorcad` → **0 个外部调用方**。
* 不是 L4：基准行 0；报告的绿道是自往返（`write_iges`→`read_iges`）且属旧能力。
* 若把 scope 收窄到「无裁剪 patch set → 水密网格」，新审计/修复/门机件可达 L2（本审计记为 `level_if_scoped_to_untrimmed_route: 2`）；但那不是本锚点。

**机器门预期**：声明 L1 → `judge` 走 `SCHEME[1] = ["tests"]` 通过 → `gate_level=1` → 原单 m17250 为 `succeeded` 且 `gate_level≤1` → 触发返修轮（正确方向：补 142 参数空间路径 + 144/143 端到端测试 + 有效基准行）。若判级被机器复核降到 0（invalid），则说明 evidence 键不被识读，问题在本审计的证据形态而非被审对象。

**返修最小集**（建议写进 successor 单）：① entity 142 的 (u,v) 参数空间路径，用缝合闭合（free edges）作约定自校验；② 用一把真实 trimmed 文件（哪怕再导出一份）走通 4/4 面并给出体积/面积与解析值对照；③ 在 `tensor_benchmark_cases` 注册基准行（metric=面成功率/体积误差，tolerance 显式）；④ 给 `q`（采样密度）做收敛表，让 1.5% 系统偏差变成可关闭的误差项；⑤ 提交入库（当前 0 commit）。

## 8. 结论

m17250 交的是一条**组件级成立、锚点级未达成**的管道：解析器、反演、三角化、修复、7 门都是真代码且可复现，但「裁剪面 → 水密网格」这一句在实际文件上 0/4 面通过，报告用「4 faces tessellated」把它写成了通过，且全部工件未入库、基准未注册、往返 1.5–1.6% 体积偏差未被任何门或文字覆盖。故判 **L1**，`reopened` 方向（返修）为机制上正确的结果。
## 9. 落表与交付物

* verdict 已落库：`tensor_maturity_verdicts` → `('TensorCAD','m17250 GAP-CAD-T1-IGES-PIPE', 17306, 1, …, 'board-audit')`，SQL 见 `docs/audit_m17306/verdicts.sql`。
* 机器门自检：直接调用 `/root/tp-agent/maturity_gate.py:judge(1, evidence, 'TensorCAD', 17306)` → 返回 `(1, [])`，即无缺项、`gate_level=1`；原单 m17250 为 `succeeded` 且 `gate_level≤1`，机制将开返修轮（`reopened`）。
* 证据快照：`artifacts/audit_m17306/snapshot/`（11 个文件，含 `SHA256SUMS.json`）——因被审交付件在该工作树内无 commit，字节级快照用于防证据蒸发。
* 原始数据：`artifacts/audit_m17306/checks.json`、`gold_quick.json`、`logs/pytest_iges.log`、`logs/repro_hulls.log`、`figure.png`。
* 平台报告：`/tensor/reports` id=16305。
