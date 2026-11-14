## 1. 案例介绍：m20416 交差核验（XVAL-GEOM-VOLCROSS）

**物理背景。** 排水体积 V 是所有水动力量的归一化分母（Cb、阻力系数、载重、稳性）。
单一积分器无法自证：两个共享同一错误概念的实现会互相一致而与真值无关。锚点
XVAL-GEOM-VOLCROSS 要求的正是"同体积、异方法类别、异代码库"的交叉互证。

**被审对象** 是 m20416 在 TensorCAD 工作树中产出的证据包（自研 NURBS/散度体积积分
vs 独立实体引证）。本单**只核验不重做**：审计在只读方式下对 m20416 工作树取证，
修复性探针一律在 `/tmp/m20464_probe` 的**一次性副本**上执行。

**几何与参数**（4 条船，参考值证据分级各不相同）：

| 船 | L (m) | B (m) | T (m) | Cb | 参考体积 (m³) | 参考证据级 |
|---|---:|---:|---:|---:|---:|---|
| Wigley | 10.0 | 1.0 | 0.625 | — | 2.7777777778（解析闭式 L·B·T·n·m/((n+1)(m+1))） | REAL（代码内求值，非抄录） |
| KCS | 230.0 | 32.2 | 10.8 | 0.651 | 52030.0 | WITHHELD（台账 recalled literature, verified=false） |
| JBC | 280.0 | 45.0 | 16.5 | 0.858 | 178369.9 | WITHHELD（同上） |
| Series60 | 100.0 | 13.3333 | 5.3333 | 0.6 | 4266.6293334 | WITHHELD（closure Cb·L·B·T，非外部验证） |

几何先行校验（Cb/L/B）用台账字段与实际网格体积互算：生产分辨率 241×49 下
KCS 网格 V/LBT = 0.6453（台账 Cb 0.651，差
0.88 %），
JBC 0.8486（0.858，差
1.10 %），
Series60 0.5949（0.6，差
0.86 %）——
即"偏差"全部落在放样多边形化，而不是台账几何本身，与文献参考值的 3~4 位有效数字
不确定度同量级。

## 2. 怎么做的：审计流程与关键设置

1. **取证**（只读）：读 m20416 工作树三件套（impl / entry / run_artifacts / tests），
   `git status`、`git cat-file HEAD:tensorcad/volume_xval.py`、`py_compile`、
   `pytest --collect-only`，并查 `tensor_dev_missions.id=20416` 的 result/evidence。
2. **裁判**：`SELECT case_name,metric,tolerance FROM tensor_benchmark_cases
   WHERE anchor_slug='XVAL-GEOM-VOLCROSS'` → **0 行**（全表
   131 行）⇒ **无机器裁判行**，按仓内既有纪律回落到
   模块自声明带宽 + 工件内参考值，并显式标注"无裁判"。
3. **可复现性探针**：把交付树拷到 `/tmp/m20464_probe`，做**最小机械修补**（删 1 个游离
   token、删 1 处重复 guard、把 20 行块回退 4 空格；不增删任何逻辑），再跑
   `pytest tests/test_volume_xval.py -q` 与 `tools/volume_xval_study.py`，把重生成工件与
   原工件做 560 个叶子的逐项比对。
4. **集成核查**：全树 grep `volume_xval|cross_validate_*`（排除三件套自身）+ `tensorcad/api/`，
   判 L3。
5. **判级**：L0–L5 按"功能实质"而不是文档厚度。
6. **落表**：`tensor_maturity_verdicts`（repo=TensorCAD, scope=m20416 XVAL-GEOM-VOLCROSS,
   mission_id=20464, level=2, verdict_by=board-audit）。

工具链：psql(tensorplatform) / python3.11 / pytest / numpy / trimesh 5.0.0 / matplotlib。
审计脚本：`docs/audit_m20464/audit_m20416.py`、`docs/audit_m20464/land_verdict.py`
（本报告数字全部由这两步产出的 JSON 注入）。

## 3. 计算结果：三件套状态与逐项误差

### 3.1 三件套存在性（机器判定）

| 项 | 判定 | 证据 |
|---|---|---|
| impl 文件 | 存在，**1072 行，但 `py_compile` 失败**：`Sorry: IndentationError: unexpected indent (volume_xval.py, line 1053)` | `tensorcad/volume_xval.py:1043-1073` |
| 版本控制 | **未入仓**：6 个交付文件全部 `git status ??`；`head_commit=null`, `branch=null`；`HEAD:tensorcad/volume_xval.py` 不存在 | `git status --porcelain` |
| run_artifacts | 存在：`artifacts/volume_xval_study.json`（25787 B）、md（9002 B）、png（179994 B） | 见 §3.3 |
| benchmark ref vs got | **无裁判行**（0/131），只能用工件内参考值 | `tensor_benchmark_cases` |
| repro_cmd | 交付态**不可执行**：pytest `--collect-only` rc=2，0 tests；study 工具 ImportError。副本修补后 pytest **13 passed, 1 warning in 22.18s**、study rc=0 | §3.4 |
| 平台报告 | `tensor_reports WHERE mission_id=20416` → **0 行**（未发） | DB |
| 集成调用方 | **0 处**（全树 grep + api 目录） | §3.5 |

### 3.2 参考值证据级

| 船 | V_ref (m³) | 证据级 |
|---|---:|---|
| Wigley | 2.7777777777777777 | REAL |
| KCS | 52030.0 | WITHHELD |
| JBC | 178369.9 | WITHHELD |
| Series60 | 4266.6293334 | WITHHELD |

### 3.3 逐路线数字（工件原值，程序化注入）

| 船 | 路线 | 角色 | 证据 | V (m³) | |Δ| vs 参考 |
|---|---|---|---|---:|---:|
| Wigley | self-surface-sum (divergence) | self | REAL | 2.7712 | 0.2360 % |
| Wigley | trimesh (third-party, file round trip) | external | REAL | 2.7712 | 0.2360 % |
| Wigley | ray-stabbing (this module, numpy) | external | REAL | 2.7701 | 0.2778 % |
| Wigley | B-rep solid (gmsh/OCP) | external | WITHHELD | — | — |
| KCS | self-surface-sum (divergence) | self | REAL | 51610.3354 | 0.8066 % |
| KCS | trimesh (third-party, file round trip) | external | REAL | 51610.3354 | 0.8066 % |
| KCS | ray-stabbing (this module, numpy) | external | REAL | 51585.4984 | 0.8543 % |
| KCS | B-rep solid (gmsh/OCP) | external | WITHHELD | — | — |
| JBC | self-surface-sum (divergence) | self | REAL | 176417.3901 | 1.0946 % |
| JBC | trimesh (third-party, file round trip) | external | REAL | 176417.3914 | 1.0946 % |
| JBC | ray-stabbing (this module, numpy) | external | REAL | 176369.0045 | 1.1218 % |
| JBC | B-rep solid (gmsh/OCP) | external | WITHHELD | — | — |
| Series60 | self-surface-sum (divergence) | self | REAL | 4230.0569 | 0.8572 % |
| Series60 | trimesh (third-party, file round trip) | external | REAL | 4230.0569 | 0.8572 % |
| Series60 | ray-stabbing (this module, numpy) | external | REAL | 4227.9711 | 0.9061 % |
| Series60 | B-rep solid (gmsh/OCP) | external | WITHHELD | — | — |

引擎互证（独立代码库 trimesh 5.0.0，二进制 STL 往返）与仓内散度求和在**同一批三角形**上
一致到 1.0e-8 相对（float32 STL 精度极限）；ray-stabbing（异方法类别：不打面积分，靠射线
内部弦长配对）在 256 ray/side、自测标准误 ×3 的带宽内一致。

### 3.4 放样收敛（KCS，lon/lat 同向加密）

| nx×nz | 三角形 | V (m³) | 相对参考 |
|---|---:|---:|---:|
| 61×13 | 2998 | 50107.609 | -3.6948 % |
| 121×25 | 11758 | 51118.012 | -1.7528 % |
| 241×49 | 46558 | 51610.335 | -0.8066 % |
| 481×97 | 185278 | 51850.742 | -0.3445 % |

观测阶 **1.18**（log-log 拟合），Richardson 外推 V→h→0 = **52041.79 m³**
对参考 52030.0 差 **+0.0227 %** ⇒ 生产分辨率的 −0.81 % 是
放样亏损，不是积分错误。

### 3.5 判级（L0–L5，按功能实质）

| 级 | 内容 | 判定 | 依据 |
|---|---|---|---|
| L0 | 文档成果 | 过程 | 8 节 md + png 存在但本身不计级 |
| L1 | 数学/单元验证 | **通过（修补后 13/13；交付态 collect error）** | 精确立方体、内切球、实测不确定度带、NURBS 求积收敛阶、Richardson 极限入带（13 项测试，修补后 13 passed） |
| L2 | 独立功能实现（输入→计算→输出） | **未落地（实质具备）** | 台账船体 → 4 路引擎 → 逐路线体积表 + 归因；工件 179/179 数值叶子可重生成 |
| L3 | 系统集成 | **未达** | 0 调用方、无 REST/solver 接入 |
| L4 | 基准验证（完整算例+误差分析+可复现） | **未达** | 无裁判行；4 个参考中 3 个自声明非外部（WITHHELD）；**交付态不可复现**（未入仓 + 不可导入） |
| L5 | 科研工程成果 | 未达 | 无超出上述的成果 |

**结论：level = 1（返修），不是 2。** 判据按"已落地功能"：L2 的文字判据是"真代码：输入→计算→输出"，而交付态的实现无法 import（`volume_xval.py:1053` IndentationError）、六个交付文件全未入仓（`head_commit=null`/`branch=null`，`HEAD:tensorcad/volume_xval.py` 不存在）——仓库今天不存在可运行的体积互证能力，把"6 行修好就能跑"记成 L2 会把一次坏交付销账。故按交付记 L1，同时把实质证据留全（修补 6 行后 13/13 绿、工件 179/179 数值叶子复现、全部有基准量 <3 %），返修成本极低且路径明确：修补+入仓→L2，补裁判行→L4 候选。

## 4. 结果图

左：四条船、三条真实路线对各自参考值的相对偏差（对数轴，虚线=模块自声明 0.3 % 带）。
右：KCS 放样收敛（log-log，误差随 1/nx 单调下降，Richardson 外推入带）。
B-rep 路线（gmsh/OCP）为 NaN/WITHHELD，故不入图——缺路就是缺路，不用别的路线补画。

## 5. 注意点（踩过的坑/局限/适用范围）

1. **"已 succeeded" 不等于交付可运行**：m20416 的 result 里 `finished_status=succeeded`
   且 `accepted=true`，同一段文本却自认 `IndentationError ... 如实记录为 FAIL`，
   `fail_stage=budget_exhausted+forced_finish`、`server_quality_audit.gate_valid=false`。
   状态位是平台给的，可交付性是树给的——审计只能认后者。
2. **自报成绩会随最后一次坏补丁失效**：自称 `13 passed / 19.78s` 指一个**已不存在**的
   代码状态；本次修补后实测 `13 passed / 22.18s`，可复现但不是在交付态上。
3. **无机器裁判行**：131 行的 `tensor_benchmark_cases` 里没有 XVAL-GEOM-VOLCROSS，本轮
   只能用模块自声明带宽 0.3 % 与工件内参考值；这类"自设带宽"必须同时看参考值证据级，
   否则就是把自证当互证。
4. **3/4 参考值非外部**：KCS/JBC 是 recalled literature（verified=false），Series60 是
   closure Cb·L·B·T（算术，不是测量）。因此 −0.81 %/−1.09 %/−0.86 % 不能归因为几何误差，
   模块自己的 `locate_discrepancy` 也如实输出 `reference not external (cannot attribute)`。
5. **B-rep 路线缺席**：loongarch64 上 gmsh/OCP/cadquery 均不可装（配置索引无对应轮子），
   工件标 WITHHELD 而非悄悄丢掉，这是正确做法；但也意味着"独立实体引证"这一半的锚点
   只由 trimesh + ray-stabbing 承担。
6. **适用范围**：结论只覆盖 m20416 交付的这套实现与其工件；不评价 trimesh 自身正确性
   （只到 1e-8 往返一致），也不评价 KCS/JBC 文献值的真伪。审计未修改被审工作树。

## 6. 误差判定（逐项显式回答"是否 <3%"）

| 量 | 参考 | got | 相对误差 | <3 % ? | 备注 |
|---|---|---:|---:|---|---|
| Wigley V（自研散度） | 2.7777777778 | 2.7712221 | 0.2360 % | **是** | 亦 <0.3 % 自声明带 |
| Wigley V（trimesh 往返） | 同上 | 2.7712221 | 0.2360 % | **是** | 与自研差 ~1e-8 |
| Wigley V（ray-stabbing） | 同上 | 2.7700607 | 0.2778 % | **是** | 在自身采样带内 |
| KCS V（自研/ trimesh） | 52030.0 | 51610.335 | 0.8066 % | **是** | 参考非外部；>0.3 % 带 |
| KCS V（ray） | 52030.0 | 51585.498 | 0.8543 % | **是** | 采样带内 |
| KCS V（Richardson） | 52030.0 | 52041.79 | +0.0227 % | **是** | 外推后入 0.3 % 带 |
| JBC V（自研/ trimesh） | 178369.9 | 176417.390 | 1.0946 % | **是** | 参考非外部 |
| JBC V（ray） | 178369.9 | 176369.004 | 1.1218 % | **是** | 采样带内 |
| Series60 V（自研/ trimesh） | 4266.6293334 | 4230.057 | 0.8572 % | **是** | 参考是 closure，非外部 |
| Series60 V（ray） | 同上 | 4227.971 | 0.9061 % | **是** | 采样带内 |
| 引擎互证（self vs trimesh） | — | 相对差 ~1e-8 | — | **是** | float32 STL 精度极限 |
| B-rep 实体体积 | — | WITHHELD | — | **无基准** | loongarch64 无 gmsh/OCP 轮子，不声明通过 |

全部有基准的量均 <3 %（最大 1.1218 %），
本项目级判据满足；但**模块自设的 0.3 % 带**在 KCS/JBC/Series60 三条船上被超出
（−0.81 %/−1.09 %/−0.86 %），这是工件的主动发现（归因为放样多边形化，外推后消除），
不是被容忍的状态。B-rep 一路无基准，按 WITHHELD 记账，不计入通过项。

**诚实结论**：m20416 的算法与数字站得住（可复现、误差 <3 %、外推后 <0.1 %），但**交付形态
不成立**——未入仓、交付态 `IndentationError`、无裁判行、报告未发、0 集成调用方。
故裁定 **L1**，`maturity_gate` 据此应走 reopened（返修轮），而不是 confirmed（消账）。（改判记录：本单第一步按"工件实质"记 L2，随后按 gate 的机器判据复核——L2 键要求 impl 在仓内可运行——确认该实现不在仓内且不可导入，当场改判 L1；两次判定的数字证据完全相同，差别只在"落地与否"。）
