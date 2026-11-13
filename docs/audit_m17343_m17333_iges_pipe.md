# 交差核验：m17333 · GAP-CAD-T1-IGES-PIPE（IGES/STEP → 水密网格进口管道）

**审计单** m17343 · **被审单** m17333（`succeeded`，返修轮 1）· **仓** TensorCAD · **判定：L1**
**审计口径**：只核验不重做。所有数字来自 `artifacts/audit_m17343/checks.json`（`tools/audit_m17343_crosscheck.py` 程序化采集），零手抄；图 `artifacts/audit_m17343/figure.png` 由 `make_figure.py` 从同一 JSON 直接绘制。
**独立性限制**：审计者与被审者同为平台智能体；判据只用可复跑命令与机器表（`tensor_benchmark_cases`、pytest、trimesh、git、`psql`），不采纳任何报告自述。

## 1. 案例（被审对象是什么）

锚点要求的产业能力是：**裁剪 B-rep（IGES 144/143、STEP）→ 水密三角网格**。TensorCAD 原先只会 sew 无裁剪 NURBS patch-set（entity 128），`geomquality/io.py:652` 明写 `IGES trimmed surface (type 144) is out of scope`；ship-platform 的几何输入正卡在这里。
m17333 的交付 = 提交 `51eab47 m17333: matured mission deliverables`（branch `agent/m17333`），内容为 3 个自写 IGES 输入 + provenance + `reproduce.py` + `runs/gap_t1_iges/*`（3 report/3 STL/manifest/evidence）+ 文档与测试改动。它自报 **L2**（`docs/gap_cad_t1_iges_pipe.md:4`）。

## 2. 方法（怎么核的）

1. **取证**：`git show --stat 51eab47`、`git merge-base --is-ancestor 51eab47 main`、`git status`；逐件 sha256。
2. **基准裁判**：`SELECT count(*) FROM tensor_benchmark_cases WHERE anchor_slug='GAP-CAD-T1-IGES-PIPE'` → **0 / 52**（表内无该锚点行）。
3. **复跑**：在 m17333 树内跑它自己的 `python artifacts/gap_t1_iges/reproduce.py`（3.4 s），比对 3 份 report 的体积/三角数与 3 个 STL 的 sha256。
4. **独立量测**：用 `trimesh`（不是被审方自己的 `audit_mesh`）重新判水密/体积/Euler 数；用 IGES D 段逐记录统计实体类型（文件是 patch-set 还是 trimmed）。
5. **能力面探针**：把全生态唯一第三方 trimmed IGES（KP505 螺旋桨，`/root/tp-agent/work/m15295/.../KP505_modelScale_250mm.igs`，sha256 `0bb54ae7…`，D 段 `{100:16,102:24,108:4,110:24,120:4,124:16,126:106,128:8,142:16,144:16}`）送进同一份 `read_iges_trimmed`/`trim_mesh`，逐环量 trim 反演残差（预算 25 s）。
6. 集成面 `grep`；测试面 `pytest tests/test_iges_trim_import.py`（被审树 vs 审计树两版）。

## 3. 计算结果（数字表）

### 3.1 三件套存在性

| 核验项 | 实测 | 结论 |
|---|---|---|
| 是否入库 | `51eab47` 在 branch `agent/m17333`；`merge-base --is-ancestor` → **NOT_IN_MAIN** | 已提交分支、未进 main ✅/⚠ |
| impl_files | `tensorcad/geomquality/iges_trim.py`(1203 L) `trim_mesh.py`(1020 L) `import_cli.py`(131 L)；**本次交付未改一行 `tensorcad/`**（`git show --stat 51eab47` 里 `tensorcad/` 0 个文件，改的只有 docs/tests/artifacts） | 代码真、但属 m17250 旧交付 ✅ |
| entry | `python -m tensorcad.geomquality.import_cli`；未被 `tensorcad/geomquality/__init__.py` 导出（grep 0），`tensorcad/cli.py` 无 iges 入口（grep 0） | 模块级入口 ✅，公开 API ❌ |
| run_artifacts | `runs/gap_t1_iges/{series60_cb060,series60_cb065,kvlcc2_class}.{report.json,stl}` + `manifest.json` + `evidence.json`；复跑后体积/三角数/STL sha256 **逐位一致**（3/3） | ✅ 真实、可复现 |
| benchmark | 表行 **0/52**；`evidence.json` 里的"benchmark" 是 `Cb·L·B·T`（与写 IGES 的 `baselines` 同一个数生成器）对自写文件的往返差 | 自洽，非基准 ❌ |
| repro_cmd | `python artifacts/gap_t1_iges/reproduce.py` → rc=0，3.39 s，`all_gates_passed=True` | ✅ |

### 3.2 输入到底是什么（entity census，本次审计实测）

| 文件 | 字节 | D 段实体 | provenance | 说明 |
|---|---|---|---|---|
| series60_cb060.iges | 149,121 | `{128: 8}` | ✅ sha 一致 | 4 面 NURBS patch-set，**无 144/143** |
| series60_cb065.iges | 194,724 | `{128: 8}` | ✅ sha 一致 | 同上 |
| kvlcc2_class.iges | 131,301 | `{128: 8}` | ✅ sha 一致 | 同上 |
| wigley_iii.iges | 6,075 | `{128: 8}` | **❌ 不在 provenance、不在 manifest** | 提交进库但未参与运行；实测 2255 三角、**不水密**（G3/G4/G6 全 FAIL） |

三个"真实 IGES"全部由仓内 `write_iges` 从 `geomquality.baselines` 的数学船体写出（`reproduce.py:35-48`），因此是**自往返**，不是第三方工业导出。

### 3.3 精度（唯一有数字的轴；对自洽基准）

| 文件 | 网格体积 | Cb·L·B·T | 偏差 | trimesh 独立判定 |
|---|---|---|---|---|
| series60_cb060 | 3.693121 m³ | 3.7500 m³ | **1.517 %** | watertight=True, winding=True, Euler=2, V=3.693121 |
| series60_cb065 | 3.983409 m³ | 4.0606 m³ | **1.901 %** | watertight=True, Euler=2 |
| kvlcc2_class | 307,508.176 m³ | 312,621.670 m³ | **1.636 %** | watertight=True, Euler=2 |

三例同号偏低（−1.5~−1.9 %），是 u=24/v=12 采样弦高在接缝/艏艉的系统性丢体积；无门约束、无收敛表。

### 3.4 锚点能力探针：第三方 trimmed IGES（KP505）

| 量 | 实测 | 判据 |
|---|---|---|
| 面读入 | 4/4（`BoundedParametricSurface`） | — |
| trim 环反演残差（8 环，实测 4 环后触时间预算） | **126.06、133.97、105.97** mm；仅 face0-loop1 = 3.43e-6 mm | 容差 1.518e-4 mm（=1e-6×extent 151.77） |
| 超容差倍数 | **8.3e5×、8.8e5×、7.0e5×** | — |
| 面通过/落网格 | 被审方前序运行记录（`/root/tp-agent/work/m17250/.../KP505…report.json`）：`faces_tessellated=0, faces_refused=4`，残余 1.001e2~1.456e2，**0 三角面** | 锚点要求 4/4 面 → 未达成 |

m17333 的输入集里**没有 KP505**：唯一能把 trimmed 路线证伪/证成的第三方文件被换成了一个自写 patch-set。

### 3.5 测试面

| 版本 | 结果 | 说明 |
|---|---|---|
| 被审树 `tests/test_iges_trim_import.py`（`51eab47`） | 14 passed, 1 xfailed | 新增了 provenance/hash-pinned/体积回归三条集成测试 ✅；但同时**删掉了** main 版里的断言（`assert len(triangles) >= 4`、`assert len(triangles) == len(L_shape)-2`、`_signed_area(outer)>0`、`assert broken.free_edge_count == 3`）与 xfail 的文档指针 |
| 审计树（main）同文件 | 9 passed, 2 skipped, 1 xfailed | 两条 skip 是等 `runs/gap_t1_iges/evidence.json`（仓库不带工件） |
| 「144/143 走成网格」的测试 | grep 命中 **1**（仅 import 语句），无一条端到端 | trimmed 分支仍无一次成功的输入→输出回归 |

### 3.6 集成面

`grep -rn 'import_cad|trim_mesh' --include=*.py tensorcad/`（除模块自身）→ **0**；`geomquality/__init__.py` 未导出；`tensorcad/cli.py` 无 iges 子命令。→ 不是 L3。

## 4. 图

`artifacts/audit_m17343/figure.png`（`make_figure.py` 从 `checks.json` 直绘）：
左=三例体积偏差 vs ±3 % 线（1.517/1.901/1.636 %，标注"自洽基准"）；
中=KP505 各 trim 环反演残差（log 轴）vs 容差线，3 环高 6 个数量级，绿柱为唯一正常环；
右=成熟度阶梯（L1 达成；L2 仅"无裁剪路线"阴影部分成立；L3–L5 红）。

## 5. 注意点（坑与局限）

* **"真实 IGES"的口径**：自写器写出的 patch-set 也是合法 IGES，但它只证明"自己的写→自己的读"，不证明工业输入能进；判定必须看 D 段实体类型（本例全部 `{128: 8}`，没有 144/143/142）。
* **换掉失败样本 ≠ 修好失败样本**：m17250 的 KP505（0/4 面）在 m17333 的输入集里消失，门全绿；同机复测残差照旧 1e2 量级。
* **测试改写方向要核**：本次提交同时增加集成测试（好）与删除单元断言（坏），只看"14 passed"会漏掉后者；`git show 51eab47:tests/test_iges_trim_import.py` 与 main 版逐行对比即可见。
* **自报 L2 的机制缺口**：`JUDGE_SCHEME[2]` 只核 `impl_files/entry/run_artifacts` 存在性，自写输入的往返能通过机器门；锚点语义（trimmed）不在机器判据里 → 本条必须由审计侧写清。
* **提速**：KP505 的 4 面全量反演 >2 min（face2-loop1 单环 77 s），本单只测到 4/8 环并标 `not_run_time_budget`，剩余环沿用被审方前序报告的逐面记录（同为 1e2 量级，可互证）。

## 6. 误差判定（逐项，显式回答"是否 <3%"）

| 关键量 | 基准 | 实测 | 误差 | <3 % ? |
|---|---|---|---|---|
| series60_cb060 体积 | 3.7500 m³（自洽 baselines） | 3.693121 | 1.517 % | **是**（但无外部基准、无门约束） |
| series60_cb065 体积 | 4.0606 m³（同上） | 3.983409 | 1.901 % | 是（同上） |
| kvlcc2_class 体积 | 312,621.670 m³（同上） | 307,508.176 | 1.636 % | 是（同上） |
| 网格水密性（3 例） | 水密+流形 | trimesh: 3/3 watertight, winding 一致, Euler=2 | 0 | 是 |
| **KP505 trimmed 面成功率** | 4/4 面 → 水密网格 | **0/4**，残差 1.0e2~1.5e2 mm vs 容差 1.5e-4 mm | 功能未达成（非精度问题） | **否** |
| 锚点基准误差 | 表内 0 行 | 无 ref/got/tolerance | — | **不适用（无基准≠通过）** |

## 7. 判级与理由（L0–L5）

**L1。**
* 不是 L0：解析器/反演/三角化/修复/7 门是真代码真执行，复跑逐位一致（体积、三角数、STL sha256），组件数学有闭式单元测试。
* 不是 L2：锚点所指的**功能实质是裁剪 B-rep → 水密网格**。全生态唯一第三方 trimmed 文件 4/4 面被拒、0 三角面（本次独立复测残差 8.3e5 倍于容差），整仓**没有一条把 144/143 走成网格的测试**；三项"真实 IGES"全是自写 patch-set（`{128:8}`），只算自往返。且本次交付**未改一行 `tensorcad/`**，锚点能力相对 m17250 无任何增量。
* 不是 L3：`tensorcad` 内 0 个调用方，未导出，主 CLI 无入口。
* 不是 L4：`tensor_benchmark_cases` 该锚点 0/52 行；"benchmark" 是同一 `baselines` 生成器的自洽差；文档自称外锚 `DTMB5415_Hull_exact` 实际未使用（`evidence.json:gold_standard` 自注 "proxy: Series60/KVLCC2 analytic volumes"）。
* **若把 scope 收窄到"无裁剪 patch-set → 水密网格"**：impl+entry+run_artifacts+可复现+独立水密复核齐备，可达 **L2**（记为 `level_if_scoped_to_untrimmed_route: 2`，与 m17306 前次审计一致）；但那不是本锚点。**`confirmed` 会消掉一个仍然存在的产业缺口，故必须 reopened。**

## 8. 返修最小集（交给下一轮）

1. entity 142/141 的参数空间路径（当前 3-D 反演在 CATIA 导出上残差 ~1e2 mm，只有 1/8 环可反演）；给逐面残差门而不是整文件级。
2. 4/4 面成网格的真实 trimmed 样本一条（KP505 或任何第三方导出），出体积/面积/水密三量。
3. `tensor_benchmark_cases` 注册该锚点行（metric=面成功率/体积误差，tolerance 显式）——否则 L4 不可能。
4. 收敛表：u×v 加密让 −1.5~−1.9 % 的系统偏差变成可关闭误差项。
5. 恢复被删的单元断言；把 `wigley_iii.iges` 纳入 provenance+manifest 或删除；STEP 分支接线或显式 NOT_SUPPORTED。