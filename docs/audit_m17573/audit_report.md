# 交差核验报告 m17517 — GAP-CAD-A1-NURBS-INTERSECT（NURBS 曲面求交 + 裁剪）

审计单 m17573 · 审计席 board-audit · 对象单 m17517（原单 m17346 的返修轮 1）· 仓 TensorCAD · 被审分支 `agent/m17517 @ ae0b100` · 审计工作树 `m17573 @ 338625f`（main）· 元平台 `loongarch64` / Python 3.11.6

判据来源：m17517/m17346 任务目标文本（球×平面交线半径误差 < 1e-8·R；球×球两例解析对照；与 OpenCASCADE/pythonocc 对拍 ≥3 例；退化相切 case 如实报 limitation）。**基准裁判表 `tensor_benchmark_cases` 对本锚点为 0 行**（见 §6.2），因此所有权重都落在目标文本 + 本次审计的独立闭式复核上，报告好看不算数。

## 1. 案例介绍

被测能力：两个参数曲面（NURBS patch 或解析 patch 适配器）求交，返回交线顶点、两曲面上的 uv 参数、闭合判定，并给出两曲面的裁剪域（trim loop / 面积）。几何与参数全部取自被测算例本身：

| 项 | 值 |
|---|---|
| 球 R | 2.0（×平面）/ 1.0 与 1.0、1.2（球×球）|
| 平面 z0 | 0.5 / 相切 case z=1.0（R=1 球顶）|
| 球心距 d | 1.2 / 1.0 |
| 柱 | NURBS 圆柱 R=1, h=4, base_z=-2，绕轴周期 u |
| 斜面 | 法向绕 x 倾 30°，patch 半跨 4.0 m |
| 采样 | samples=(64,16) 单元 + refine_levels=2（默认）|
| 求交内核 | marching-squares 播种 + 4 维 Gauss-Newton 投影精化（`tensorcad/nurbs_intersect.py:1063`）|

## 2. 怎么做的（核验方法，非重做）

1. **取证**：读 m17517 的 `runs/nurbs_intersect_run.json`、`artifacts/nurbs_intersect/*`（6 个工件）与工作树 `work/m17517`；`git status` 显示这些交付物是 **未跟踪（??）** 且 `git check-ignore` 退出 1（非忽略），即分支上没有它们 → 审计把工件按字节快照进 `artifacts/audit_m17573/m17517/`，sha256 记 `artifacts/audit_m17573/manifest.json`。
2. **从干净源码重跑**：本审计工作树是 main `338625f`（不含上述未跟踪文件），用 `examples/nurbs_intersect/run_example.py` 重新生成全部结果到 `artifacts/audit_m17573/rerun/`，与主体 JSON 逐叶比对（脚本 `artifacts/audit_m17573/verify_and_plot.py`）。
3. **独立闭式复核**：审计脚本 `artifacts/audit_m17573/independent_check.py` 在同进程自己算参考值（√(R²−h²)、(d²+R1²−R2²)/2d、2πr、椭圆周长用 Simpson 数值积分而非被测代码里的 AGM），不接受被测 JSON 里的任何 ref 值，并测每个顶点到两曲面的残差、平面外偏差、相切/重合退化行为。
4. **判级**：按 maturity_gate 的 JUDGE_SCHEME 逐键核 evidence，并离线跑了一次 `judge()` 干跑。
5. **测试**：`python3 -m pytest tests/test_nurbs_intersect.py -q` 在审计树重跑。

## 3. 计算结果

### 3.1 半径（审计席独立闭式，容差 1e-8·R）

| 算例 | ref 半径 | got 半径 | 误差/R | 判据 | 判定 |
|---|---|---|---|---|---|
| sphere R=2 x plane z=0.5 (circle) | 1.936491673 m | 1.936491673..1.936491673 m | 5.551e-16 | 1e-8 R | PASS |
| sphere x sphere R=1,1 d=1.2 | 0.800000000 m | 0.800000000..0.800000000 m | 4.441e-16 | 1e-8 R | PASS |
| sphere x sphere R=1,1.2 d=1.0 | 0.960000000 m | 0.960000000..0.960000000 m | 5.551e-16 | 1e-8 R | PASS |

### 3.2 弧长 / 周长（生态系统精度线 < 3%）

| 算例 | ref | got | 相对误差 |
|---|---|---|---|
| sphere x plane arc length | 12.167336028 m | 12.167258038 m | 6.410e-06 |
| sphere x sphere R=1,1 d=1.2 | 5.026548246 m | 5.026516027 m | 6.410e-06 |
| sphere x sphere R=1,1.2 d=1.0 | 6.031857895 m | 6.031819232 m | 6.410e-06 |
| cylinder R=1 x plane 30 deg (ellipse perimeter, Simpson) | 6.777917613 m | 6.777877630 m | 5.899e-06 |

### 3.3 收敛与二阶外推（主体 results.json `/studies/convergence`，审计树重跑同值）

| u 单元 | 交点数 | 弧长 | 相对误差 | 测得阶 |
|---|---|---|---|---|
| 32 | 33 | 6.081174137 | 4.099e-04 | 1.999 |
| 64 | 65 | 6.083044188 | 1.025e-04 | 2.000 |
| 128 | 129 | 6.083512034 | 2.564e-05 | 2.000 |
| 256 | 257 | 6.083629018 | 6.410e-06 | - |

### 3.4 与另一独立算法的对拍（`tensorcad.surface_intersection` 隐式水平集 marcher）

| 算例 | 本算法长度 | 隐式 marcher 长度 | 相对差 |
|---|---|---|---|
| sphere R=2 x plane z=0.5 | 12.167258038 | 12.167258035 | 2.635e-10 |
| sphere R=2 x sphere R=1.2 at z=1 | 5.729732211 | 5.729732210 | 2.761e-10 |
| NURBS cylinder x plane z=0.3 | 6.283145033 | 6.283145032 | 2.493e-10 |

顶点到两曲面的最大残差：球×平面 1.110e-15 m，球×球 6.661e-16 / 1.110e-15 m；平面外偏差 6.661e-16 m（patch 尺度 12.167 m 量级的环）。

### 3.5 诚实失败面（本轮审计实测，非读取）

- 相切（R=1 球顶与 z=1 平面）：曲线 0 条，`tangential_candidate=True` → 如实标记而不编造退化曲线，**PASS**。
- 重合两 patch：抛 `NURBSIntersectionError`（"the two patches coincide over 100% of the sampled domain: the intersection is a surface, not a curve"）而不是返回点云，**PASS**。
- 小 patch 截断（审计席边界扫描，平面半跨 0.6→4.0 m）：整环周长 6.777918 m；`closed=True` 只在半跨 ≥1.2 m（长度比 0.999994）出现，半跨 ≤1.1 m 时一律报 open（长度比 0.07–0.80），半跨 0.6 m 时返回 0 条曲线而不是编一条环；「closed 却短」的违例数 = 0 → **PASS**。

| 平面半跨 (m) | 分支数 | closed | 总弧长 (m) | 占整环比例 |
|---|---|---|---|---|
| 0.6 | 0 | [] | 0.000000 | - |
| 0.8 | 4 | [False, False, False, False] | 0.473427 | 0.069848 |
| 0.9 | 4 | [False, False, False, False] | 1.786369 | 0.263557 |
| 1.0 | 2 | [False, False] | 4.583255 | 0.676204 |
| 1.1 | 2 | [False, False] | 5.416389 | 0.799123 |
| 1.2 | 1 | [True] | 6.777878 | 0.999994 |
| 1.5 | 1 | [True] | 6.777878 | 0.999994 |
| 2.0 | 1 | [True] | 6.777878 | 0.999994 |
| 4.0 | 1 | [True] | 6.777878 | 0.999994 |

## 4. 结果图

- `artifacts/audit_m17573/audit_figure.png`（审计席自绘）：(a) 各算例「got vs 闭式」相对误差与 1e-8·R / 3% 两条判据线；(b) 球×平面交线顶点与解析圆重叠，标题带实测最大半径误差。
- 主体交付的三张图经审计树重跑后 **逐字节相同**（sha256 见 `verify_repro.json`）：`curves_vs_analytic.png`、`convergence.png`、`trim_and_radius_error.png`

## 5. 注意点（坑与局限）

1. **审计第一轮把柱×斜面算例做砸了**：审计自建平面 patch 只给 du=4 m（< 椭圆周长 4.30 m），得到开放曲线 2.764 m，误差 35.8% —— 这是**审计者自己的构造错误**（patch 太小截断曲线），不是被测实现的缺陷；把 patch 半跨放到 4.0 m 后交线闭合、误差回到 1e-5 量级。记下来防止下一位审计者误判。
2. **交付物未入库**：m17517 的 artifacts/、runs/、repro 脚本在工作树里是未跟踪文件，分支推上去也不会带；本轮靠审计快照 + 从 main 源码重建兜住（重建后 3/3 PNG 逐字节一致、JSON 2307/2309 叶相同）。
3. **无基准裁判行**：`tensor_benchmark_cases` 对本锚点 0 行，容差只存在于任务目标文本；这既挡住 L4 认证，也让「误差是否达标」缺少板级登记。
4. **OpenCASCADE 对拍缺失**：目标里写明的 ≥3 例对拍是 WITHHELD —— 本机 OCP/OCC/cadquery 三个绑定都 ImportError（审计独立探测）；`tools/xcheck_occt_nurbs_intersect.py` 写了但从未产出数字。主体没有拿「in-repo 隐式 marcher 对拍」冒充 OCCT 对拍，这点是诚实的。
5. **适用范围**：只对 C0/C1 参数 patch、交线非相切且不重合的常规 case 成立；相切/重合只能报 flag/异常，不能给退化交线；一个 patch 内的多分支由 marching 播种决定，极端拓扑（交线穿过 patch 角点）未验证。

## 6. 误差判定（逐项，显式回答 <3%）

### 6.1 关键量 vs 基准

| 关键量 | 基准来源 | 相对误差 | <3% ? |
|---|---|---|---|
| 球×平面交线半径 | 闭式 √(R²−h²)（审计自算） | 5.551e-16 | 是（≈机精度，且远优 1e-8 判据）|
| 球×球半径（等半径）| 闭式 | 4.441e-16 | 是 |
| 球×球半径（不等半径）| 闭式 | 5.551e-16 | 是 |
| 球×平面弧长 | 2πr（审计自算） | 6.410e-06 | 是 |
| 球×球弧长（两例） | 2πr | 6.410e-06 | 是 |
| 柱×斜面椭圆周长 | Simpson 数值积分（审计自算，积分自检圆周长/2π = 1.000000000） | 5.899e-06 | 是 |
| 与隐式 marcher 对拍 | 仓内独立算法 | 2.761e-10 | 是 |
| OpenCASCADE 对拍 | 厂商内核 | **WITHHELD** | **无法判定**（绑定不可用，如实记）|

结论：所有可判定的解析基准项都 **< 3%**（实际在 1e-5 以下，半径项 ≈ 机精度）；唯一无法判定的是厂商对拍项，如实标 WITHHELD，不粉饰。

### 6.2 基准裁判表原样回执

```sql
SELECT case_name, metric, tolerance FROM tensor_benchmark_cases
 WHERE anchor_slug='GAP-CAD-A1-NURBS-INTERSECT';  -- 0 rows
```
同表 TensorCAD 其它锚点有 2 行（`GAP-CAD-T1-BASELINE-FLEET`、`GAP-CAD-B2-MESH-FORMAT-MATRIX`），说明表是活的，本锚点确实没登记。

## 7. 判级与裁定

| 级 | 判据（maturity_gate JUDGE_SCHEME）| 本轮核验结果 |
|---|---|---|
| L1 | tests | 18 passed in 73.67s（审计树重跑；主体声称 18/18 一致）→ **过** |
| L2 | impl_files + entry + run_artifacts | `tensorcad/nurbs_intersect.py` 1238 行，placeholder 标记 0 个；入口 `nurbs_surface_intersection`（包级导出 + `examples/nurbs_intersect/run_example.py`）；工件 6 件真实非空（PNG 79–262 kB，results.json 101 kB）+ 审计重建工件 → **过** |
| L3 | + integration_callers | **不过**：全仓 grep 只命中定义/导出/examples/tools/tests，`tensorcad/cli.py` 无该子命令，无 adapter/主求解器调用 |
| L4 | + benchmark + repro_cmd | **不过**：`tensor_benchmark_cases` 0 行（无板级容差）；OCCT 对拍 WITHHELD。注：数值质量本身达到 L4 水准（误差 ~1e-16 R、二阶收敛、可复现）|
| L5 | + external_anchor | 不过：无厂商/实验外部锚 |

**审计判级：L2**（离线 `maturity_gate.judge()` 干跑：claimed 2 → gate_level 2，`why=[]`）。verdict 已落表 id=590（`tensor_maturity_verdicts`, repo=TensorCAD, scope='m17517 GAP-CAD-A1-NURBS-INTERSECT', verdict_by='board-audit'），等 maturity_gate 周期裁定（confirmed=消账）。

### 发现（按严重度）

| id | 严重度 | 发现 | 证据 |
|---|---|---|---|
| F1 | medium | 交付物未跟踪，分支不携带 | `git status --short` in work/m17517；`git check-ignore` 退出 1 |
| F2 | low | 文档渲染于 f43e888，run json 记 ae0b100；数值一致 | `nurbs_intersect_run.json:/revision` vs `docs/nurbs_intersect.md:7` |
| F3 | medium | 目标里的 OCCT ≥3 例对拍未完成（WITHHELD）| `occt_xcheck.json`；审计独立 import 探测 |
| F4 | low | 锚点无基准裁判行 | 上述 SQL 0 行 |

**处置：approve with notes**（L2 confirmed；建议另开接线任务给 L3/L4 补 integration caller + 基准登记行，并在有 OCP 的镜像上补 OCCT 对拍）。

---

审计工件：`artifacts/audit_m17573/`（independent_check.py/.json、verify_and_plot.py、verify_repro.json、manifest.json、audit_figure.png、m17517/ 快照、rerun/ 重建、logs/ 三份日志、insert_verdict.py）。本页由 `artifacts/audit_m17573/render_report.py` 程序化渲染，数字均来自上述 JSON，零手抄。
