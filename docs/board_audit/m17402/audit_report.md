# m17350 交差审计：GAP-CAD-B2-MESH-FORMAT-MATRIX（TensorCAD meshio 全格式矩阵）

**审计单**: m17402（board-audit） · **被审单**: m17350（succeeded, auto-approve）
**审计人**: board-audit · **日期**: 2026-11-13 · **基准裁判**: tensor_benchmark_cases（PostgreSQL tensorplatform）

---

## ① 案例介绍（物理背景/几何/参数）

m17350 属 owner 2026-11-13 TensorCAD 专项攻关批 B 战线（互操作）：验证 TensorCAD 表面三角网格
经 meshio 5.3.5 在 6 格式（STL/OBJ/VTK/SU2/CGNS/Gmsh）× 3 几何（cube/sphere/wigley）上的
**往返保真**（体积/拓扑/取向不变）。验收标准原文：「往返保真：体积/拓扑/取向不变（6 格式×3 几何=18 组合全绿表）」。

| 几何 | 顶点/面 | 报告标称解析体积 | 实测有符号体积 |
|---|---|---|---|
| cube（单位立方体，外向） | 8v/12f | 1.0 m³ | +1.000000 |
| sphere（16×32 经纬球） | 513v/992f | 4π/3≈4.18879 | −4.121942（法向朝内，符号一致） |
| wigley（Wigley III 半船体，L=2,B=0.3,D=0.125，无甲板开面） | 341v/678f | ≈0.075 | −0.034432（开面，非封闭体） |

## ② 怎么做的（方法/流程/工具链）

本单为**交差核验**（非重做），独立取证路径：
1. 读 m17350 工作树 `/root/tp-agent/work/m17350`（branch agent/m17350）：实现
   `tensorcad/meshio_exchange.py`（182 行）、测试 `tests/test_meshio_format_matrix.py`（324 行）、
   报告 `docs/meshio_format_matrix_report.md`、`runs/m17350_meshio_matrix/{run_artifacts.json,repro_cmd.sh}`、
   2 张 PNG 工件；
2. 基准裁判查询：`SELECT case_name,metric,tolerance FROM tensor_benchmark_cases WHERE
   anchor_slug='GAP-CAD-B2-MESH-FORMAT-MATRIX'` → **0 行**（全表 53 行，唯一 CAD 锚点为
   GAP-CAD-T1-BASELINE-FLEET id=52）；
3. 实跑 repro：`python -m pytest tests/test_meshio_format_matrix.py -v` → **53 passed, 0 failed, 1.29s**（rc=0）；
4. 独立复算（不复用被测代码的断言，直接调 `round_trip_meshio` + 散度定理体积）：12 组合
   体积相对误差 **0.000e+00**（<1e-15），拓扑 8/8·12/12、513/513·992/992、341/341·678/678 全保；
   OBJ/VTK/Gmsh 逐点逐索引 exact=True；STL 重排（exact=False，体积/取向仍保）；SU2/CGNS 6/6
   正确抛 `MeshioFormatError`；
5. 平台报告核验：tensor_reports id=16345，六段（## 1–## 6）正则全命中，figures=2；
6. git 取证：`git status --short` 显示 5 个交付文件全部 `??`（untracked）；`git log agent/m17350`
   head=505d10d（=main）；`git ls-remote origin | grep 17350` → 无远端分支。

## ③ 计算结果（数字表，程序化注入零手抄）

| 项 | 报告声称 | 本单实测 | 判定 |
|---|---|---|---|
| 测试通过 | 53/53 | 53/53（1.29s, rc=0） | ✅ 一致 |
| cube 体积（4 格式） | 1.000000 vs 1.000000 | 1.000000 vs 1.000000, err=0.000e+00 | ✅ 一致 |
| sphere 体积（4 格式） | −4.121942 vs −4.121942 | −4.121942 vs −4.121942, err=0.000e+00 | ✅ 一致 |
| wigley 体积（4 格式） | −0.034432 vs −0.034432 | −0.034432 vs −0.034432, err=0.000e+00 | ✅ 一致 |
| 拓扑（12 组合） | 计数不变 | 计数不变（上表） | ✅ 一致 |
| 无损格式 exact | OBJ/VTK/Gmsh 逐点逐索引 | exact=True（atol 1e-12 + 索引全等） | ✅ 一致 |
| SU2/CGNS 拒绝 | 正确拒绝 | 6/6 抛 MeshioFormatError | ✅ 一致 |
| 工件 sha256 | — | run_artifacts.json=29dce299…；repro_cmd.sh=b91b6b56…；volume_error_comparison.png=bdb87c90…；topology_preservation.png=e149e3eb…；meshio_exchange.py=f8d8517c…；test_meshio_format_matrix.py=20a0a2d7… | 已钉 |

## ④ 结果图

被审工件 2 张 PNG（2250×600，有效 PNG 数据）随本报告内嵌：
- `volume_error_comparison.png`：12 组合体积误差对比（全部 <1e-15，远低于 0.1% 判据线）；
- `topology_preservation.png`：往返前后顶点/面计数对比。
（审计单自身无新计算场量，图=被审工件原图，非本单生成。）

## ⑤ 注意点（坑/局限/适用范围）

- **F1（高）交付物全部 untracked**：5 个交付文件（模块/测试/报告/runs/工件）在 m17350 工作树
  均为 `??` 未跟踪；branch agent/m17350 与 main 同 head（505d10d），无 commit、无远端分支。
  任务 result 自证：`branch: null, head_commit: null`，`server_quality_audit.gate_valid=false`，
  靠 TP_AUTO_APPROVE_DEV_ACCEPTANCE 自动放行。工作树一清，能力即蒸发。
- **F2（高）基准裁判缺席**：tensor_benchmark_cases 无本锚点行 → 平台机器裁判无从裁定，
  L4 不可判过（唯一裁判缺席=不能判过）。
- **F3（中）平台报告图未内嵌**：tensor_reports id=16345 两条 figure 的 png_base64 为空
  （仅 name/caption），PNG 仅存工作树。
- **F4（低）报告自相矛盾**：§6 写「4/4 unsupported combinations」，矩阵实为 6 组合
  （SU2/CGNS×3 几何）；§1 表 wigley「Analytic Volume ≈0.075」与实测 −0.034432 并存
  （开面半船体，散度定理值非船体排水体积），易误读。
- **F5（信息）验收口径**：任务要「18 组合全绿」，实际 12 PASS + 6 NOT_SUPPORTED
  （meshio SU2/CGNS writer 仅体单元，属格式本质限制，非 TensorCAD 缺陷，文档化处理诚实）。
- 适用范围：本模块仅表面三角网格；SU2/CGNS 需先四面体化（模块已给出明确报错指引）。

## ⑥ 误差判定（逐项 + 显式回答）

| 关键量 | 判据 | 实测 | 是否 <3% |
|---|---|---|---|
| 体积相对误差（12 组合） | <0.1%（任务判据） | 0.000e+00（<1e-15） | **是**（远小于 3%） |
| 拓扑计数（12 组合） | 全等 | 全等 | **是**（0 偏差） |
| 取向（12 组合） | 符号不变 | 符号不变 | **是** |
| 无损格式逐点/逐索引 | 全等 | exact=True | **是** |
| 外部参考值误差 | 无外部基准（自指往返保真，无 ref 值可对标） | 不适用 | 无基准，原因=往返保真类任务无外部参考解 |

**结论**：功能实质=真代码独立功能（输入网格→格式转换→输出网格，53 项不变量断言全绿，
本单独立复算逐位一致）→ **L2**。不达 L3：全仓 grep 无任何集成调用者（0 命中）。
不达 L4：基准案号未注册（F2）+ 交付物未入库（F1）+ 无外部参考值。
**verdict = L2**；F1/F2 为返修项（commit 入库 + 注册基准案号），gate 机器复核。
