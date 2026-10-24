# TensorCAD 项目分析（0.3.0）

- 分析日期：2026-10-24
- 分析对象：`cxs503/TensorCAD`，`main` @ `4c37665`（0.3.0，pyproject / `__version__` / CHANGELOG 三处一致）
- 分析方式：静态阅读全部源码与文档 + 本地全量构建/测试/lint/mypy/示例 + 独立数值复核 + 边界与失败路径实测。**未修改任何源码。**
- 环境：CPython 3.11.2，`pip install -e ".[test,dev]"`（numpy / pytest / ruff / mypy / pytest-cov）。**未安装 torch**，因此神经几何的*训练*路径未重跑；但抽取器（numpy-only）、bundle 校验、CLI `info`、以及真实权重包的结构全部实测过，见 §4。
- 上一份报告 `docs/analysis/2026-10-24-project-analysis.md` 针对 0.2.0：其中 **1 项已修复**（覆盖率门禁真正生效）、**1 项部分修复**（神经几何线已并进 README 与格式文档，schema 0.2 已对齐加载器）、**7 项仍成立**。本文 §5 逐条标注了状态。

**一句话结论：** 这仍是一个工程纪律远高于同阶段研究原型的几何基础库——"承诺即测试"在 0.3.0 扩展到了事务与元数据层；但 **0.3.0 新引入的神经几何线有 4 处真实缺陷（其中 2 处会造成活动版本被污染或派生件被静默删除）**，0.2.0 遗留的 `estimate_discretisation_error` 越界崩溃依旧可复现，并且**首次被测量的通用网格静水力性能问题（O(V·S)）恰好卡在两条线的汇合点上**。

### 问题速览

| 编号 | 位置 | 后果 | 状态 |
|---|---|---|---|
| P1-1 | `hull.py:405-412` | 收敛估计越界崩溃；`tensorcad-hull` 留下无 metadata 的半成品目录 | 0.2.0 遗留，已复现 |
| P1-2 | `cli.py:72-80` | 导出/元数据阶段的异常直接变成 traceback，`--json` 无结构化错误 | 0.2.0 遗留，已复现 |
| P1-3 | `cli_neural.py:280` | 无 torch 时 `verify` 整体崩溃（前 4 项检查本不需要 torch）；且从不校验 `derived` | 0.3.0 新增，已复现 |
| P1-4 | `hydrostatics.py:194/211` | `max_section_area` 近似 O(V²)；64³ 神经网格静水力 14.1 s，96³ 约一分钟 | 0.3.0 首次测量 |
| P2-1 | `objects.py:529-543` | 发布前把派生网格写进**活动** bundle：失败提交会污染活动版本 | 0.3.0 新增 |
| P2-2 | `objects.py:671` | `extract` 不登记 `derived[]`；下一次 `edit` 会静默删除这些文件 | 0.3.0 新增 |
| P2-3 | `neural_geometry.py:353` | 读端不校验 `derived` 记录，`../` 路径逃逸被接受 | 0.3.0 新增，已复现 |
| P2-4 | `hull.py:generate_hull` | 网格粗到 −35.7 % 体积仍 `warnings == []` | 0.2.0 遗留，已复现 |
| P2-5 | README / CHANGELOG | 测试数、覆盖率数字与实际不符 | 已复现 |
| P3 ×10 | 见 §5 | py.typed、OBJ 语义丢失、水线语义、重复计算、`write_stl` 无测试…… | 多为遗留 |

---

## 1. 项目定位

TensorCAD 不是桌面 CAD 的复刻，而是 Tensor 科研生态里的**几何前置层**：

| 层 | 角色 | 状态 |
|---|---|---|
| TensorCAD | 参数化几何、验证、网格、元数据与版本化 | 已实现（0.3.0） |
| TensorShipDesign | 基于 TensorCAD 的船舶设计应用 | 计划中 |
| TensorLBM / FEM / DEM / BEM | 消费几何与网格的求解器 | 外部 |
| TensorPlatform | 编排 agent 任务、仿真、版本与迭代 | 外部 |

0.3.0 的核心主张值得认真对待：**几何的权威表示不应该是网格文件，而应该是"连续场 + 可编辑状态 + 版本/事务元数据"，OBJ/STL 只是派生件。** 这是本项目与普通参数化建模脚本最本质的区别，也是它把"可复现"从一句口号变成一个可被哈希校验的目录结构的唯一办法。

---

## 2. 代码结构与依赖分层

```
tensorcad/
  __init__.py          67   公共门面（__all__ 23 项）
  geometry.py         183   Mesh / BoundingBox / MeshBuilder（精确坐标共享顶点）
  parameters.py       122   ParameterSpec / RatioSpec，一次报全部错误
  mesh.py             297   边统计、水密、定向、焊接、退化面、积分
  validation.py       292   Check / ValidationReport / inspect_geometry / 尺寸比对
  hydrostatics.py     215   散度定理面元积分：体积、湿面积、水线面、LCB/LCF/VCB、Cb/Cp
  export.py           206   OBJ/CSV/JSON 写出 + OBJ 读回 + 二进制 STL
  hull.py             436   唯一的船型专用模型
  cli.py               97   tensorcad-hull
  neural_geometry.py  722   ★ NGO 0.2 格式：清单、哈希、原子发布、编辑事务（纯 stdlib）
  neural/extract.py   371   marching tetrahedra → 闭合网格（仅需 numpy）
  neural/field.py     525   ConditionalTSDF 训练/度量（需 torch）
  neural/objects.py   719   fit / edit / extract 的对象级编排（需 torch）
  cli_neural.py       349   tensorcad-neural：fit/edit/extract/info/verify
  adapters/            16   OCP/gmsh 接缝（空包 + 边界文档）
tests/              1,900+  13 个文件
examples/                    4 个可执行基准（自带断言，失败即非零退出）
docs/                       格式说明 + 3 份研究文档 + 2 次真实实验的完整产物
schemas/                     manifest 0.1 / 0.2 JSON Schema
```

源码 **4,633 行**（核心 ~1,930，神经几何线 ~2,700）/ 测试 **1,942 行**。**依赖方向是单向且可执行的**：

```
hull ──▶ mesh ──▶ geometry          （纯 stdlib）
  │        └──▶ validation
  ├──▶ hydrostatics ──▶ mesh
  ├──▶ parameters
  └──▶ export ◀── neural.objects ──▶ neural.field (torch) / neural.extract (numpy)
                     │
                     └──▶ neural_geometry（纯 stdlib：清单/哈希/事务）
```

三条独立承诺都有测试兜底：`import tensorcad` 不引入 numpy/torch（`test_public_surface.py` 重新 import 后断言）、核心不 import `adapters`、`neural_geometry` 不反序列化可执行对象。

---

## 3. 实测状态

### 3.1 门禁

| 项目 | 结果 | 说明 |
|---|---|---|
| `pytest`（默认，无 torch） | **175 passed，2 module skipped**，3.7 s | 跳过的是 `test_cli_neural.py` / `test_neural_pipeline.py`（需要 `neural` extra） |
| `pytest`（含 neural extra，CI） | 175 + 14 = **189** | 与 CHANGELOG 的 "189 tests" 一致 |
| 覆盖率 | **94.47 %**（阈值 90 %，已真实生效） | 0.2.0 报告指出的"覆盖率门禁不存在"**已修复** |
| `ruff check .` | All checks passed | |
| `mypy`（strict） | Success，16 个源文件 | |
| 示例基准 | 全部检查 PASS，退出码 0 | 含五级网格收敛研究与解析体积对照 |

各模块覆盖率：`geometry/hull/hydrostatics/mesh/parameters/validation/cli` 均 98–100 %；`neural_geometry` 89 %、`neural/extract` 93 %、**`export.py` 84 %**（原因见 P3-6）。

### 3.2 船体数值基线（默认 L=10, B=2, T=0.8）—— 与文档逐位一致

```
41×13   1001 顶点 /  1998 面   体积 7.446554（解析 7.512698，−0.880 %）
161×41 12961 顶点 / 25918 面   体积 7.502488（−0.136 %）
401×101 80401 顶点/160798 面   体积 7.510140（−0.034 %）
LCB/LCF/VCB = 5.0000 / 5.0000 / −0.3170 ; Cb/Cp = 0.4654 / 0.7031
```

误差随加密单调下降，示例脚本自己会断言这一点——**这是仓库里最有价值的工程习惯：数字不是写进 README 的，是算出来并被断言住的。**

### 3.3 性能（纯 Python，本机 2 vCPU）

| 网格 | 顶点/面 | 生成 + 验证 + 静水力 | 峰值 RSS |
|---|---|---|---|
| 41 × 13 | 1,001 / 1,998 | 36 ms | — |
| 161 × 41 | 12,961 / 25,918 | 617 ms | — |
| 401 × 101 | 80,401 / 160,798 | 4.15 s | — |
| 801 × 201 | 320,801 / 641,598 | **25.2 s** | **919 MB** |

按此线性外推，参数上限 `4001 × 513`（约 410 万顶点 / 820 万面）需要 **约 5 分钟、10 GB 以上内存**——即"参数范围允许，实现上不可用"。对 agent 调用的默认分辨率（41×13，36 ms）完全够用；问题只出现在上限区间，建议把文档里的 `nx/nz` 上限改成"推荐上限 + 硬上限"两档，或在超过 801×201 时给出显式 warning。

---

## 4. 独立复核（不信任文档，自己算）

1. **仓库里 3 个真实权重包全部通过 `load_bundle()` 校验**（SHA-256、元数据 JSON、路径安全），0.1.0 旧包确实仍可解析——"升级不打开旧文件"的承诺成立。
2. **0.3.0 的参考实验产物自洽**：`summary.json` 里的 96³ 抽取（97,332 顶点 / 194,660 面、水密、V−E+F=2、体积 0.5154 m³、面积 3.1756 m²）与 manifest 的 `derived[]`、`validation.json` 逐项一致；`history.json` 三条记录（committed / committed / rolled_back），失败项 `inf/nan` 已序列化为 `null`。
3. **抽取器用纯 numpy 独立复测**（解析球 `|x|−0.5`，bounds ±0.58）：96³ 用时 **1.28 s**，96,614 顶点 / 193,224 面，**水密、边流形、V−E+F=2**，体积误差 **−0.029 %**、面积误差 **−0.015 %**。这个模块的实现质量很高（Kuhn 六四面体划分 + 逐三角形局部定向，无全局修补）。
4. **两条线在验收层其实已经兼容**：把抽取出的网格喂给 0.2.0 的 `inspect_geometry()`，10 项检查全 PASS（含顶点流形 warning 项）。也就是说"神经场 → 网格 → 复用船体验收管线"只差一个显式调用，**但性能上目前过不去**——见 P1-4。
5. `tensorcad-hull` 两次同参数字节级一致（`hull.obj`、`metadata.json` 均 `cmp` 通过，元数据里没有时间戳）。

---

## 5. 问题清单

严重度：**P1 = 会产生错误结果、崩溃或污染数据**；**P2 = 承诺与实现不一致、可被合法输入触发**；**P3 = 整洁性/健壮性**。P1-1…P1-4 与 P2-3/P2-4/P2-5 都在本机跑出了证据；P2-1/P2-2 需要 torch 才能走完训练路径，故由代码路径直接判定（并给出相应的复核方法）。

### P1-1（0.2.0 遗留，仍成立）`estimate_discretisation_error` 在加密越界时抛异常，且留下半成品导出

`hull.py:405-412` 的守卫只检查 `nx*nz > 400_000`，但加密后是 `2*nx-1 × 2*nz-1`，会越过 `nx ≤ 4001 / nz ≤ 513`：

```bash
$ tensorcad-hull --nx 4001 --nz 3 --no-export --json
ValueError: nx must be <= 4001 stations, got 8001      # exit 1，完整 traceback

$ tensorcad-hull --nx 2002 --nz 3 --output /tmp/partial
ValueError: nx must be <= 4001 stations, got 4003      # exit 1
$ ls /tmp/partial
hull.obj  offsets.csv                                  # ← metadata.json 缺失：半成品导出
```

三个独立触发点（`2002×3`、`3×258`、`4001×3`）都实测崩溃，而**模型本身完全合法且已生成成功**——崩溃只发生在"收敛估计"这一步。更糟的是 `export_hull` 先写 OBJ/CSV、最后写 metadata，所以失败会留下一个**没有 metadata 的资产目录**，下游无法判断它对应什么参数。

现有测试 `test_hull.py:137` 用的是 `nx=1500, nz=270`（`nx*nz = 405,000`），恰好被体积守卫拦住，所以这个边界完全没有覆盖。

修复（约 8 行，同时修掉丢 `water_density` 的隐式默认）：

```python
refined_nx, refined_nz = 2 * p.nx - 1, 2 * p.nz - 1
if (p.nx * p.nz > MAX_MESH_SIZE_FOR_CONVERGENCE
        or refined_nx > MAX_LONGITUDINAL_STATIONS
        or refined_nz > MAX_VERTICAL_LEVELS):
    return {"status": "skipped",
            "reason": f"refined grid {refined_nx}x{refined_nz} exceeds the "
                      f"{MAX_LONGITUDINAL_STATIONS}x{MAX_VERTICAL_LEVELS} parameter limits"}
refined = HullParameters(p.length, p.beam, p.draft, refined_nx, refined_nz,
                         p.water_density)          # ← 第 405 行丢掉了它
```

回归测试建议参数化 `(4001,3) / (3,513) / (2002,3) / (1500,270)`。

### P1-2（0.2.0 遗留）`tensorcad-hull` 的 try 只包住 `generate_hull`

`cli.py:72-80`：`generate_hull` 失败返回 2，但 `export_hull(...)` 与 `metadata(with_discretisation=True)` 的异常全部冒泡成 traceback（P1-1 就是现成例子），`--json` 模式下也没有结构化错误。对"agent-callable"的 CLI 来说，这是最关键的一条错误契约。

```python
try:
    ...
except ValueError as exc:
    print(f"tensorcad-hull: invalid model: {exc}", file=sys.stderr)
    if args.json:
        print(json.dumps({"error": str(exc), "stage": stage}, sort_keys=True))
    return 2
```

顺带建议：先算 metadata（含收敛字段）再落盘，避免像 P1-1 那样留下半成品目录。

### P1-3 `tensorcad-neural verify` 在没有 torch 时直接崩溃

```bash
$ tensorcad-neural verify --object docs/experiments/.../bundle
Traceback ... from .neural.objects import build_field
ModuleNotFoundError: No module named 'torch'          # exit 1，完整 traceback
```

`cli_neural.py:280` 的 `from .neural.objects import build_field` 位于整个 `try` **之外**。而 verify 的前 4 项检查（清单解析、revision、权重哈希、history 对账）**根本不需要 torch**；只有 `field_query` 需要。这与"元数据层纯 stdlib"的分层承诺直接冲突，也让 CI 里那条"torch-free 检查"只覆盖到 `--help`。

建议：把 import 移进 `try`，`field_query` 在缺少 extra 时记录为 `"skipped: requires the neural extra"` 而不是让整条命令无法运行；同时**补上 `derived` 完整性检查**（P2-3）。

### P1-4 `max_section_area` 是 O(V × 站点数)，使"神经网格 + 静水力"实际上不可用

`hydrostatics.py:211-215` 对每个 `x` 站点调用一次 `section_area_at`，后者又遍历**全部顶点**。船体网格的 x 只有 `nx` 个取值，所以问题被掩盖；**marching tetrahedra 的网格每个顶点 x 几乎都不同**，于是退化为近似 O(V²)：

```
48³ 抽取网格：24,182 顶点 / 48,360 面 / 4,665 个不同 x
compute_hydrostatics 4.40 s，其中 max_section_area 4.23 s（96 %）
64³ 抽取网格：            85,536 面 / 8,163 个不同 x
compute_hydrostatics 14.09 s，其中 max_section_area 13.79 s（98 %）
对照：161×41 船体网格 25,918 面，同一函数 0.18 s
```

按 `V × 站点数` 外推，96³（约 19 万面、~1.9 万个不同 x）需要**一分钟量级**，而这正是 §4 第 4 点里"神经场产出 → 复用船体静水力验收"的必经之路。

修法很直接：一次遍历按 `x` 分桶（`dict[float, list[float]]`），把每站的 `(z, |y|)` 收集好再积分——O(V log V)，顺带让 §4 第 4 点的合流真正跑得动。

### P2-1 `apply_dent_edit` 在**发布之前**就往活动目录里写派生网格

`objects.py:529-543`：编辑通过门禁后，`_extract_and_write(directory=store.root / "derived", ...)` 直接把网格写进**当前活动版本目录**，随后才 `store.publish()`（内部再 stage → commit）。这条路径破坏了 `ObjectStore` 自己声明的两条不变量：

- 抽取本身是秒级的，而这段时间活动版本里已经出现 manifest 未记录的 `derived/`；
- 一旦 `commit` 失败（`expected_parent` 过期、staged manifest revision 不符、staging 写入失败），**活动版本不再是"逐字节不变"**——多了一套不属于任何已发布 revision 的网格。

`examples/neural_geometry/run_object_pipeline.py` 打印的 "manifest bytes unchanged" 只覆盖"门禁失败"分支（该分支不做抽取），所以 CI 抓不到。

修法：在临时目录里抽取，把字节直接放进 `files` 字典交给 `publish`；`derived[]` 记录用内存里的 `sha256_bytes` 计算，不要在活动目录落盘。

### P2-2 `extract` 不登记派生件，并会被下一次 `edit` 静默删除

`extract_object()`（`objects.py:671-700`）写出 `derived/surface.obj|stl` 并返回记录，但**从不更新 manifest 的 `derived[]`**（只有 `edit --resolution` 那条路径会写）。后果是两个可观测的异常行为：

1. `tensorcad-neural info` 对一个刚 extract 过的对象仍显示 `derived ... none`；
2. `commit` 用目录级 rename 替换整个 bundle 根，**staging 里没有的文件会随之消失**——所以一次 `extract` 之后再跑一次 `edit`（默认 `--resolution 0`），派生网格被**静默删除**。

文档 `docs/neural-geometry-format.md` §5 写的"两者都以哈希记录在 manifest 的 derived[] 里"只对 edit 路径成立。建议：`extract` 通过一次"仅元数据更新"的事务把 `derived[]` 写回（revision 不变或 +1 需明确），或在文档里明确它是一次非权威的旁路导出。

### P2-3 `parse_manifest` 不校验 `derived` 记录，路径逃逸能被接受

写端 `derived_record()` 强制包内相对路径，但读端（`neural_geometry.py:353-356`）只判断"是不是 dict"：

```text
# 把 manifest 的 derived[0].path 改成 ../../../../etc/passwd 后：
load_bundle(...)  →  接受，man.derived[0]["path"] == "../../../../etc/passwd"
```

`info --json` / `summary()` 会把它原样打印出去。这与格式文档"路径必须是包内相对路径，`../` 与绝对路径一律拒绝"矛盾，也和 `load_bundle` 对 `weights/state/training/validation` 的严格校验不对称（`test_neural_store.py:196` 只测了写端）。建议在读端对每条记录复用 `_safe_relative_path` + `^[0-9a-f]{64}$` + `bytes >= 0`，并加一条"读端也拒绝"的测试。

同样的不对称还有：`ObjectStore.history()` 用的是裸 `json.loads`（允许 `NaN`），而其余读取路径都带 `parse_constant=_reject_constant`。

### P2-4（0.2.0 遗留）没有任何"网格太粗"的告警

```python
m = generate_hull(HullParameters(nx=3, nz=3))
m.validation.warnings                              # []      ← 空
m.hydrostatics.displacement_volume                 # 4.8284（解析 7.5127，−35.7 %）
```

`validation.warnings` 是文档指给用户的"下一步看什么"的入口，但它只覆盖比例、偶数列 `nx`、顶点流形。粗网格造成的数值失真**零成本可知**（`TAPER_INTEGRAL` 就在同一个文件里），却完全没有提示。建议在 `generate_hull` 增加一条 warning：网格体积与解析体积的相对误差 > 2–5 % 时提示提高 `nx/nz` 并指向 `estimate_discretisation_error`。

### P2-5 文档数字漂移

| 位置 | 声称 | 实测 |
|---|---|---|
| `README.md:112` | "129 tests" | 默认 175（含 extra 189）——**0.2.0 的陈旧数字** |
| `CHANGELOG.md` 0.3.0 | "currently 96 %" | 94.47 %（`export.py` 只有 84 %） |
| `README.md` 导出章节 | "every export carries units, coordinate frame" | `read_obj` 把坐标系丢弃为 `"unspecified right-handed frame"`（见 P3 第 2 条） |

建议把这三处改成由脚本生成的数字，或在 CI 里加一条"README 里的测试数 == `pytest --collect-only` 输出"的廉价断言。

### P3 级

1. **缺 `tensorcad/py.typed`**（0.2.0 遗留）：`mypy --strict` 只惠及仓库内部，PEP 561 下游拿不到类型信息。
2. **`read_obj` 丢弃语义元数据**（0.2.0 遗留）：顶点逐位往返没问题，但 `# units:` / `# coordinate frame:` 注释不被解析；实测把 `coordinate_frame` 换成任意值，读回来都是默认值（`units` 只是因为默认恰好是 `"m"` 才"看起来"往返成功）。英尺制模型会被静默标成米制。
3. **`compute_hydrostatics(mesh, waterline_z=...)` 的语义可被误用**（0.2.0 遗留）：`waterline_z` 只用于筛选湿面，体积始终是整个闭合网格的积分。实测对默认船体传 `waterline_z=-0.4`：体积纹丝不动（4.8284），湿面积/水线面却变成 12.15 / 7.07 —— 一个"一半正确"的结果，且无任何提示。建议断言 `waterline_z == max(z)`，或真正实现水线切割。
4. **重复实现与重复计算**（0.2.0 遗留）：`mesh._tet_volume`（×6 体积）与 `hydrostatics._tetra_volume`（已 ÷6）是同一公式两份实现；`compute_hydrostatics` 每面依次调用 `face_centroid` / `triangle_area` / `face_normal`，同一个叉积算 2–3 次，湿面再来一次。
5. **`write_stl` 完全没有测试**：`grep -rn write_stl tests/` 只命中 `test_public_surface.py` 的 `__all__` 断言。它是公开 API、二进制格式、给下游查看器/打印用的。我手工验证过它是**正确**的（`84 + 50*n` 字节、声明的三角形数、法向全为单位向量），但"正确且无人守"是最容易腐化的状态。
6. **`export.py` 84 % 覆盖率**就是上一条的直接后果，而 `fail_under = 90` 只对总量生效，单个模块掉下去不会报警。
7. **`expected_parent` 不是原子的**：`commit` 先 `active_revision()` 再 `os.replace`，两个写者可以同时通过检查、后者覆盖前者。文档"两个人不会互相覆盖"的表述比实现强——准确说法是"乐观并发，窗口极小但存在"。真要保证需要 `flock` 或"带版本号的目录 + 唯一符号链接"。
8. **无输入大小上限**：`extract_object` 与 `_weights_bytes` 一次性 `read_bytes()` 整个权重文件；`sha256_file` 已改为分块（好），但解析路径仍是全量读入。对不受信任的 `.ngo` 建议加一个 `max_payload_bytes`。
9. **`schemas/…-v0.1.json` 的 `format_version` 正则** `^0\.[0-9]+\.[0-9]+$` 允许 `0.01.0`，而加载器拒绝前导零。0.2 的 schema 已改对，0.1 的仍与唯一事实来源（加载器）不一致。
10. **`beam_matches_request` 与 `check_consistency` 的"采样亏欠允许带"重叠**（0.2.0 曾指出）：功能正确，但两块逻辑说的是同一件事，长期应收敛成一处并加注释。

---

## 6. 做得好的地方（都有证据）

1. **承诺即测试**：stdlib-only、闭合实体、OBJ 17 位往返、字节级可复现、失败编辑不发布——每一项都有对应测试或 CI 步骤，且我逐项复核过。
2. **失败的方式正确**：参数错误一次报全；几何不可用直接 `raise`；warning/error 分级；示例脚本自己会失败；rolled-back 的编辑退出码是 1 而不是 0。
3. **事务是真的**（`ObjectStore`）：隐藏 staging + 哈希校验 + 目录级 rename；`history.json` 记录每一次尝试（含失败项与实测指标），`inf/nan` 写成 `null` 保证报告本身永远可解析。这个不变量在 `test_neural_store.py`（17 个测试）里被认真覆盖。
4. **数值诚实**：明确区分"网格积分值"与"解析值"，把 0.88 % / 0.14 % 写进元数据；收敛估计标注为"误差估计，不是收敛证明"。
5. **抽取器实现质量高**：Kuhn 六四面体 + 局部定向 + 按整数边键共享交点，避免了 marching cubes 的裂缝与容差焊接；我独立复测的结果（水密、V−E+F=2、体积误差 0.03 %）比文档还好一点。
6. **CI 维度合理**：3.10–3.13 矩阵、ruff + mypy strict + 真覆盖率门禁、CLI 产物、可复现性、基准、以及独立的 30 分钟神经几何实验流水线（带 `--fail-under=85`）。

真正的短板不在工程质量，而在**能力面**：只有一个前后对称、无方尾、无干舷、龙骨带 knuckle 的解析船型；没有 B-rep/STEP（issue #4 仍未落地）；没有面网格质量指标（长宽比/最小角/雅可比）；神经场只有球族；两条线还没有在"产出 → 验收"上真正合流（§4 第 4 点说明只差一步，但 P1-4 卡着）。

---

## 7. 建议路线（按投入产出排序）

**0.3.1（半天到一天，全是小改动、全是真实缺陷）**

1. 修 P1-1：加密越界返回 `skipped` + 显式传 `water_density` + 4 组参数化回归测试。
2. 修 P1-2：CLI 全流程 try + `--json` 结构化错误 + 先算 metadata 再落盘（消灭半成品导出）。
3. 修 P1-3：`verify` 的 import 移进 try、`field_query` 缺 extra 时标注 skipped、补 `derived` 完整性检查。
4. 修 P2-1 / P2-2：派生件只在 staging 内产生；`extract` 登记 `derived[]`（或明确文档化为旁路）。
5. 修 P2-3：`parse_manifest` 校验 `derived[]`；`history()` 改用严格 JSON 解析。
6. 修 P2-4：粗网格数值告警。
7. 加 `tensorcad/py.typed`；同步 README/CHANGELOG 数字；给 `write_stl` 补一条字节级测试。

**0.3.2（性能与数值边界）**

8. 修 P1-4：`max_section_area` 单次分桶；随后把"神经网格 → `inspect_geometry` + `compute_hydrostatics`"写成一条正式的验收路径（`examples/` 里加一个断言版基准）。
9. 合并四面体体积/面几何的重复实现，加一次静态性能剖析基线（避免回归）。
10. 把 `nx/nz` 上限拆成"推荐 / 硬上限"，或在 801×201 以上给 warning。
11. 面网格质量检查（长宽比、最小角、退化阈值、法向一致性）并入 `ValidationReport`——这是"求解器能不能用"的第二道门。

**0.4（能力扩展）**

12. `tensorcad.adapters.ocp`：B-rep 构建、内核有效性、STEP 往返（核心仍 stdlib-only）。
13. 船型 v2：方尾、干舷、前后非对称、真实剖面放样。
14. 神经对象扩展到椭球/圆柱/船体，把"场误差"与"提取网格误差"分开报告（与 `docs/research/neural-geometry-object.md` 的验收标准对齐）。
15. `expected_parent` 用文件锁或版本目录落地真正的并发保证。

---

## 8. 附：本次复现命令

```bash
python -m venv .venv && .venv/bin/pip install -e ".[test,dev]"

.venv/bin/python -m pytest -q                                   # 175 passed, 2 skipped
.venv/bin/python -m pytest -q --cov --cov-report=term-missing   # 94.47 %（阈值 90 生效）
.venv/bin/ruff check .                                          # clean
.venv/bin/mypy                                                  # clean (16 files)
.venv/bin/python examples/parametric_hull/run_example.py --output /tmp/hull-demo

# P1-1 / P1-2
.venv/bin/tensorcad-hull --nx 4001 --nz 3 --no-export --json    # ValueError, traceback
.venv/bin/tensorcad-hull --nx 2002 --nz 3 --output /tmp/partial # 留下 hull.obj + offsets.csv
# P1-3
.venv/bin/tensorcad-neural verify --object docs/experiments/2026-10-24-neural-geometry-object/bundle
# P2-4
.venv/bin/python -c "from tensorcad import *; m=generate_hull(HullParameters(nx=3,nz=3)); print(m.validation.warnings)"
# P1-4（marching tetrahedra 网格上的静水力）
.venv/bin/python -c "import numpy as np; from tensorcad.neural.extract import extract_surface; \
  from tensorcad import compute_hydrostatics; m=extract_surface(lambda p: np.linalg.norm(p,axis=1)-0.5, (-0.58,)*3, (0.58,)*3, 48); \
  compute_hydrostatics(m, 0.58, closed=True)"

# §4：真实权重包校验 + 抽取器独立复测（均不需要 torch）
.venv/bin/python -c "from tensorcad.neural_geometry import load_bundle; \
  m=load_bundle('docs/experiments/2026-10-24-neural-geometry-object/bundle'); print(m.revision, m.summary())"
.venv/bin/python -c "import numpy as np; from tensorcad.neural.extract import extract_surface, summarise; \
  m=extract_surface(lambda p: np.linalg.norm(p,axis=1)-0.5, (-0.58,)*3, (0.58,)*3, 96); r=summarise(m,96); print(r.to_dict())"
```

> 注：`.venv/`、`results/` 已在 `.gitignore` 中。本文未修改任何源码，只新增了这一份分析文档。
> 本次未覆盖：(1) torch 未安装，故 `fit` / `edit` 的**训练**未重跑——相关数字来自仓库内真实实验产物，我复核了产物自身的哈希与指标一致性，但**没有**重新训练来独立确认门限值是否仍成立；(2) `mypy` / `ruff` 的 3.10–3.13 矩阵未在本机全跑，只在 3.11 上验证；(3) 未跑 `tests/test_cli_neural.py` 与 `tests/test_neural_pipeline.py`（需 `neural` extra）。
