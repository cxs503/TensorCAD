# TensorCAD 项目分析报告

- 分析日期：2026-10-24
- 分析对象：`cxs503/TensorCAD`，`main` 分支 `febcb9c`（0.2.0）
- 分析方式：静态阅读 + 本地全量构建/测试/基准复现（未修改源码）
- 结论一句话：**这是一个工程质量明显高于同阶段研究原型的几何基础库；核心承诺（stdlib-only、确定性、闭合实体、可测量验证）都有代码/测试支撑，但存在 1 个可复现的功能缺陷、若干门禁与文档不一致，以及"只有一种简化船型、无 B-rep"的能力边界。**

---

## 1. 项目定位

TensorCAD 不是一个桌面 CAD 的复刻，而是 Tensor 科研计算生态里的**几何前置层**：

| 层 | 角色 |
|---|---|
| TensorCAD | 参数化几何、验证、网格与预处理接口 |
| TensorShipDesign（计划中） | 基于 TensorCAD 的船舶设计应用 |
| TensorLBM / FEM / DEM / BEM | 消费几何与网格的求解器 |
| TensorPlatform | 编排 agent 任务、仿真、评估与迭代 |

设计原则在代码中基本都落到了实处：headless-first、可复现、先几何后优化、复用成熟内核（OCP 作为可选适配器）、边界清晰、**示例必须跑真实代码并测量**（最后由脚本自己的断言兜底）。

---

## 2. 代码结构与分层

```
tensorcad/
  __init__.py        65   公共 API 门面（__all__ 共 22 个符号）
  geometry.py       183   Mesh / BoundingBox / MeshBuilder（精确坐标去重）
  parameters.py     122   ParameterSpec / RatioSpec，单位与范围，一次报全部错误
  mesh.py           297   通用网格算法：边统计、水密、定向、焊接、退化面、积分
  validation.py     292   Check / ValidationReport / inspect_geometry / 尺寸比对
  hydrostatics.py   215   面元积分：排水体积、湿表面积、水线面、LCB/LCF/VCB、Cb/Cp
  export.py         167   OBJ/CSV/JSON 写出，OBJ 读回，17 位有效数字往返
  hull.py           436   唯一的船型专用模型：参数→截面环→放样→封盖→验证→静水力
  cli.py             97   tensorcad-hull 入口
  neural_geometry.py197   NGOM 清单安全加载（只读元数据，不反序列化可执行对象）
  adapters/          16   OCP/gmsh 适配器占位（核心禁止 import）
examples/parametric_hull/  可执行基准，自带断言
examples/neural_geometry/  torch 条件神经场研究与"球形凹陷"编辑实验
docs/research/             连续几何场 / NGOM 元数据 / 神经几何对象 三份设计文档
schemas/                   NGOM manifest JSON Schema v0.1
tests/                     9 个文件，1,247 行
```

源码 2,087 行 / 测试 1,247 行（测试/源码比 ≈ 0.60），模块总数 11，函数与类 115 个。**依赖方向是单向的**：`hull → mesh/geometry/validation/hydrostatics/parameters`，`cli → hull/export`，核心不反向依赖 CLI，也不依赖任何第三方库（有一条测试专门 reload 模块后断言 `numpy/scipy/OCP/trimesh` 等未被加载）。

### 生成流水线

```
HullParameters.validate()
      │
      ├─ build_mesh: 每站生成一条闭合截面环（右侧水线→龙骨→左侧水线）
      │              相邻站放样成四边形，水线处加甲板封盖
      │              MeshBuilder 按精确坐标共享顶点（不做容差焊接）
      │
      ├─ meshlib.orient_outward: 广度优先统一绕向 + 按连通分量体积翻转
      │
      ├─ 验证三件套：inspect_geometry（拓扑/朝向/退化）
      │              check_dimensions（L/B/T/水线面与请求比对）
      │              check_consistency（偏移表、比例、nx 奇偶）
      │              → 任何 error 直接 raise，绝不静默返回坏网格
      │
      ├─ compute_hydrostatics: 散度定理面元积分（闭曲面）
      │
      └─ HullModel → export_hull: hull.obj / offsets.csv / metadata.json
```

关键设计选择是**用构造保证拓扑，而不是事后修补**：一条截面环共享龙骨点，避免了镜像焊接导致的绕向反转与容差伪影；这正是 0.2.0 修掉的两个真实缺陷的根因。`z = -draft * (1 - fraction)` 的写法也是为了让水线严格落在 `z = 0`，而不是靠浮点抵消。

---

## 3. 实测结果（本次在本机复现）

### 3.1 门禁与测试

| 项目 | 结果 |
|---|---|
| `pytest` | **142 passed**，1.49 s |
| 分支覆盖率 | **96%**（geometry/hull/hydrostatics/mesh/parameters/validation 均 98–100%；`neural_geometry` 79%） |
| `ruff check .` | All checks passed |
| `mypy`（strict，packages=tensorcad） | Success: no issues found in 11 source files |
| 示例基准 | 全部检查 PASS，退出码 0 |
| 核心运行时依赖 | **0**（stdlib），`test`/`dev`/`neural` 为可选 extra |

示例脚本 `examples/parametric_hull/run_example.py` 不只是"跑通"：它比较解析体积、检查 LCB = LCF、做五级网格收敛研究，并在任一断言失败时返回非零。这是仓库里最有价值的工程习惯。

### 3.2 数值基线（默认 L=10, B=2, T=0.8）

```
网格 1001 顶点 / 1998 面
体积       7.446554 m^3   （解析 7.512698，-0.88%）
湿表面积  22.838010 m^2
水线面积  14.061669 m^2
LCB/LCF/VCB  5.0000 / 5.0000 / -0.3170 m
Cb / Cp      0.4654 / 0.7031
nx=161,nz=41 → 体积 7.502488（-0.136%），误差随加密单调下降
```

### 3.3 性能（纯 Python）

| 网格 | 顶点/面 | 生成+验证耗时 |
|---|---|---|
| 41 × 13 | 1,001 / 1,998 | 28 ms |
| 161 × 41 | 12,961 / 25,918 | 434 ms |
| 401 × 101 | 80,401 / 160,798 | 3.7 s |

复杂度是 O(F)，但常数因子偏大（验证阶段每个面被多次遍历）。默认导出（含一次加密重算）约 126 ms，对 agent 调用足够快；`401×101` 已经接近人可接受的交互上限，而参数上限 `4001×513` 在纯 Python 下实际上不可用。

---

## 4. 发现的问题

### P1 — `estimate_discretisation_error` 在加密后越界时抛异常，而不是返回 `skipped`（已复现）

`tensorcad/hull.py:405-412`：

```python
if p.nx * p.nz > MAX_MESH_SIZE_FOR_CONVERGENCE:   # 400_000
    return {"status": "skipped", ...}
refined = HullParameters(p.length, p.beam, p.draft, 2 * p.nx - 1, 2 * p.nz - 1)
...
fine_mesh, _ = build_mesh(refined)                # ← 这里再次 validate()
```

`build_mesh` 会校验 `nx ≤ 4001`、`nz ≤ 513`，但加密后的参数没有被 `MAX_MESH_SIZE_FOR_CONVERGENCE` 这个"体积守卫"覆盖：

- `nx ≥ 2002` → `2*nx-1 = 4003 > 4001`
- `nz ≥ 258` → `2*nz-1 = 515 > 513`
- 只要 `nx*nz ≤ 400_000`（例如 `2002×3 = 6006`），守卫就不会拦住。

复现：

```bash
$ tensorcad-hull --nx 4001 --nz 3 --no-export --json
...
  File "tensorcad/hull.py", line 412, in estimate_discretisation_error
    fine_mesh, _ = build_mesh(refined)
ValueError: nx must be <= 4001 stations, got 8001
exit=1
```

影响面比看上去大：CLI 的**默认路径**（`--no-convergence` 未给出）和 `--json` 都会调用 `metadata(with_discretisation=True)`；而这个 `ValueError` 发生在 `generate_hull` 的 `try` 块之外，所以用户看到的是一段完整 traceback，而不是"参数非法"的友好错误。更反直觉的是：**模型本身完全合法且已生成成功**，只是收敛估计失败。

修复建议（约 6 行）：

```python
refined_nx, refined_nz = 2 * p.nx - 1, 2 * p.nz - 1
if (p.nx * p.nz > MAX_MESH_SIZE_FOR_CONVERGENCE
        or refined_nx > MAX_LONGITUDINAL_STATIONS
        or refined_nz > MAX_VERTICAL_LEVELS):
    return {
        "status": "skipped",
        "reason": f"refined grid {refined_nx}x{refined_nz} would exceed the "
                  f"{MAX_LONGITUDINAL_STATIONS}x{MAX_VERTICAL_LEVELS} parameter limits",
    }
```

并补一条参数化回归测试（`nx=4001, nz=3`、`nx=3, nz=513`、`nx=2002, nz=3`）。

附带一个潜在陷阱：`HullParameters(p.length, p.beam, p.draft, 2*nx-1, 2*nz-1)` 丢掉了 `p.water_density`，重新落到默认 1025。当前只比较体积/面积，所以结果不受影响；但这是典型的"今天无害、明天有害"的隐式默认值，应显式传入。

### P2 — 覆盖率"门禁"并不存在

- `.github/workflows/tests.yml:25` 的步骤名是 *Tests with coverage gate*，实际命令是 `pytest --cov --cov-report=term-missing`；
- `pyproject.toml` 的 `[tool.coverage.report]` 没有 `fail_under`，CI 也没有 `--cov-fail-under`；
- `CHANGELOG.md:58` 声称 "coverage (90 % gate, currently 99 %)"，`README.md` 声称 129 个测试。

实测：总覆盖率 96%，`neural_geometry` 79%，测试数 142。也就是说**门禁名称、文档数字、真实状态三者不一致**。今天如果直接加 `fail_under = 90`，CI 会因为 `neural_geometry` 的 79% 而失败——正确顺序是先把该模块的错误分支补齐，再开门禁。

### P2 — 没有任何"网格太粗"的告警

```python
m = generate_hull(HullParameters(nx=3, nz=3))
m.validation.warnings          # []
m.hydrostatics.displacement_volume  # 4.8284，解析值 7.5127（-35.7%）
```

`validation.warnings` 是文档给出的"下一步该看什么"的入口，但它只覆盖比例失真、偶数列 `nx` 的梁宽亏欠和顶点流形。粗网格带来的**数值失真没有任何提示**，而生成器其实可以零成本地知道这一点：`TAPER_INTEGRAL` 就在 `hull.py` 里，解析体积可用 O(1) 算出来，和网格积分一比即可。建议在 `generate_hull` 中增加一条 warning（例如相对误差 > 5% 时提示提高 `nx/nz`，并引用 `estimate_discretisation_error`），这比让用户自己发现 `metadata.json` 里的收敛字段友好得多。

### P2 — 缺少 `py.typed`

`mypy strict` 只对仓库内部生效；`pyproject.toml` 里 `packages.find` 正常，但没有 `tensorcad/py.typed`。作为"agent-callable"的基础库，下游按 PEP 561 安装后拿不到类型信息，strict 模式的收益无法传递给使用者。一行文件即可解决。

### P3 — 静水力 API 的水线语义可以被误用

`compute_hydrostatics(mesh, waterline_z=...)` 的 `waterline_z` 只用于**分类哪些面是湿面**；排水体积始终是**整个闭合网格**的散度定理积分。docstring 已声明"网格是由平面水线封盖的湿表面、没有面跨越水线"，但 API 并不阻止调用者传入任意 `waterline_z`——比如对一个有干舷的船体传真实水线，会得到"湿面积正确、体积错误"的混合结果，且没有任何错误。建议要么断言 `waterline_z == max(z)`，要么显式实现水线切割。

### P3 — OBJ 往返丢失单位与坐标系

`write_obj` 把 `units`/`coordinate_frame` 写成注释，`read_obj` 返回 `Mesh(..., units="m", coordinate_frame="unspecified right-handed frame")`（`export.py:101`），并没有解析这些注释。顶点表确实逐位往返（有测试锁定），但**语义元数据没有往返**，导入一个英尺制或局部坐标系的模型会被静默标成米制。建议 `read_obj` 解析 `# units:`/`# coordinate frame:`，或至少提供显式参数。

### P3 — 重复实现与重复计算

- `mesh._tet_volume`（`mesh.py:88`）与 `hydrostatics._tetra_volume`（`hydrostatics.py:177`）是同一公式的两份实现；
- `compute_hydrostatics` 对每个面依次调用 `face_centroid`、`triangle_area`、`face_normal`（`hydrostatics.py:104-116`），每个函数都重新做一遍叉积，湿面还要再来一次。一个返回 `(area, normal, centroid)` 的辅助函数可以把大网格上的静水力耗时明显压下来，同时消除重复代码。

### P3 — CLI 只有一半的异常处理

`cli.py` 的 `try` 只包住 `generate_hull`（返回码 2），`export_hull` / `metadata(with_discretisation=True)` 的异常直接冒泡成 traceback（P1 就是现成例子）。`--json` 模式下错误也没有结构化输出。建议把整段运行包起来，stdout 输出 `{"error": ..., "stage": ...}`，stderr 保留人类可读信息，返回非零码。

### P3 — 神经几何研究线与主线是"两个仓库"

`docs/research/` 三份文档 + `examples/neural_geometry/` + 独立的 30 分钟 CI workflow 构成了一条研究支线，但：

- README 完全没有提及它，新读者会以为仓库只有船体；
- `tensorcad.neural_geometry` 是纯 stdlib 的清单加载器，却不在 `__all__` 里（可能是有意的"研究态 API"，但没有说明）；
- `schemas/neural-geometry-manifest-v0.1.schema.json` 的 `format_version` 允许 `0.01.0`（正则 `^0\.[0-9]+\.[0-9]+$`），而 `parse_manifest` 禁止前导零并要求语义化三段式；schema 的 `additionalProperties: true` 与加载器的严格字段检查也不完全对应。建议把加载器的规则视为唯一事实来源，反向修 schema。

安全方面这条线做得不错：拒绝 `NaN`/`Infinity` 常量、SHA-256 校验权重、限制路径逃逸、拒绝反序列化可执行对象。尚未覆盖的是**文件大小上限**（`payload.read_bytes()` 一次性读入内存）和重复 JSON key 检测，属于低优先级的健壮性项。

### 说明：一处"看起来像问题但其实是设计"的地方

`check_dimensions` 里的 `beam_matches_request` 带一段"采样亏欠允许带"逻辑，但在 `_dimension_target` 中目标梁宽已经被乘上"网格可达的最大 taper"，所以实测亏欠通常是 0，这段允许带基本不生效；真正的偶数列梁宽提示在 `check_consistency` 的 warning 里。功能没错，但两块逻辑存在重叠，长期看应该合并或加注释说明分工。

---

## 5. 架构与工程评价

**做得好的：**

1. **承诺即测试**：stdlib-only 有测试、闭合实体有测试、OBJ 17 位往返有测试、同参数同字节有 CI 步骤。这比大多数 0.2.0 的研究库扎实。
2. **失败的方式正确**：参数错误一次报全；几何不可用直接 raise；warning 与 error 分级；示例脚本自己会失败。
3. **分层可执行**：`adapters` 是空包 + 文档边界，核心不 import 它，`pyproject` 用 extras 表达可选后端。
4. **数值诚实**：明确区分"网格积分值"与"解析值"，并把 0.88%/0.14% 的离散误差写进元数据，而不是宣称精确。
5. **CI 维度合理**：3.10–3.13 矩阵、ruff + mypy strict + 覆盖率 + CLI 产物 + 可复现性 + 基准，覆盖面完整。

**真正的短板：**

1. **能力面窄**：只有一个前后对称、无方尾、无干舷、龙骨带 knuckle 的解析船型。它证明的是"几何基础设施可用"，而不是"能设计船"。
2. **无 B-rep / STEP**：issue #4 仍未落地，`adapters/ocp` 只是文档里的接缝。对工程用户来说，STEP 才是有交换价值的东西。
3. **无面网格质量指标**：长宽比、最小角、雅可比、法向跳变都没有，而这些才是"求解器能不能用"的第二道门（第一道水密已经解决）。
4. **纯 Python 的可扩展性上限**：`4001×513` 的参数上限在实现上不可达；缺少可选的 numpy 加速路径或分块流式校验。
5. **研究线与工程线未合流**：神经场实验产出的表面还没有走 `inspect_geometry`/`compute_hydrostatics` 这条已验证的验收管线，这正是两条线最自然的汇合点。

---

## 6. 建议路线（按投入产出排序）

**立刻（0.2.1，半天量级）**

1. 修 P1：加密网格越界时返回 `skipped`，并显式传递 `water_density`；补回归测试。
2. 加 `tensorcad/py.typed`。
3. CLI 全流程异常处理 + 结构化错误输出。
4. 同步 README/CHANGELOG 的测试数与覆盖率数字；要么删掉"90% gate"的表述，要么真正加上 `fail_under`（先补 `neural_geometry` 的错误分支测试）。

**短期（0.3）**

5. 粗网格数值告警（用解析体积做 O(1) 对照）。
6. 面网格质量检查（长宽比、最小角、退化三角形阈值、法向一致性），并入 `ValidationReport`。
7. OBJ 单位/坐标系往返；`read_obj` 对非有限坐标直接拒绝。
8. 为 `metadata.json` 补一份 JSON Schema（目前只有神经几何有 schema，船体元数据 schema_version 2.0 只是字符串）。
9. 合并重复的四面体体积/面几何计算，做一次性能剖析。

**中期（0.4，对应 issue #4 / #6）**

10. `tensorcad.adapters.ocp`：B-rep 构建、内核有效性检查、STEP 往返验证；核心仍保持 stdlib-only。
11. 船型 v2：方尾、干舷、前后非对称、真实剖面放样，让静水力从"数值演示"变成"可用于设计比较"。
12. 神经场 → 网格提取 → 复用 `inspect_geometry` / `compute_hydrostatics` 的验收管线；把"场误差"和"提取网格误差"分开报告（这与 `docs/research/neural-geometry-object.md` 的验收标准一致）。
13. 用一次性体积网格交接（gmsh 适配器）打通"几何 → 求解器"的最后一段。

---

## 7. 附：本次复现命令

```bash
python -m venv .venv && .venv/bin/pip install -e ".[test,dev]"
.venv/bin/python -m pytest -q                       # 142 passed
.venv/bin/python -m pytest -q --cov --cov-report=term-missing   # 96%
.venv/bin/ruff check .                              # clean
.venv/bin/mypy                                      # clean
.venv/bin/python examples/parametric_hull/run_example.py --output /tmp/hull-demo
.venv/bin/tensorcad-hull --nx 4001 --nz 3 --no-export --json   # 复现 P1
```

> 注：`.venv/` 与 `results/` 已在 `.gitignore` 中；本次分析未修改任何源码。
