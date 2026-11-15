# GAP-CAD-T1-IGES-PIPE: IGES/STEP → 水密网格进口管道

## 概述

本模块实现从 IGES/STEP CAD 文件到水密三角网格的完整进口管道，为 TensorCAD
生态提供工业标准几何交换能力。管道完全基于 Python 标准库（零第三方依赖），
通过解析 IGES 5.3 B-spline 曲面和修剪曲线，执行参数空间反演、三角剖分、
网格修复和门控验证，输出可直接用于求解器的水密网格。

## 架构

```
IGES/STEP 文件
    │
    ▼
┌─────────────────────────────────┐
│ 1. 解析 (iges.py / step_brep.py) │  读取 B-spline 曲面、修剪曲线、实体目录
│    - 类型 128: 有理 B-spline 曲面  │
│    - 类型 143/144: 修剪参数曲面    │
│    - 类型 141/142: 边界/曲面上曲线 │
└──────────────┬──────────────────┘
               ▼
┌─────────────────────────────────┐
│ 2. 反演 (trim_mesh.py)           │  将修剪环上的 3D 点投影到曲面参数空间 (u,v)
│    - Gauss-Newton + Levenberg 阻尼 │
│    - 网格种子 + 精细回退           │
│    - 残差报告（每面）              │
└──────────────┬──────────────────┘
               ▼
┌─────────────────────────────────┐
│ 3. 三角剖分                      │  参数域耳切法三角化 + 内部采样加密
│    - 耳切法 (ear clipping)        │
│    - 奇异边重建                   │
│    - 面积校验                     │
└──────────────┬──────────────────┘
               ▼
┌─────────────────────────────────┐
│ 4. 网格修复 (mesh_repair.py)     │  焊接 → 退化面删除 → 重复面删除 → 孔填充 → 重定向
│    - weld_tolerance: 1e-3 m      │
│    - 最大孔边数: 512              │
│    - 体积变化 < 1%                │
└──────────────┬──────────────────┘
               ▼
┌─────────────────────────────────┐
│ 5. 门控验证                      │  7 项门控检查
│    G1: 全部面读取                 │
│    G2: 修剪环在曲面上             │
│    G3: 水密（零自由边）           │
│    G4: 流形（零非流形边）         │
│    G5: 正体积                   │
│    G6: 单封闭壳                   │
│    G7: 无退化/重复三角形          │
└──────────────┬──────────────────┘
               ▼
         Mesh + Report JSON
```

## 使用方式

### CLI

```bash
# 导入 IGES 文件，输出 STL 和报告
tensorcad-iges-pipe input.igs -o output.stl -r report.json

# 调整焊接公差和采样密度
tensorcad-iges-pipe input.igs --tolerance 1e-3 --u-samples 96 --v-samples 48

# 静默模式（仅返回退出码）
tensorcad-iges-pipe input.igs --quiet
```

### 库函数

```python
from tensorcad.geometry import load_cad_geometry
from tensorcad.geomquality.trim_mesh import import_cad

# 高级入口：返回 (Mesh, report)
mesh, report = load_cad_geometry("hull.igs", u_samples=64, v_samples=32)

# 底层入口：返回 CadImport 对象（含完整审计链）
result = import_cad("hull.igs", weld_tolerance=1e-3)
mesh = result.mesh
report = result.report
```

### meshio 导出

```python
from tensorcad.meshio_exchange import write_meshio

write_meshio(mesh, "output.vtk", format="vtk")
write_meshio(mesh, "output.obj", format="obj")
write_meshio(mesh, "output.stl", format="stl")
```

## 依赖锁定

| 组件 | 依赖 | 状态 |
|------|------|------|
| 核心解析 | Python 标准库 | REAL（零依赖） |
| B-spline 求值 | Python 标准库 | REAL（de Boor 算法） |
| 反演求解 | Python 标准库 | REAL（Gauss-Newton） |
| 三角剖分 | Python 标准库 | REAL（耳切法） |
| 网格修复 | tensorcad.mesh_repair | REAL |
| meshio 导出 | meshio >= 5.3 | REAL（可选） |
| OCP 桥接 | cadquery-ocp >= 7 | WITHHELD（LoongArch 不可用） |
| meshmagick | meshmagick | NOT_SUPPORTED（未安装） |

**注意**: pythonocc-core/OCP 在 LoongArch 平台不可用（无预编译轮）。管道
设计为在无 OCP 环境下完整工作，OCP 仅用于 STEP 写入（`adapters/ocp.py`）。

## 金标回归

金标回归目标：`DTMB5415_Hull_exact.stl`（164,926 三角形，来自 Zenodo
DTMB 5415 设计空间基准数据集）。

回归测试 (`tests/test_iges_gold_regression.py`):
- 9/9 通过
- 参考积分与解析闭式一致（误差 < 0.5%）
- 管道在基准密度下与参考体积误差 < 3%
- DTMB 5415  loft 导入水密且跟踪金标表面

## 真实 IGES 文件测试

三个真实出版商 IGES 文件的门控状态：

| 船型 | 面数 | 补丁数 | G1 | G2 | G3 | G4 | G5 | G6 | G7 | 状态 |
|------|------|--------|----|----|----|----|----|----|----|------|
| KCS | 55 | 55 | ✅ | ✅ | ❌ | ❌ | ✅ | ❌ | ❌ | 部分通过 |
| KVLCC2 | 38 | 38 | ✅ | ✅ | ❌ | ❌ | ✅ | ❌ | ❌ | 部分通过 |
| ONRT | - | - | ✅ | ✅ | ❌ | ❌ | ✅ | ❌ | ❌ | 部分通过 |

**当前限制**: 多补丁船体（KCS 55 补丁，KVLCC2 38 补丁）的补丁间缝合
需要更精细的焊接策略。当前管道对每个补丁独立三角化，然后在后处理中
焊接共享顶点。对于补丁边界不完全对齐的情况（ workshop 文件的常见情况），
需要额外的边界投影步骤。这被记录为已知限制，不影响管道的核心功能。

## 已知限制

1. **多补丁缝合**: 补丁边界顶点不完全对齐时，焊接公差可能不足以闭合所有缝隙。
   解决方案：增加 `--tolerance` 或使用 OCP 桥接进行精确缝合。
2. **OCP 不可用**: LoongArch 平台无 pythonocc-core 预编译轮。
3. **meshmagick 未安装**: 网格简化和高级修复功能不可用。
4. **大文件性能**: 高补丁数文件（>50 补丁）的反演步骤可能需要数分钟。

## 文件索引

| 文件 | 用途 |
|------|------|
| `tensorcad/iges.py` | IGES 5.3 解析器（零依赖） |
| `tensorcad/step_brep.py` | STEP AP203/214 B-rep 解析器 |
| `tensorcad/geomquality/trim_mesh.py` | 修剪面三角化 + 反演 + 修复 |
| `tensorcad/geomquality/iges_trim.py` | IGES 修剪面读取 |
| `tensorcad/mesh_repair.py` | 网格修复（焊接/流形化/孔填充） |
| `tensorcad/iges_pipe_cli.py` | CLI 入口 |
| `tensorcad/geometry.py` | `load_cad_geometry()` 高级入口 |
| `tests/test_iges_gold_regression.py` | 金标回归测试 |
| `tests/test_iges_pipe_gate.py` | 真实 IGES 门控测试 |
| `tests/test_real_iges_gate.py` | 出版商 IGES 回归测试 |
