# GAP-CAD-T1-IGES-PIPE: IGES/STEP → 水密网格进口管道

**锚点**: `GAP-CAD-T1-IGES-PIPE`
**状态**: ✅ 已实现（v0.25.0）
**金标回归**: `DTMB5415_Hull_exact.stl`

## 概述

本管道实现从工业标准 CAD 格式（IGES 5.3 / STEP AP203）到水密三角形网格的完整导入链，
服务于 CFD/水动力学求解器的几何输入需求。管道完全用 Python 标准库实现（零第三方依赖），
可选的 OpenCASCADE 后端通过 `tensorcad[ocp]` 提供。

## 架构

```
IGES/STEP 文件
    │
    ├─ trimmed route (iges_trim.py + trim_mesh.py)
    │   ├─ 解析 IGES 5.3 B-spline 曲面 + 裁剪曲线
    │   ├─ Gauss-Newton 反演裁剪点到曲面参数空间
    │   ├─ 参数域三角剖分（带孔洞的耳切法）
    │   └─ 评估曲面生成三角形
    │
    ├─ patch-set fallback (exchange.py)
    │   └─ 未裁剪的 128 型曲面 → 均匀采样 → 三角化
    │
    ├─ STEP B-rep route (step_brep.py)
    │   └─ 解析 ISO-10303-21 → 解析曲面 → 参数域三角剖分
    │
    ▼
mesh_repair.py 修复链
    ├─ weld（顶点合并）
    ├─ drop_degenerate_faces（删除零面积面）
    ├─ drop_duplicate_faces（删除重复面）
    ├─ manifoldize（非流形顶点分裂）
    ├─ orient_consistently（一致定向）
    └─ fill_holes（平面孔洞填充）
    ▼
水密网格 + 审计报告
```

## 七道门（Gate Verdicts）

| 门 | 判据 | 通过条件 |
|---|---|---|
| G1 | 全部面读取 | `faces_tessellated == faces_in_file` |
| G2 | 裁剪在曲面上 | 最大反演残差 ≤ 容差（默认 1e-6） |
| G3 | 水密 | 自由边数 = 0 |
| G4 | 流形 + 一致定向 | 非流形边 + 定向冲突 = 0 |
| G5 | 正体积 | 有符号体积 > 0 |
| G6 | 单壳 | 组件数 = 1 且封闭组件数 = 1 |
| G7 | 无退化 | 退化面 + 重复面 = 0 |

## CLI 用法

```bash
# 安装后可用
tensorcad-import hull.igs --out-dir runs/import

# 或直接作为模块运行
python -m tensorcad.geomquality.import_cli hull.igs --out-dir runs/import

# 批量导入
tensorcad-import hull1.igs hull2.step hull3.iges --out-dir runs/batch

# 允许开放曲面（螺旋桨、升力面）
tensorcad-import propeller.igs --allow-open --out-dir runs/prop
```

## Python API

```python
from tensorcad.geometry import load_cad_geometry

# 导入 IGES/STEP 文件
mesh, report = load_cad_geometry("hull.igs", u_samples=96, v_samples=48)

# 检查门 verdicts
for gate_name, verdict in report["gates"].items():
    print(f"{gate_name}: {'PASS' if verdict['passed'] else 'FAIL'}")

# 修复后指标
after = report["after_repair"]
print(f"watertight: {after['watertight']}")
print(f"volume: {after['signed_volume']} m³")
print(f"triangles: {after['triangles']}")
```

## 金标回归：DTMB5415

`DTMB5415_Hull_exact.stl` 是平台的标准船体参考网格（164,926 三角形，83,350 顶点）。
管道通过以下步骤验证：

1. 从 STL 读取三角形
2. 用 NURBS 插值生成 IGES B-rep（`gold_loft.py`）
3. 通过管道导入 IGES
4. 验证：水密 + 单壳 + 边界框与金标一致（误差 < 1% 船长）

测试结果：`test_dtmb5415_loft_imports_watertight_and_tracks_the_gold_surface` ✅ PASS

## 基准测试

```bash
# 运行完整基准（生成 IGES 文件 + 导入 + 体积对比）
python benchmarks/verified/iges_import/run_benchmark.py --out-dir runs/gap_t1_iges

# 输出：
# - artifacts/gap_t1_iges/inputs/*.iges  (生成的 IGES 文件)
# - runs/gap_t1_iges/*.report.json       (每文件审计报告)
# - runs/gap_t1_iges/*.stl               (修复后网格)
# - runs/gap_t1_iges/manifest.json       (汇总)
# - runs/gap_t1_iges/evidence.json       (证据文件，用于回归测试)
```

## 测试覆盖

| 测试文件 | 状态 | 说明 |
|---|---|---|
| `test_iges_gold_regression.py` | ✅ 9/9 PASS | 金标回归 + 参考积分验证 |
| `test_iges_baseline.py` | ✅ 基础通过 | IGES 解析器单元测试 |
| `test_iges_trim_import.py` | ✅ 基础通过 | 裁剪导入 + 修复链测试 |
| `test_iges_units_and_patch_quarantine.py` | ✅ 基础通过 | 单位处理 + 歧义面隔离 |
| `test_exchange_verification.py` | ✅ 通过 | STEP 面片交换验证 |
| `test_mesh_repair.py` | ✅ 通过 | 网格修复链测试 |

## 依赖锁定

核心包（`tensorcad`）零依赖声明成立：
- `import tensorcad` 在无 numpy/torch 环境下可导入
- IGES 解析器（`iges.py`）纯标准库
- STEP B-rep 解析器（`step_brep.py`）纯标准库
- 裁剪导入（`trim_mesh.py`）纯标准库
- 网格修复（`mesh_repair.py`）纯标准库

测试依赖（`numpy`）仅在 `test` extra 中声明。

## 已知局限

1. **IGES 支持范围有限**：仅支持 128 型有理 B-spline 曲面 + 143/144 型裁剪面
2. **STEP 仅支持解析 B-rep**：不支持 FACETED_BREP（需通过 exchange.py 的 patch-set 路线）
3. **孔洞填充为平面**：曲面孔洞用平面片填充，偏差在报告中记录
4. **不自交修复**：自交面拒绝而非修复（需 fTetWild 级重网格）
5. **无 OpenCASCADE 时不支持复杂 B-rep**：ADVANCED_FACE/CYLINDRICAL_SURFACE 等需 `tensorcad[ocp]`

## 变更日志

- **v0.25.0** (2026-11-13): 新增 `tensorcad-import` CLI 入口；完善 IGES/STEP 管道文档
- **v0.24.0**: 金标回归测试通过（DTMB5415）
- **v0.23.0**: 裁剪导入管道完整实现
- **v0.22.0**: STEP B-rep 解析器实现
- **v0.21.0**: 网格修复链实现（GAP-CAD-009）
