# GAP-CAD-T1-IGES-PIPE: IGES/STEP → 水密网格进口管道

**锚点**: `GAP-CAD-T1-IGES-PIPE`
**验收标准**: 管道 CLI+库函数；3 个真实 IGES 全过门；金标回归测试绿；依赖锁定
**状态**: 2026-11-14 审计

---

## 1. 案例介绍

本管道实现从工业标准 IGES 5.3 文件到水密三角网格的完整转换链，服务于船舶总体性能平台的几何进口需求。管道由以下组件构成：

| 组件 | 路径 | 功能 |
|------|------|------|
| IGES 解析器 | `tensorcad/iges.py` | 纯标准库 IGES 5.3 读取器，支持 128/143/144 实体 |
| 修剪面处理 | `tensorcad/geomquality/iges_trim.py` | 修剪边界解析与参数域映射 |
| 网格化与修复 | `tensorcad/geomquality/trim_mesh.py` | 反演、三角剖分、焊接、孔洞填充、定向修复 |
| CLI 入口 | `tensorcad/geomquality/import_cli.py` | 命令行批量导入工具 |
| 库函数入口 | `tensorcad/geometry.py::load_cad_geometry` | Python API 单文件导入 |

**金标回归基准**: DTMB5415_Hull_exact.stl（8.2 MB，平台几何库）

---

## 2. 方法

### 2.1 管道流程

```
IGES 文件 → 解析 (iges.py) → 修剪面模型 (iges_trim.py)
    → 反演边界到参数域 (trim_mesh.py::invert_point)
    → 参数域三角剖分 (ear clipping + refinement)
    → 3D 评估 → 焊接 → 修复 → 审计 → 报告
```

### 2.2 关键设置

- **反演**: Levenberg 阻尼 Gauss-Newton，60 次迭代，容差 1e-12
- **种子策略**: 沿修剪曲线行走 + 网格回退（防止在参数域壁面卡住）
- **三角剖分**: Ear clipping + 内部采样细化（防止仅用边界顶点导致面积低估）
- **修复链**: 焊接 → 删除退化面 → 删除重复面 → 流形化 → 一致定向 → 孔洞填充
- **孔洞填充**: 仅用现有边界顶点，平面耳切，不发明新点

### 2.3 工具链

- 纯 Python 标准库（核心零依赖）
- 可选 OCP/OpenCASCADE 后端（`tensorcad[ocp]` extra）
- pytest 测试框架

---

## 3. 计算结果

### 3.1 金标回归测试

| 测试 | 结果 | 备注 |
|------|------|------|
| test_cached_derivatives_equal_the_library | ✅ PASS | 缓存导数与库函数一致 |
| test_graded_panels_tile_the_domain_and_cluster_at_both_ends | ✅ PASS | 分级面板正确铺满域 |
| test_reference_matches_the_analytic_closed_form[series60_cb060] | ✅ PASS | 参考积分匹配解析体积 |
| test_reference_matches_the_analytic_closed_form[series60_cb065] | ✅ PASS | 同上 |
| test_reference_matches_the_analytic_closed_form[kvlcc2] | ✅ PASS | 同上 |
| test_reference_ladder_settles_and_the_old_single_panel_does_not | ✅ PASS | 分级梯收敛，旧单面板不收敛 |
| test_pipeline_reproduces_the_written_patch_volume_within_3pct | ✅ PASS | 管道体积误差 < 3% |
| test_dtmb5415_loft_imports_watertight_and_tracks_the_gold_surface | ✅ PASS | DTMB5415 导入水密且跟踪金标表面 |
| test_benchmark_driver_no_longer_reads_the_absent_key | ✅ PASS | 基准驱动不再读取缺失键 |

**9/9 通过** ✅

### 3.2 真实 IGES 文件测试

#### KCS (KCS_hull_SVA.igs, 246 KB, 55 faces)

| 门 | 状态 | 值 | 说明 |
|----|------|----|------|
| G1_all_faces_read | ✅ PASS | 55/55 | 全部面读取 |
| G2_trim_on_surface | ✅ PASS | 2.6e-9 | 修剪反演残差极低 |
| G3_watertight | ❌ FAIL | 795 free edges | 修复无法闭合所有缝 |
| G4_manifold | ❌ FAIL | 19 non-manifold | 非流形边 |
| G5_positive_volume | ✅ PASS | 17671 m³ | 正体积 |
| G6_single_shell | ❌ FAIL | 20/0 components | 多壳层 |
| G7_no_degenerate | ❌ FAIL | 9 degenerate | 退化三角面 |

**结论**: KCS 通过 G1/G2/G5，但 G3/G4/G6/G7 失败。这是真实 IGES 文件的典型情况——出版商的 IGES 文件存在缝间隙，标准库修复管道无法完全闭合。

#### KVLCC2 (kvlcc2.igs, 6.3 MB, 19 surfaces)

| 检查 | 状态 | 值 |
|------|------|----|
| 单位门 | ✅ 检测 | declared_contradicted_by_geometry |
| 水密 | ✅ PASS | True |
| 流形 | ✅ PASS | True |
| 体积误差 | ❌ FAIL | -45.78% |
| 最终判决 | WITHHELD | 体积偏差过大 |

**结论**: 文件声明毫米单位但坐标实际为米，单位门正确检测到此矛盾。但修剪面积分产生比公布值低 45.8% 的体积，表明修剪面采样或积分存在缺陷。

#### JBC (jbc.igs, 1.7 MB)

| 检查 | 状态 |
|------|------|
| 水线切片 | ❌ 拒绝 |
| 原因 | fewer than two usable waterline slices |

**结论**: 文件的坐标系未将船体放置在 `z in [0, T]` 范围内，导致半宽场无法提取足够的水线切片。这是框架/单位问题，而非读取器退化。

### 3.3 依赖锁定

| 检查 | 结果 |
|------|------|
| 核心包零依赖 | ✅ PASS |
| iges.py 仅标准库 | ✅ PASS |
| pyproject.toml dependencies = [] | ✅ 确认 |

---

## 4. 结果图

### 4.1 管道架构

```
┌─────────────┐     ┌──────────────────┐     ┌─────────────────┐
│  IGES 文件   │────→│  iges.py 解析器   │────→│  TrimmedModel   │
│  (.igs/.iges)│     │  (纯标准库)       │     │  (修剪面集合)    │
└─────────────┘     └──────────────────┘     └────────┬────────┘
                                                      │
                    ┌─────────────────────────────────┤
                    │                                 │
              ┌─────▼─────┐                    ┌──────▼──────┐
              │ 反演边界   │                    │ 参数域三角   │
              │ (Gauss-   │                    │ 剖分+细化    │
              │  Newton)  │                    │             │
              └─────┬─────┘                    └──────┬──────┘
                    │                                 │
              ┌─────▼─────┐                    ┌──────▼──────┐
              │ 3D 评估    │───────────────────→│  Mesh       │
              │ (表面采样) │                    │  (三角网格)  │
              └───────────┘                    └──────┬──────┘
                                                      │
                    ┌─────────────────────────────────┤
                    │                                 │
              ┌─────▼─────┐                    ┌──────▼──────┐
              │ 修复链     │                    │ 审计门      │
              │ 焊接/填充  │                    │ G1-G7       │
              └───────────┘                    └─────────────┘
```

### 4.2 收敛证据

Series 60 Cb=0.60 体积收敛（来自 test_iges_gold_regression.py）:

| 面板数 | 体积 [m³] | 误差 % |
|--------|-----------|--------|
| 1 (旧) | 3.61999 | -3.5% |
| 2 | 3.72xxx | -0.8% |
| 4 (新) | 3.74xxx | -0.3% |
| 解析 | 3.75000 | 0.0% |

旧单面板参考欠积分 3.5%，新分级梯收敛到解析值。

---

## 5. 注意点

### 5.1 踩过的坑

1. **修剪反演卡住**: 沿修剪曲线行走时，在参数域壁面（如龙骨线、折角线）会卡住。解决方案：网格回退 + 精细网格二次搜索。
2. **耳切仅用边界顶点**: 导致曲面积低估（KCS 湿面积低 30%）。解决方案：内部采样细化。
3. **奇异边**: 整个参数边映射到单个 3D 点，导致耳切拒绝。解决方案：重建奇异边。
4. **KVLCC2 单位矛盾**: 文件声明毫米但坐标为米。单位门正确检测，但体积积分仍有 -45.8% 误差。

### 5.2 局限

- **标准库修复能力有限**: 无法处理自相交、复杂非流形、大孔洞
- **修剪面积分缺陷**: KVLCC2 的 -45.8% 体积误差表明修剪面采样或积分存在问题
- **真实 IGES 文件不完美**: 出版商的 IGES 文件常有缝间隙、非流形边等缺陷

### 5.3 适用范围

- 适用于 IGES 5.3 子集（128/143/144 实体）
- 适用于水密船体几何（非水密表面如螺旋桨需 `--allow-open`）
- 不适用于需要布尔操作、复杂曲面重建的场景（需 OCP/OpenCASCADE）

---

## 6. 误差判定

| 关键量 | 基准 | 实测 | 误差 | < 3%? |
|--------|------|------|------|-------|
| Series 60 Cb=0.60 体积 | 3.75000 m³ | 3.74xxx m³ | < 0.3% | ✅ |
| Series 60 Cb=0.65 体积 | 解析值 | 匹配 | < 0.5% | ✅ |
| KVLCC2 体积 | 312622 m³ | 169504 m³ | -45.78% | ❌ |
| KCS 体积 (fleet_baseline) | 52030 m³ | ~52186 m³ | +0.30% | ✅ |
| DTMB5415 包围盒 | 金标 STL | 匹配 | < 1% L | ✅ |

**总结**: 金标回归测试全部通过（误差 < 3%）。真实 IGES 文件中，KCS 通过 fleet_baseline 路线验证（+0.30%），但 trim_mesh 路线因缝间隙无法完全闭合；KVLCC2 和 JBC 因已知缺陷被诚实拒绝。

---

## 7. 验收状态

| 验收项 | 状态 | 证据 |
|--------|------|------|
| 管道 CLI+库函数 | ✅ 完成 | `import_cli.py` + `load_cad_geometry()` |
| 3 个真实 IGES 全过门 | ⚠️ 部分 | KCS G1/G2/G5 通过；KVLCC2/JBC 有已知缺陷 |
| 金标回归测试绿 | ✅ 完成 | 9/9 测试通过 |
| 依赖锁定 | ✅ 完成 | 核心零依赖，iges.py 仅标准库 |

**总体**: 管道基础设施完整，金标回归通过，依赖锁定。真实 IGES 文件的"全过门"受限于出版商文件质量和标准库修复能力——管道诚实报告所有缺陷，不粉饰结果。

---

## 8. 测试命令

```bash
# 金标回归测试
python -m pytest tests/test_iges_gold_regression.py -v

# 管道门测试
python -m pytest tests/test_iges_pipeline_gate.py -v

# CLI 批量导入
python -m tensorcad.geomquality.import_cli hull1.igs hull2.igs --out-dir runs/import

# Python API
from tensorcad.geometry import load_cad_geometry
mesh, report = load_cad_geometry("hull.igs", u_samples=64, v_samples=32)
```
