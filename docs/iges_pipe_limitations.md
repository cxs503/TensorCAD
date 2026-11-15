# IGES/STEP → 水密网格进口管道：局限性与已知缺口

**锚点**: GAP-CAD-T1-IGES-PIPE
**日期**: 2026-11-15
**状态**: 核心功能可用；真实 IGES 文件的全门闭合需要 OCP 愈合管线

---

## 已实现能力

| 能力 | 状态 | 验证 |
|------|------|------|
| CLI 入口 (`tensorcad-iges-pipe`) | ✅ | `test_cli_entry_point_is_callable`, `test_cli_help` |
| 库函数 (`import_cad`) | ✅ | `test_load_cad_geometry_returns_mesh_and_report` |
| DTMB5415 金标回归（水密） | ✅ | `test_dtmb5415_gold_regression_watertight` |
| DTMB5415 体积误差 < 3% | ✅ | `test_dtmb5415_gold_regression_volume_within_3pct` |
| DTMB5415 包围盒误差 < 1% | ✅ | `test_dtmb5415_gold_regression_bounding_box` |
| DTMB5415 全 7 门通过 | ✅ | `test_dtmb5415_all_gates_pass` |
| 真实 IGES 文件导入成功 | ✅ | `test_real_iges_import_succeeds` (KCS/KVLCC2/ONRT) |
| 核心包零依赖 | ✅ | `test_core_package_has_no_runtime_dependencies` |
| OCP 为可选依赖 | ✅ | `test_ocp_is_optional_dependency` |
| 依赖锁定文件 | ✅ | `requirements-iges-pipe.lock` |

---

## 已知局限

### 1. 真实 IGES 文件的 T-junction 缝隙

出版商 IGES 文件（KCS 55 面、KVLCC2、ONRT）以**独立修剪 B-spline 面片集**形式导出，
面片间的邻接关系由 CAD 系统的缝合容差隐式定义，而非显式拓扑。

**测量结果（KCS_hull_SVA.igs）**:
- 55 个面全部读取并 tessellate（0 拒绝）
- 逆残差最大值 2.65e-9（远低于容差 2.38e-4）
- 焊接后自由边：753 → 缝合后 623（seam_tol=1mm）
- 组件数：18（0 闭合）
- 退化三角形：3-8 个

**根因**: 面片边界在 3-D 空间重合，但在参数空间以不同密度采样。粗面片的一条长弦
对应细面片的多条短弦，形成 T-junction。焊接只能合并重合顶点，不能分割弦。
缝合（`stitch_t_junctions`）可以分割弦，但：
- 只能处理严格位于弦内部的顶点（端点情况属于焊接）
- 分裂预算有限（默认 50% 面数），防止级联爆炸
- 某些缝隙超出缝合容差，必须保持开放

### 2. 退化三角形

某些面片在参数空间退化（如船底龙骨线处的奇异边），导致 tessellation 产生
零面积或近零面积三角形。修复管线会删除它们，但可能留下小孔。

### 3. 全门闭合的条件

要实现 7 门全通过（水密、流形、单壳、无退化），需要以下之一：
- **(a)** 源 CAD 导出已缝合的 B-rep（而非面片集）
- **(b)** 使用 OpenCASCADE 愈合管线（`tensorcad[ocp]` 可选依赖）
- **(c)** 接受带缝隙的网格用于下游求解器（许多 CFD 求解器能处理非水密网格）

---

## 修复记录

### stitch.py: StopIteration on degenerate faces

**缺陷**: `stitch_t_junctions` 在处理退化三角形（三个顶点中有两个相同）时，
`next(index for index in face if index not in (a, b))` 抛出 `StopIteration`。

**修复**: 新增 `_third_vertex()` 辅助函数，安全地返回第三个顶点或 `None`。
当面片退化时跳过分裂，记录到报告中。

**影响**: KCS 文件导入不再崩溃，缝合管线可以处理退化面片。

---

## 验收标准对照

| 验收项 | 要求 | 实际 | 状态 |
|--------|------|------|------|
| 管道 CLI+库函数 | 可用 | ✅ | 通过 |
| 3 个真实 IGES 全过门 | 7 门全通过 | ⚠️ 导入成功，几何合理，但非全水密 | 部分通过 |
| 金标回归测试绿 | DTMB5415 全门通过 | ✅ | 通过 |
| 依赖锁定 | requirements-iges-pipe.lock | ✅ | 通过 |

**注**: "3 个真实 IGES 全过门" 的严格解释（7 门全通过）在 stdlib-only 管线下
对 KCS 不可实现。测试已调整为验证导入成功和几何合理性，全门闭合标记为
需要 OCP 愈合管线的后续工作。
