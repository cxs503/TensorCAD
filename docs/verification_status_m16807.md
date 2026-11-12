# TensorCAD 验证状态综合报告

**Mission**: m16807 (selfstart generalist seat)  
**Date**: 2026-11-12  
**Purpose**: 综合审计验证覆盖度，识别下一步可行动项

## 执行摘要

TensorCAD 已达到高度成熟的验证覆盖：2791 个测试全部通过，所有主要模块均有解析解/文献/独立数值路径的验证文档。本报告显示验证状态并识别剩余缺口。

## 1. 已验证模块（全部测试通过）

### 1.1 船体型族

| 模块 | 验证基准 | 关键结果 | 状态 |
|------|---------|---------|------|
| `tensorcad.hull` | 解析体积 | 体积收敛至解析解（161×41 时 0.14%） | ✅ 完整 |
| `tensorcad.wigley` | 精确解析静水力（Cb=4/9, VCB=5T/8） | 所有系数精确匹配 | ✅ 完整 |
| `tensorcad.series60` |  published parent forms (Todd 1951) | Cwp 匹配 published 值 <0.01%（fidelity fix 后） | ✅ 完整 |
| `tensorcad.suboff` | DARPA offset table (Groves et al. 1989) | 体积 0.699211 m³ vs 文献 0.699211 m³ (0.000%) | ✅ 完整 |
| `tensorcad.planing` | 常升角船型闭式静水力 | 湿面积和体积匹配闭式 | ✅ 完整 |
| `tensorcad.asymmetric_hull` | 左右非对称 loft | 水密、流形、尺寸正确 | ✅ 完整 |

### 1.2 静水力与稳性

| 模块 | 验证基准 | 关键结果 | 状态 |
|------|---------|---------|------|
| `tensorcad.hydrostatics` | 散度定理 on mesh | 箱体精确至舍入误差 | ✅ 完整 |
| `tensorcad.stability` | 水线面矩（Green 定理） | 箱体/半球锚点匹配闭式 | ✅ 完整 |
| `tensorcad.righting_arm` | 箱形船 wall-sided 公式、圆柱 | GZ(φ) 匹配闭式 <1e-15 m | ✅ 完整 |
| `tensorcad.gz` | 大倾角稳性 with heeled displacement | 体积守恒 <1e-13 相对误差 | ✅ 完整 |
| `tensorcad.trim` | TPC/MCT 双路径 | 箱体精确；参数化船体 <0.6% | ✅ 完整 |
| `tensorcad.damage` | 损失浮力法 | 舱室进水匹配解析 | ✅ 完整 |
| `tensorcad.free_surface` | 精确倾角相关液舱修正 | 小倾角匹配一阶 | ✅ 完整 |
| `tensorcad.bonjean` | 沿船体剖面属性 | Simpson 积分对线性面积精确 | ✅ 完整 |
| `tensorcad.intact_stability` | IS Code 2008 part A | 6 条通用准则 + 天气准则 | ✅ 完整 |

### 1.3 几何与网格质量

| 模块 | 验证基准 | 关键结果 | 状态 |
|------|---------|---------|------|
| `tensorcad.intersect` | 精确有理算术 oracle | 4000 对中 0 逻辑错误 | ✅ 完整 |
| `tensorcad.nurbs` | 有理 B 样条理论 (Piegl & Tiller) | 圆柱面积精确至 1e-15；球体积 O(h²) 收敛 | ✅ 完整 |
| `tensorcad.voxel` | SDF 体素化 | SUBOFF 体积 0.698982 m³ vs 文献 0.699211 m³ (-0.033%) | ✅ 完整 |
| `tensorcad.curvature` | 微分几何 | 原语主曲率匹配解析 | ✅ 完整 |
| `tensorcad.mesh_conditioning` | 网格质量度量 | 薄片检测和报告 | ✅ 完整 |
| `tensorcad.high_order_geometry` | 解析曲面积分 | 球/柱/环面积误差 <1e-12 | ✅ 完整 |

### 1.4 专用模块

| 模块 | 验证基准 | 关键结果 | 状态 |
|------|---------|---------|------|
| `tensorcad.icechannel` | 球面三角学 | 1061 km 尺度走廊畸变 <0.5% | ✅ 完整 |
| `tensorcad.pc6` | PC6 船型参数 | L/B/Cb 门控验证 | ✅ 完整 |
| `tensorcad.statutory` | 完整稳性准则 | 准则评估和报告 | ✅ 完整 |
| `tensorcad.resistance` | ITTC-57/Schoenherr/Prandtl-Schlichting | 摩擦线匹配 <2.5%（Re≥1e7） | ✅ 完整 |
| `tensorcad.csg` | 布尔运算（并/交/差） | 6 个闭式体积匹配 <0.02% | ✅ 完整 |
| `tensorcad.volume_mesh` | 体网格生成 | 箱体精确至舍入；船体匹配面网格体积 | ✅ 完整 |
| `tensorcad.fillet` | 倒角/倒棱 | 面积变化闭式精确 | ✅ 完整 |
| `tensorcad.girder` | 剖面属性 | 矩形/工字梁精确至 1e-15 | ✅ 完整 |
| `tensorcad.shear_flow` | 薄壁剪流 | I 梁 VQ/I 精确至舍入 | ✅ 完整 |
| `tensorcad.sketch` | 2D 约束求解 | 10 个夹具 100% 收敛，位置 <1e-10 | ✅ 完整 |
| `tensorcad.design` | 设计循环 | Newton 收敛至 1e-7 相对 | ✅ 完整 |

## 2. 验证方法论

### 2.1 多路径验证

关键量由两个或更多独立路径计算：
- **体积**：mesh 散度定理 vs 解析闭式 vs tanh-sinh 求积
- **水线面矩**：Green 定理 on 水线面环 vs 站 sum
- **GZ 曲线**：heeled displacement 积分 vs wall-sided 公式（箱形船）
- **TPC/MCT**：水线面面积/矩 vs 位移有限差分

### 2.2 解析锚点

每个有闭式解的模块都对其验证：
- **Wigley**: Cb=4/9, Cp=Cm=Cwp=2/3, VCB=5T/8（精确有理数）
- **Series 60**: Published parent forms (Todd 1951) with Cwp fidelity fix
- **SUBOFF**: DARPA offset table (Groves et al. 1989)
- **Box barge**: Wall-sided formula exact to round-off
- **Sphere**: Area=4πR², Volume=4/3πR³

### 2.3 收敛研究

网格积分验证在加密下收敛：
- **参数化船体**: 体积误差 5.29% → 0.14%（nx 11 → 161）
- **Wigley**: 所有系数收敛至解析值
- **NURBS 球**: 体积 O(h²) 收敛
- **右臂**: GZ(φ) 随网格加密收敛

## 3. 剩余缺口与未来工作

### 3.1 高优先级

1. **OCP B-rep 后端** (GAP-CAD-001)
   - **缺口**: 无工业标准 B-rep 几何（STEP/IGES 导入/导出）
   - **影响**: 无法与 CAD 系统接口或执行布尔操作
   - **工作量**: 大（需要 OCP 集成、多片 NURBS、裁剪曲面）
   - **建议**: 专用 mission + CAD 内核专业知识

2. **精确曲面-曲面求交** (GAP-CAD-002)
   - **缺口**: intersect.py 只做三角-三角相交，无曲面级求交
   - **影响**: 船体/结构前处理阻塞
   - **工作量**: 大
   - **建议**: 专用 mission

3. **多片 NURBS 缝合**
   - **缺口**: 仅单片 NURBS；无裁剪片或多片模型
   - **影响**: 无法表示复杂船体型（球鼻艏、艉板角）
   - **工作量**: 中大
   - **建议**: 在现有单片基础设施上构建

### 3.2 中优先级

4. **网格平滑/简化** (GAP-CAD-009)
   - **缺口**: 有去退化面、Taubin 平滑、抽稀；但检测到自相交后无修复路径
   - **影响**: 改善网格质量 → 求解器稳定性
   - **工作量**: 中
   - **建议**: 集成网格处理库

5. **多船型族扩展**
   - **现状**: Wigley (0.20.0) 和 Series 60 已补；缺 Savitsky 计划船型
   - **影响**: 实际船舶设计
   - **工作量**: 中
   - **建议**: 在 hull.py 新增 Savitsky 分布

### 3.3 低优先级

6. **草图约束求解器扩展** (GAP-CAD-006)
   - **现状**: sketch.py 已实现基础求解器（10 个夹具验证通过）
   - **缺口**: 无 3D 约束、无不等式/尺寸方案排序、无约束优先级
   - **工作量**: 中
   - **建议**: 按需扩展

## 4. 对 TensorShipDesign 的接棒位评估

### 4.1 参数化船体族 → scan_runner 参数扫描线

**现状**:
- ✅ 参数化接口完整（HullParameters → generate_hull → export_hull）
- ✅ 静水力自动计算（排水体积/水线面/LCB/Cb 等）
- ✅ 收敛研究自动化
- ✅ 多船型族（Wigley, Series 60, SUBOFF, planing, asymmetric）
- ⚠️ 缺 Savitsky 计划船型

**建议**: TensorShipDesign 可从 TensorCAD 参数化船体起步。

### 4.2 L2 几何编码器的几何源

**现状**:
- ✅ 统一场协议（signed_distance/gradient/bounds/state）
- ✅ 解析场与神经场可互换
- ✅ 场 → 网格 → 验证管线完整
- ✅ SDF 体素网格导出（to_voxel_grid）
- ✅ CSG 布尔运算（union/intersection/difference）
- ⚠️ 神经场精度有限（~1%）
- ⚠️ 仅球+凹坑编辑

**建议**: L2 可直接使用 GeometryField 协议，需扩展神经场形状族。

### 4.3 OCP B-rep 后端优先级

**优先级**: 🔴 **高** — TensorShipDesign 需要真实船型，B-rep 是工业标准。

## 5. 结论

TensorCAD 已达到高度成熟的验证覆盖，所有主要模块均有解析解/文献/独立数值路径的验证。剩余缺口主要是大型 CAD 内核功能（B-rep、曲面求交），需要专用 mission 和 CAD 内核专业知识。

**下一步建议**:
1. 优先推进 OCP B-rep 后端（GAP-CAD-001），这是 TensorShipDesign 的核心依赖
2. 补充 Savitsky 计划船型，完善多船型族
3. 继续推进开源 CAD/几何内核技术盘点战役（FEAT-OSS-INVENTORY-CAD）

---

**验证命令**:
```bash
python3 -m pytest tests/ -v
# 结果: 2791 passed
```

**工件索引**:
- 验证状态汇总: 本文档
- 各模块详细验证: 见 docs/ 目录下各模块文档
- 测试套件: tests/（2791 个测试）
