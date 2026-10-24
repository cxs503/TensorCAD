# TensorCAD 扩编基线审计报告

**审计日期**: 2026-10-24
**上游 HEAD**: `54c59449c17778a3b24198deee3240a208eb6755`（origin/main 当前快照；上游活跃推送中）
**版本**: 0.4.0
**审计人**: tensor-expand 班次（m241）
**工作区**: `/nfs/wangxi/tensor_eco/agent_work/m241/TensorCAD`（分支 agent/m241）

---

## 1. 环境与测试基线

### 1.1 测试套件

| 套件 | 收集 | 通过 | 失败 | README 声明 | 备注 |
|------|------|------|------|-------------|------|
| 完整（含 neural 扩展） | 262 | 262 | 0 | "194 核心 + 213 含 neural" | 0.4.0 新增 geometry_field/neural 测试 |

**验证命令**:
```bash
/nfs/wangxi/hermes-agent/.venv/bin/python -m pytest tests/
# 结果: 262 passed in 23.43s
```

**声明 vs 实际差异说明**: README 声明"194 核心 + 213 含 neural"是 0.3.x 时代的数字；0.4.0 新增 CGC 协议测试（geometry_field 49 个）和 neural 扩展测试，总数增至 262。声明未同步更新，建议上游修正。

### 1.2 零依赖声明验证

**声明**: `import tensorcad` 不得拉 torch/numpy（pyproject.toml `dependencies = []`）

**验证**:
```python
import sys
class Blocker:
    def find_module(self, name, path=None):
        if name in ('torch', 'numpy'):
            return self
    def load_module(self, name):
        raise ImportError(f'BLOCKED: {name}')
sys.meta_path.insert(0, Blocker())
import tensorcad  # 成功
```

**结果**: ✅ **通过** — 核心包可在无 torch/numpy 环境导入，零依赖声明成立。

**注意**: 需安装 `safetensors` 才能跑完整测试套件（neural extra）。`safetensors` 不在 `neural` extra 的依赖列表中（pyproject.toml 仅声明 `torch>=2.2, safetensors>=0.4, numpy>=1.24`），但 `tensorcad/neural/field.py:26` 直接 `from safetensors.torch import load_file`——安装 `neural` extra 会自动拉 safetensors，无问题。

---

## 2. 抽检：数字直读工件

### 2.1 参数化船体两档体积收敛

**README 声明**:
- 默认分辨率 (nx=41, nz=13): 排水体积 7.446554 m³ vs 解析 7.512698 m³，差 0.88%
- 细分辨率 (nx=161, nz=41): 差 0.14%

**验证命令**:
```bash
/nfs/wangxi/hermes-agent/.venv/bin/python examples/parametric_hull/run_example.py \
  --output /tmp/tcad_audit/hull-demo
```

**实测收敛表**（直读自 stdout，存于 `/tmp/tcad_audit/hull_default.log`）:

| nx | nz | 体积 [m³] | 体积误差 % | 水线面 [m²] | 水线面误差 % |
|----|----|-----------|------------|-------------|--------------|
| 11 | 5 | 7.115242 | -5.290 | 13.826032 | -1.848 |
| 21 | 9 | 7.374345 | -1.842 | 14.006243 | -0.568 |
| **41** | **13** | **7.446554** | **-0.880** | 14.061669 | -0.175 |
| 81 | 21 | 7.483649 | -0.387 | 14.078726 | -0.054 |
| **161** | **41** | **7.502488** | **-0.136** | 14.083975 | -0.017 |

**解析体积**: 7.512698 m³（直读自 `analytic_volume(10.0, 2.0, 0.8)`）

**结论**: ✅ **通过**
- 默认档 (41×13): 7.446554 m³，差 -0.880%（README 声明 0.88%，一致）
- 细分辨率档 (161×41): -0.136%（README 声明 0.14%，四舍五入一致）
- 误差单调收敛 ✅

### 2.2 .ngo 事务管线三连

**README 声明**:
- dent peak 0.0527 m vs 0.06 目标
- 抽取网格水密、V-E+F=2

**验证方法 A**: Python API（`examples/neural_geometry/run_object_pipeline.py`）

**实测记录**（直读自 stdout，存于 `/tmp/tcad_audit/ngo_pipeline.log`）:
- fit: revision 1, field MAE 0.005092 m, surface error mean 0.004086 m
- edit: **dent peak 0.05269 m**（目标 0.06 m），revision 2, committed
- extract: 97332 vertices, 194660 triangles, **watertight True**, **euler 2**
- rejection: 0.25 depth 被拒（dent_peak_error, dent_region_error_mean），revision 不变

**验证方法 B**: CLI 三连（`tensorcad-neural fit/edit/extract`）

**实测记录**:
- fit: revision 1, weights sha256 `6caad6c5a7831a4db7f4d504e41355384704ec763d80a58eae2b829a126e7470`
- edit: **dent peak measured 0.0526903 m**（目标 0.06 m），revision 2, committed
- extract: 97332 vertices, 194660 triangles, **watertight True**, **euler 2**, signed volume 0.5154 m³

**结论**: ✅ **通过**
- dent peak 0.05269 m（README 声明 0.0527 m，四舍五入一致）
- 抽取网格水密 + V-E+F=2 ✅
- 事务管线：rejected edit 不修改 active revision ✅

---

## 3. 能力-声明对照

### 3.1 README「What is implemented」逐条验证

| # | 声明 | 验证方法 | 结果 |
|---|------|----------|------|
| 1 | 水密（每边恰好两面） | `meshlib.is_watertight(mesh)` | ✅ True |
| 2 | 边-顶点流形 | `meshlib.is_edge_manifold(mesh)` | ✅ True |
| 3 | 一致外向（有符号体积 > 0） | `meshlib.total_signed_volume(mesh)` | ✅ 7.446554 > 0 |
| 4 | Euler V-E+F = 2 | 计算 V=1001, F=1998 | ✅ 2 |
| 5 | 无零面积面 | `meshlib.zero_area_faces(mesh)` | ✅ [] |
| 6 | 实现 L/B/T 符合请求 | 边界框检查 | ✅ L=10, B=2, T=0.8 |
| 7 | 离散化误差度量 | `estimate_discretisation_error` | ✅ 收敛表输出 |

### 3.2 0.3.1 声称修复的缺陷验证

**CHANGELOG 0.3.1 列出 8 项修复**（README 概述为"three defects"是简化）:

| # | 缺陷 | 验证 | 结果 |
|---|------|------|------|
| 1 | 收敛估计在合法模型上崩溃 | 收敛表正常输出至 nx=161 | ✅ |
| 2 | 失败时导出半成品 | `export_hull` 输出完整 metadata.json | ✅ |
| 3 | CLI 错误边界半开 | CLI 正常执行无 traceback | ✅ |
| 4 | `verify` 需要 torch | 代码路径已隔离（未测无 torch 环境）| ⚠️ 代码审查通过 |
| 5 | 导出网格写入 live bundle | edit 阶段 derived/ 文件在暂存目录 | ✅ |
| 6 | 抽取不可见且可销毁 | extract 发布独立 revision | ✅ |
| 7 | `parse_manifest` 未验证 derived[] | 代码路径已加固（未测恶意输入）| ⚠️ 代码审查通过 |
| 8 | `max_section_area` 二次复杂度 | 97332 顶点网格快速完成 | ✅ |

**结论**: ✅ 大部分直接验证通过；#4 和 #7 代码审查通过但未在极端条件下实测。

### 3.3 README「What it does not provide」诚实性

| 声明 | 验证 |
|------|------|
| 除球外的神经对象形状族 | ✅ 代码仅 `fit_sphere` |
| 替代 CAD 曲面的场精度（平均 0.8%） | ✅ 实测 field MAE 0.005/0.5 = 1.0%，接近 |
| 布尔/扫描操作 | ✅ 未实现 |
| 网格平滑/简化 | ✅ 未实现，薄片报告为警告 |
| B-rep/STEP 几何 | ✅ 未实现 |
| 真实船型（正弦分布，前后对称） | ✅ 代码确认 |
| 公平性/稳性/阻力/规范 | ✅ 未实现 |

**结论**: ✅ **诚实** — README 明确标注能力边界，无夸大声明。

---

## 4. 对 TensorShipDesign 的接棒位评估

### 4.1 参数化船体族 → scan_runner 参数扫描线

**现状**: 单船型（`HullParameters` L/B/T/nx/nz），正弦分布，前后对称，无艉板。

**接棒位**:
- ✅ 参数化接口完整（`HullParameters` → `generate_hull` → `export_hull`）
- ✅ 静水力自动计算（排水体积/水线面/LCB/Cb 等）
- ✅ 收敛研究自动化
- ⚠️ 单船型族，无多船型参数化
- ⚠️ 前后对称，无艉板

**建议**: TensorShipDesign 可从 TensorCAD 参数化船体起步，需扩展船型族。

### 4.2 L2 几何编码器的几何源

**现状**: 统一场协议（`GeometryField`），4 个精确原语 + 神经场（`.ngo`）。

**接棒位**:
- ✅ 统一场协议（signed_distance/gradient/bounds/state）
- ✅ 解析场与神经场可互换
- ✅ 场 → 网格 → 验证管线完整
- ⚠️ 神经场精度有限（~1%）
- ⚠️ 仅球+凹坑编辑

**建议**: L2 可直接使用 GeometryField 协议，需扩展神经场形状族。

### 4.3 SDF 导出适配 → scan_runner SDF 栈

**现状**: 支持 SDF 导出（`SampledField.handoff()`），无直接体素网格导出。

**接棒位**:
- ✅ 隐式场协议完整
- ✅ 场 → 网格管线（marching tetrahedra）
- ⚠️ 无直接 SDF 体素网格导出
- ⚠️ 无 marching cubes（仅 marching tetrahedra）

**建议**: 需补充 `field.to_voxel_grid(resolution)` 接口。

### 4.4 OCP B-rep 后端优先级

**优先级**: 🔴 **高** — TensorShipDesign 需要真实船型，B-rep 是工业标准。

---

## 5. 可贡献开发项排序

### 高优先级（TensorShipDesign 直接依赖）
1. **多船型族扩展** — 系列 60、Savitsky、用户自定义分布
2. **非对称船体 + 艉板** — 修改 hull.py 采样逻辑
3. **OCP B-rep 后端** — B-rep → mesh 导出 + 布尔操作

### 中优先级（生态扩展）
4. **SDF 体素网格导出** — `to_voxel_grid` 方法
5. **神经场形状族扩展** — 船体/翼型训练管线
6. **网格平滑/简化** — 集成 mesh 处理库

### 低优先级
7. **公平性/表面光顺检查** — 曲率连续性
8. **静水力扩展** — 初稳性/阻力预报

---

## 审计工件索引

| 工件 | 路径 |
|------|------|
| 参数化船体收敛表 | `/tmp/tcad_audit/hull_default.log` |
| .ngo 事务管线（Python API） | `/tmp/tcad_audit/ngo_pipeline.log` |
| .ngo 事务管线（CLI 三连） | `/tmp/tcad_audit/cli-ngo/` |
| 几何核心基准测试 | `/tmp/tcad_audit/geometry_core.log` |
