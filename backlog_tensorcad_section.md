## TensorCAD
> 首波审计：m241 @ 2026-10-24，上游 54c59449 (v0.4.0)
> 审计报告：agent_work/m241/TensorCAD/expand_audit_tensorcad.md
> 测试基线：262/262 PASS（23.43s）；README 声明"194 核心 + 213 含 neural"为 0.3.x 旧数字，0.4.0 增至 262
> 零依赖验证：✅ import tensorcad 在 torch/numpy 被 meta_path 拦截时仍可导入
> 案例抽检：
>   - 参数化船体两档收敛：默认 (41×13) 7.446554 m³ / 解析 7.512698 m³ = -0.880%；细 (161×41) -0.136%
>   - .ngo 事务管线三连（Python API + CLI）：dent peak 0.05269 m / 目标 0.06 m；抽取 mesh 水密 + V-E+F=2
>   - 几何核心基准：4 个精确原语 + 神经场对比，all checks passed
> 0.3.1 声称修复 8 项缺陷：6 项直接验证通过，2 项代码审查通过但未极端实测

### 当前能力（已验证）
- 参数化船体：HullParameters(L/B/T/nx/nz) → generate_hull → export_hull（OBJ/CSV/JSON）
- 闭固体保证：水密 + 边流形 + 外向 + V-E+F=2 + 无零面积面，全部检查并在失败时报错
- 静水力：排水体积/水线面/LCB/LCF/VCB/Cb/Cp，网格积分
- 收敛研究：estimate_discretisation_error 自动细化网格，超限报告 skipped
- 连续几何协议（CGC v0.1）：GeometryField 统一 signed_distance/gradient/bounds/state
- 4 个精确原语：Plane/Sphere/Cylinder/Ellipsoid，双精度零误差
- 神经几何对象（.ngo）：fit/edit/extract 事务管线，rejected edit 不修改 active revision
- 场→网格：to_mesh(field, resolution) → marching tetrahedra → 验证管线
- 标准库-only 核心：零依赖声明成立

### 缺口（README 自声明 + 审计确认）
- ~~**单船型族**~~ → **部分已补（0.20.0）**：新增 `tensorcad.wigley`（抛物型 Wigley，指数 n/m 可调），带精确解析静水力（Cb=4/9、Cp=Cm=Cwp=2/3、VCB=5T/8）与收敛证据 `docs/wigley.md`；仍缺系列 60 站型表 / Savitsky / 带艉板非对称族
- **神经场形状族**：仅球+凹坑编辑，无船体/翼型/通用形状
- **无 B-rep/STEP**：设计原则 4 声明"evaluate OCP"但未实现
- **无布尔/扫描操作**
- **无网格平滑/简化**：2-3% 薄片报告为警告，不修复
- ~~**无 SDF 体素网格导出**~~ → **0.10.0 已补**：`to_voxel_grid`/`VoxelGrid`（docs/voxel_export.md；SUBOFF 体素体积 0.698982 m³ vs 文献 0.699211 = −0.033 %）
- **无公平性/稳性/阻力/规范合规** → 公平性 0.12 已补（docs/fairness.md）；**初稳性 KM/GM 0.14 已补**（tensorcad/stability.py + docs/stability.md：水线面惯性矩 I_T/I_L，箱体基准 2.2e-16、半球 KM=R −0.03%、船体 GM 默认网格 +1.99%/161×41 +0.35%）；大倾角 GZ 已补（0.21.0 gz.py/righting_arm.py）；**IS Code 2008 Part A 2.2 完整稳性判据 0.22.x 已补**（tensorcad/intact_stability.py + docs/intact_stability.md：6 条判据逐条判定、undetermined 三态、面积规则二阶收敛 2.000、箱驳面积 vs 墙侧闭式 1.8e-4 m·rad、KG 负控 7e-8）；自由液面修正已有 free_surface_i_x 平降，§2.3 风倾复正、阻力、第二代失效模式仍缺
- **神经场精度有限**：平均误差 ~1% 半径，最坏 ~2%

### 对 eco 价值排序（与 TensorShipDesign 接棒位）
| 优先级 | 项目 | eco 价值 | 接棒位 |
|---:|---|---|---|
| P1 | 多船型族扩展（Wigley 已补 0.20.0；缺系列 60/Savitsky） | 实际船舶设计 | scan_runner 参数扫描线几何源 |
| P0 | 非对称船体 + 艉板 | 真实船型 → 工程应用 | TensorShipDesign 刚需 |
| P0 | OCP B-rep 后端 | 工业标准几何 → 布尔/导出 | TensorShipDesign 核心依赖 |
| P1 | SDF 体素网格导出 | scan_runner SDF 栈对接 | L2 几何编码器几何源 |
| P1 | 神经场形状族扩展 | 船体/翼型近似 → 可微设计 | L2 几何编码器扩展 |
| P2 | 网格平滑/简化 | 改善 mesh 质量 → 求解器稳定性 | TensorBEM/TensorFEM 消费 |
| P3 | 公平性/光顺检查 | 曲率连续性 → 设计质量 | TensorShipDesign 后期 |
| P3 | 静水力扩展（稳性/阻力） | 完整船舶设计 → 规范合规 | TensorShipDesign 后期 |

### 可贡献开发项（按可实现性 × eco 价值）
1. **多船型族扩展**（P0）：在 hull.py 新增分布函数族（系列 60、Savitsky、用户自定义），参数化框架已有
2. **非对称船体 + 艉板**（P0）：修改 hull.py 采样逻辑，支持左右不对称 + 艉封板/艉板
3. **OCP B-rep 后端**（P0）：集成 OCP，实现 B-rep → mesh 导出 + 布尔操作（工作量大）
4. **SDF 体素网格导出**（P1）：新增 `to_voxel_grid(resolution)` 方法，对接 scan_runner SDF 栈
5. **神经场形状族扩展**（P1）：船体/翼型训练管线 + 编辑原语（工作量大）
6. **网格平滑/简化**（P2）：集成 mesh 处理库（如 trimesh/pyvista），修复薄片

### 与既有线的边界
- TensorLBM 归 benchmark_closure 线
- wedge/SI 归 si-campaign 线
- 本线（tensor-expand）只碰 TensorBEM/TensorSPH/TensorDG/TensorCAD 四新库
- TensorShipDesign 是 TensorCAD 下游应用，本线负责几何底座，应用层由 owner 另行派发
