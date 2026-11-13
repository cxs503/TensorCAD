# GAP-CAD-B3-SOLVER-INTERFACE: 统一几何库 → 三求解器端到端跑通（z-up 契约）

## ① 问题设定与基准出处
锚点要求「三求解器各一个端到端算例真跑通（坐标系契约统一 z-up）；输入=统一几何库输出」。
本仓给出 `tensorcad.solver-frame/1` 契约（x-aft / y-starboard / z-up，右手，米，adapter 不做旋转），
并把同一几何对象分别交给 **TensorBEM**（常值面元）、**TensorFEM**（Q4 Mindlin 板元）、**TensorLBM**（D2Q9/BGK）。
参考值全部是闭式解，由本仓内公式计算，不由求解器自证：
- 球体附加质量 `2/3 π ρ a³`（a=0.5 m, ρ=1000 kg/m³，球心 z=-1.5a）；体积 `4/3 π a³`、湿表面积 `4 π a²`；
- 简支方板单正弦载荷 Kirchhoff 中心挠度 `q0 a⁴/(4 π⁴ D)`，`D = E t³/(12(1-ν²))`（a=2 m, E=7e10 Pa, ν=0.3, t=0.05 m, q0=1e4 Pa，Timoshenko/Navier 解）；
- 体力驱动平面 Poiseuille 精确抛物 `u(z)=F z(H−z)/(2ν)`，`ν=(τ−1/2)/3`（τ=0.8, F=2e-5）。

## ② 方法与实现（关键代码指针）
- 契约与度量：`tensorcad/solver_interface.py` — `frame_of`（未声明坐标系直接拒绝）、`GeometryBundle`（体积/面积/边界边/欧拉数/sha256 digest）、`mesh_digest`（17 位往返）。
- 三个交接工件（写文件，不走内存对象）：`bem_panels` → `panels_<n>.json`；`fem_plate_input` → `plate_<n>x<n>.json`；`lbm_lattice` + `VoxelGrid.as_dict`（packbits + base64 占据掩膜）→ `lattice_<n>.json` / `voxels_<n>.json`。
- 端到端：`run_bem_added_mass`（TensorBEM `compute_added_mass_matrix_3d`）、`run_fem_plate`（TensorFEM `q4_mindlin_stiffness` + 本仓装配/正弦一致载荷/简支约束）、`run_lbm_channel`（TensorLBM `equilibrium`/`collide_bgk`/`macroscopic` + 本仓半程反弹 streaming）。
- 平台侧调用：`tensorcad/benchmark.py::solver_handoff_case` 成为 Case C（此前该案例自述「not implemented」）；CLI `python -m tensorcad.solver_interface`（console script `tensorcad-solvers`）。
- 2D 求解器不旋转几何，而是记录 axis map（`grid_i=+x`, `grid_j=+z`, y 不变），映射随工件落盘。

## ③ 参数表（无量纲/有量纲完整）
| 算例 | 关键参数 |
|---|---|
| BEM 球 | a=0.5 m, ρ=1000 kg/m³, 面板 16×32 与 24×48, 求积阶 6, 球心 z=-0.75 m |
| FEM 板 | a=b=2 m, E=7e10 Pa, ν=0.3, t=0.05 m, q0=1e4 Pa, D=801282 N·m, 网格 8×8 与 16×16 Q4 |
| LBM 槽道 | L=0.8 m, H=0.4 m, B=0.2 m, dx=H/n（n=8,16）, τ=0.8, ν=0.1, F=2e-5, Ma≈1.10e-02, 步数 1500/5500 |

## ④ 对照表（ref vs got，逐项误差 + 门禁）
| 求解器 | 量 | ref | got | 误差 | 门禁 |
|---|---|---|---|---|---|
| tensorbem | `enclosed_volume_16x32_m3` | 0.523599 m^3 | 0.515243 m^3 | -1.596 % | None |
| tensorbem | `enclosed_volume_24x48_m3` | 0.523599 m^3 | 0.519871 m^3 | -0.712 % | True |
| tensorbem | `wetted_area_24x48_m2` | 3.14159 m^2 | 3.13039 m^2 | -0.357 % | True |
| tensorbem | `added_mass_heave_24x48_kg` | 261.799 kg | 260.296 kg | -0.574 % | True |
| tensorbem | `added_mass_surge_24x48_kg` | 261.799 kg | 260.933 kg | -0.331 % | True |
| tensorfem | `centre_deflection_n8_m` | 0.000512478 m | 0.000509827 m | -0.517 % | True |
| tensorfem | `centre_deflection_n16_m` | 0.000512478 m | 0.000513181 m | +0.137 % | True |
| tensorlbm | `centreline_velocity_8` | 0.001575 lattice | 0.001552 lattice | -1.460 % | True |
| tensorlbm | `centreline_velocity_16` | 0.006375 lattice | 0.006352 lattice | -0.361 % | True |

门禁总计：24 条，失败 0 条（[]）。
附加：BEM 闭合残差 3.9e-17、双层面行和误差 8.4e-15；LBM 质量守恒 3.7e-13、对称性 1.5e-14、CAD→体素单元数 4096 vs 4096（精确）。

## ⑤ 收敛性/敏感性（两分辨率）
- BEM 体积误差 1.596 % → 0.712 %（面板 1.5×，O(h²)），附加质量 1.280 % → 0.574 %；
- FEM 挠度误差 -0.517 % → +0.137 %（随 P 型趋近 Kirchhoff，跨零）；
- LBM 峰值误差 1.460 % → 0.361 %（H×2，误差∝F/u_max∝1/H²，比值 0.250）；剖面**形状**与抛物差 3e-9（离散稳态本身就是抛物，剩余偏差是体力格式的常数 O(F) 偏移，已单列 `forcing_offset_lattice`）。CAD→体素对盒子精确（0 误差），对球 O(dx²)。

## ⑥ 失败模式、工件与复现
- 缺依赖 = 状态不是数字：`SolverUnavailableError` → `status=unavailable`，CLI 退出码 2，`results.json` 记录每个求解器的模块路径与版本。
- **发现（跨仓，如实记录）**：探测到所测版本的 `tensorlbm.stream` 把一个 +x 布居搬到 x−1（与自身 `d2q9.C` 反号），即逆风被动输运；本模块因此自带半程反弹 streaming，只用 TensorLBM 的平衡态/碰撞/宏观量，并把 `stream_probe` 写进工件（不静默绕过）。
- 局限：LBM 为 x-z 平面 2D，y 不变（不宣称阻力预报）；BEM 算例是球（有闭式解），不是船体；体素器为射线奇偶栅格化（盒对齐精确、曲面 O(dx²)），非 cut-cell/VOF。
- 工件：`artifacts/m17476/solver_interface/{results.json, bem/, fem/, lbm/}`（含 panels/plate/lattice/voxels/manifest/run）。
- 复现：`PYTHONPATH=. python -m tensorcad.solver_interface --out artifacts/m17476/solver_interface`；测试 `PYTHONPATH=. python -m pytest tests/test_solver_interface.py -q`（11 passed）。

## 误差判定（是否 <3%）
关键量逐项：BEM 24×48 体积 -0.712 %、面积 -0.357 %、附加质量(垂荡) -0.574 %、(纵荡) -0.331 %；FEM 16×16 挠度 +0.137 %（8×8 为 -0.517 %）；LBM 16 单元峰值 -0.361 %（8 单元 -1.460 %），剖面形状 3.8e-9。
**结论：所有在细网格上报的关键量误差均 <3%（最大 1.46 % 出现在 LBM 最粗网格 8 单元，已作为收敛阶梯的上级列出）。**

