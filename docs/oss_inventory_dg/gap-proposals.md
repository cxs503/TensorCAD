# TensorDG 缺口锚点提议（GAP-DG-*，本波 7 条）

> 程序化生成，勿手改。源：`registry.json` 的 `gaps` 段（生成时间 2026-11-13T01:51:06+00:00）。
> 现状证据一律在本仓可 grep：`TensorDG` @ `281ba5c6cf6777ced440b8df756ad4402287d696`，scope=`tensordg/`。
> 验收判据必须是机器可裁的数值门；外部锚必须落在文献/解析解，不得用本仓实现自证。

## 排序总表

| # | 锚点 | 优先级 | 参考实现（Tier） | 工作量 | 现状证据 | 判据要点 |
|---|---|---|---|---|---|---|
| 1 | `GAP-DG-001` | P0 | `nektar`(T1), `ngsolve`(T1), `mfem`(T1), `dealii`(T1) | L | ✅ absent | 在 P1/P2 下用两套非结构三角剖分（每套 >=2000 单元）解线性对流与 Poisson(SIPG)：L2 误差对解析解的收敛阶不低于 p+0.8（p=1 时 >=1.8，p=2 时 >=2.8），同一算例在三角与四边形网格上的 p=2 误差比值 <=1.3，周期边界下的总质量守恒相对误差 <=1e-12。 |
| 2 | `GAP-DG-002` | P0 | `p4est`(T2), `amrex`(T2), `philip`(T1), `trixi`(T1) | L | ✅ absent | 2D Euler（前台阶或双马赫反射）在固定自由度预算（<=均匀基线 N_dof 的 25%）下，自适应解的密度 L1 误差 <=均匀加密到 4*N_dof 解的 60% |
| 3 | `GAP-DG-003` | P1 | `trixi`(T1), `flexi`(T1), `exadg`(T1) | M | ✅ absent | 1D/2D Euler（Sod + Shu-Osher）：分裂型算子配熵守恒通量后，每步离散总熵变化 <=1e-13（符号一致，即非增） |
| 4 | `GAP-DG-004` | P1 | `trixi`(T1), `flexi`(T1), `gkeyll`(T1), `mfem`(T1) | M | ✅ partial | 2D Euler 双稀疏波（近真空，rho 与 p 起值 1e-8，p=3）：全程 rho,p >=0 且无 NaN |
| 5 | `GAP-DG-005` | P1 | `nek5000`(T1), `nekr`(T1), `pyfr`(T1), `libceed`(T2) | L | ✅ absent | 2D Euler（p=3，512x512 单元）在 4 个 MPI 秩与 1 秩上：密度 L2 差 <=1e-10（同一离散解） |
| 6 | `GAP-DG-006` | P1 | `exadg`(T1), `petsc`(T2), `sundials`(T2), `dealii`(T1) | L | ✅ partial | 2D 对流扩散（Pe=1e3，p=3）推到稳态：Newton-Krylov + p-多重网格预条件的线性迭代次数从 64x64 到 256x256 增长 <=1.5 倍，达到同一收敛判据所需时间步数 <=显式 SSP-RK3 的 1/10 |
| 7 | `GAP-DG-007` | P2 | `philip`(T1), `firedrake`(T2), `dealii`(T1) | M | ✅ partial | 2D Poisson(SIPG) 与 2D Euler 各 2 个算例上，DWR 自适应 5 轮后目标泛函（阻力/线积分）误差 <=同自由度均匀加密解的 50%，效应性指数（估计/真误差）落在 0.8~1.25 区间（4 个算例全部满足）。 |

## 逐条详述

### `GAP-DG-001` — 非结构单纯形（三角/四面体）单元与网格拓扑

- **优先级**：P0；**工作量**：L；**功能点**：`F07-simplex`, `F01-basis`
- **现状（TensorDG）**：全部 DG 算子建在张量积（Cartesian/周期）网格上：euler_dg_2d 的文档串明说 periodic Cartesian mesh，advection_1d/2d/3d、euler_dg_1d/2d、swe_1d/2d、poisson_2d(SIPG) 同构；仓内无三角/四面体单元、无面邻接拓扑、无外部网格读入（gmsh）。
- **现状证据（grep @ 281ba5c6cf67）**：
  - `triangle|tetra|simplex|gmsh|unstructured` 期望=absent → 零命中 0 处 ✔
- **参考实现**：
  - `nektar/nektar` (T1, failed) — 本次未取到该功能点的路径命中
  - `NGSolve/ngsolve` (T1, ok) — F07-simplex×1 py_tutorials/scattering_triangle.py；F01-basis×1 fem/polynomials.hpp
  - `mfem/mfem` (T1, ok) — F07-simplex×5 mesh/gmsh.cpp；F01-basis×6 fem/ceed/interface/basis.cpp
  - `dealii/dealii` (T1, ok) — F07-simplex×11 bundled/boost-1.84.0/include/boost/geometry/strategies/cartesian/side_by_triangle.hpp；F01-basis×141 bundled/boost-1.84.0/include/boost/math/interpolators/bezier_polynomial.hpp
  - `FEniCS/dolfinx` (T2, ok) — F07-simplex×1 python/dolfinx/io/gmsh.py
  - `dune-community/dune-fem` (T1, failed) — 本次未取到该功能点的路径命中
- **验收判据（机器可裁）**：在 P1/P2 下用两套非结构三角剖分（每套 >=2000 单元）解线性对流与 Poisson(SIPG)：L2 误差对解析解的收敛阶不低于 p+0.8（p=1 时 >=1.8，p=2 时 >=2.8），同一算例在三角与四边形网格上的 p=2 误差比值 <=1.3，周期边界下的总质量守恒相对误差 <=1e-12。
- **上游要点**：curated：nektar 与 hp3d 提供四面体/三角单元的基函数-求积-面拓扑三件套；ndg 与 mfem 提供可直接对照的 triangle 算子数值结果。本波未读源码。
- **外部锚**：Hesthaven & Warburton, Nodal Discontinuous Galerkin Methods, Springer 2008, Ch.6（三角单元算子与收敛阶）；Cockburn & Shu, SIAM J. Numer. Anal. 35(6) 1998。
- **影响**：船体曲面与结构前处理输出的是三角面片/四面体网格；没有单纯形单元，TensorDG 无法吃真实 CAD 网格，只能做规则区间的教科书算例。

### `GAP-DG-002` — h/p 自适应（悬挂节点界面 + p 变更 + 误差指示）

- **优先级**：P0；**工作量**：L；**功能点**：`F06-adaptivity`, `F12-error-estimation`
- **现状（TensorDG）**：无任何 AMR/hp 机制：没有四/八叉森林、没有悬挂节点（非协调）界面通量、没有单元级 p 变更记录；加密只能整网格重跑（convergence.py 是对均匀加密层做阶次拟合的工具，不是自适应）。
- **现状证据（grep @ 281ba5c6cf67）**：
  - `hanging_node|hanging_nodes|p4est|octree|quadtree|p_adapt|hp_adapt|nonconforming|\bAMR\b` 期望=absent → 零命中 0 处 ✔
- **参考实现**：
  - `cburstedde/p4est` (T2, ok) — F06-adaptivity×49 src/p4est.c
  - `AMReX-Codes/amrex` (T2, ok) — F06-adaptivity×766 Src/Amr/AMReX_Amr.H
  - `McGillCompAero/PHiLiP` (T1, ok) — F06-adaptivity×46 src/grid_refinement/field.cpp；F12-error-estimation×4 src/functional/adjoint.cpp
  - `trixi-framework/Trixi.jl` (T1, ok) — F06-adaptivity×39 src/auxiliary/p4est.jl
  - `parthenon-hpc-lab/parthenon` (T2, ok) — F06-adaptivity×12 src/amr_criteria/amr_criteria.cpp
  - `NGSolve/ngsolve` (T1, ok) — F06-adaptivity×1 py_tutorials/adaptive.py
- **验收判据（机器可裁）**：2D Euler（前台阶或双马赫反射）在固定自由度预算（<=均匀基线 N_dof 的 25%）下，自适应解的密度 L1 误差 <=均匀加密到 4*N_dof 解的 60%；跨悬挂节点界面总质量漂移 <=1e-12；p 变更后密度最小值 >=同算例均匀解的 0.99（无 Gibbs 超调）。
- **上游要点**：curated：p4est 提供 2:1 平衡的森林与悬挂节点接口构造；philip/trixi 提供 hp 判据（截断误差指示）设计；amrex/parthenon 提供块结构 AMR 的守恒回流（reflux）设计。本波未读源码。
- **外部锚**：Burstedde, Wilcox & Ghattas, SIAM J. Sci. Comput. 33(2) 2011（p4est）；Ceze & Fidkowski, J. Comput. Phys. 230 (2011)（hp 自适应判据）；Berger & Colella, J. Comput. Phys. 82 (1989)（AMR 守恒回流）。
- **影响**：船舶水动力算例里激波/自由面与边界层尺度差 2~3 个量级，全域 p 提升的代价不可支付；自适应是三维实船算例能否跑得动的前提。

### `GAP-DG-003` — 熵稳定/分裂型（split-form、两点通量）离散

- **优先级**：P1；**工作量**：M；**功能点**：`F05-entropy-stable`
- **现状（TensorDG）**：没有分裂型（skew-symmetric / split-form）算子，也没有两点熵守恒通量与熵稳定界面耗散；entropy 一词只出现在算例名（Shu-Osher 激波-熵波相互作用）与熵预算脚本里，算子本身是节点型 HLLC + TVB 限界。
- **现状证据（grep @ 281ba5c6cf67）**：
  - `split_form|split-form|entropy_stable|entropy-stable|kuscu|entropy_conserv|two-point_flux|two-point flux` 期望=absent → 零命中 0 处 ✔
- **参考实现**：
  - `trixi-framework/Trixi.jl` (T1, ok) — F05-entropy-stable×9 src/callbacks_stage/entropy_bounded_limiter.jl
  - `flexi-framework/flexi` (T1, ok) — 本次未取到该功能点的路径命中
  - `exadg/exadg` (T1, ok) — 本次未取到该功能点的路径命中
- **验收判据（机器可裁）**：1D/2D Euler（Sod + Shu-Osher）：分裂型算子配熵守恒通量后，每步离散总熵变化 <=1e-13（符号一致，即非增）；动量与能量守恒相对误差 <=1e-12；Sod 密度 L1 误差相对本仓现有 HLLC 节点基线 <=1.05 倍（精度不得退化）。
- **上游要点**：curated：trixi 的 volume/surface 分裂型内核与熵守恒通量是唯一完整的开源参考实现；flexi 给出工程尺度上的算子组织方式。本波未读源码。
- **外部锚**：Fisher & Carpenter, J. Comput. Phys. 252 (2013)（entropy stable staggered grid / split form）；Gassner, Winters & Kopriva, J. Comput. Phys. 327 (2016)（split-form DGSEM 的熵稳定性）；Chan, J. Comput. Phys. 378 (2019)（discretely entropy stable DG）。
- **影响**：耐波/操纵这类必须长时积分的时域算例，非熵稳定离散会带来缓慢的能量漂移，误差天花板由格式而非物理决定。

### `GAP-DG-004` — 二维保界（Zhang-Shu 缩放 / 子单元限界）

- **优先级**：P1；**工作量**：M；**功能点**：`F04-bound-preserving`, `F03-limiting`
- **现状（TensorDG）**：保界只做到一半：1-D Euler 有 Zhang-Shu 可容许缩放限制器 _positivity_scale（tensordg/euler_dg_1d.py:98，theta 缩放保住单元均值）；2-D 起只有事后过滤 _positivity_filter（tensordg/euler_dg_2d.py:659：不可容许的整格直接塌成单元均值——均值守恒但该单元退化为一阶，不是 theta 缩放，也没有子单元信息）。无子单元（subcell）FV 混合、无 MOOD/a-posteriori 限界、无凸限界；SWE/NS 模块没有保界路径。
- **现状证据（grep @ 281ba5c6cf67）**：
  - `_positivity_scale|zhang.?shu` 期望=present → 命中 6 处 ✔
    - `tensordg/euler_dg_1d.py:48: """Zhang-Shu admissible scaling factors ``theta`` for one modal cell block.`
    - `tensordg/euler_dg_1d.py:98: def _positivity_scale(coeff, check_points, eta=1e-10):`
    - `tensordg/euler_dg_1d.py:99: """Zhang-Shu scaling limiter: admissible cell, maximum retained high order.`
  - `subcell|\bmood\b|convex_limit|convex limiting` 期望=absent → 零命中 0 处 ✔
- **参考实现**：
  - `trixi-framework/Trixi.jl` (T1, ok) — F04-bound-preserving×25 ext/TrixiConvexECOSExt.jl；F03-limiting×18 src/callbacks_stage/entropy_bounded_limiter.jl
  - `flexi-framework/flexi` (T1, ok) — F03-limiting×2 src/filter/pplimiter.f90
  - `gkeyllorg/gkeyll` (T1, ok) — F04-bound-preserving×14 gyrokinetic/apps/gk_species_positivity.c；F03-limiting×543 gyrokinetic/apps/gk_neut_species_moment.c
  - `mfem/mfem` (T1, ok) — F03-limiting×10 fem/tmop/assemble/diag2_limit.cpp
- **验收判据（机器可裁）**：2D Euler 双稀疏波（近真空，rho 与 p 起值 1e-8，p=3）：全程 rho,p >=0 且无 NaN；总质量与总能量相对误差 <=1e-12；把 2D 的事后过滤换成 Zhang-Shu theta 缩放后，触发单元-步比例（整格塌陷计数/总单元-步）从现状 >10% 降到 <=1%；等熵涡（p=3）的 L2 误差退化 <=10%。
- **上游要点**：curated：trixi 的子单元激波捕捉（指示器 + subcell 混合）与 mfem 的公开 subcell limiter 算例是两条可搬路线；Zhang-Shu 缩放本仓 1-D 已有，2-D 化是直接增量。本波未读上游源码。
- **外部锚**：Zhang & Shu, Proc. R. Soc. A 466 (2010)（DG 保界限制器）；Dumbser, Zanotti, Loubère & Cassulli, J. Comput. Phys. 278 (2014)（a posteriori subcell limiting）；Hennemann, Rueda-Ramírez, Hindenlang & Gassner, J. Comput. Phys. 424 (2021)（子单元激波捕捉）。
- **影响**：二维近真空/空化工况是自由面与螺旋桨空泡的常态；整格塌陷会在局部把精度降成一阶，并以未受控的方式污染阻力/升力积分量。

### `GAP-DG-005` — 分布式内存与加速器后端（MPI / GPU）

- **优先级**：P1；**工作量**：L；**功能点**：`F11-parallel`
- **现状（TensorDG）**：纯单进程 NumPy/Torch CPU 实现：tensordg 内无 MPI 通信、无 CUDA/Kokkos/OpenMP 内核、无分区与 halo 交换；并行化只能靠多进程各算各的算例。
- **现状证据（grep @ 281ba5c6cf67）**：
  - `mpi4py|cuda|kokkos|cupy|numba|gpu|openacc|openmp` 期望=absent → 零命中 0 处 ✔
- **参考实现**：
  - `Nek5000/Nek5000` (T1, ok) — F11-parallel×4 core/byte_mpi.f
  - `Nek5000/nekRS` (T1, ok) — F11-parallel×217 3rd_party/adios/bindings/C/adios2/c/adios2_c_adios_mpi.cpp
  - `PyFR/PyFR` (T1, ok) — F11-parallel×44 pyfr/backends/base/generator/gpu.py
  - `CEED/libCEED` (T2, ok) — F11-parallel×122 backends/ceed-backend-list-cuda.h
  - `mfem/mfem` (T1, ok) — F11-parallel×8 general/cuda.cpp
  - `AMReX-Codes/amrex` (T2, ok) — F11-parallel×47 Src/Base/AMReX_CudaGraph.H
- **验收判据（机器可裁）**：2D Euler（p=3，512x512 单元）在 4 个 MPI 秩与 1 秩上：密度 L2 差 <=1e-10（同一离散解）；4 秩并行效率 >=0.7（相对单秩墙钟）；若实现 GPU 路径，同算例相对 CPU 单节点加速 >=3 倍且 L2 差 <=1e-10。
- **上游要点**：curated：pyfr 的 backend 抽象（同一 FR 算子多后端）、nekRS 的 matrix-free GPU 元素算子、libceed 的受限算子（restriction/prolongation）设计是三条可选路线。本波未读源码。
- **外部锚**：Fischer, Lottes, Kerkemeier et al., Nek5000 user documentation (2008)；Vincent, Witherden & Vermeire, Comput. Fluids 141 (2016)（PyFR）；Fischer et al., Parallel Comput. 109 (2022)（nekRS）。
- **影响**：三维实船/实海况 DG 算例的单进程上限在两三万单元量级；无并行则 DG 只能停留在二维验证层，无法进入工程算例。

### `GAP-DG-006` — 非线性隐式/IMEX 时间推进与 p-多重网格预处理

- **优先级**：P1；**工作量**：L；**功能点**：`F10-implicit`, `F09-time-integration`
- **现状（TensorDG）**：时间推进只有显式 SSP-RK3；隐式侧仅在 SIPG 椭圆求解里有块 Jacobi 预条件的 CG（sipg_matrixfree.py / poisson_2d.py），没有非线性 Newton-Krylov、没有 p-多重网格、没有 IMEX/ARK 分裂。
- **现状证据（grep @ 281ba5c6cf67）**：
  - `precondition` 期望=present → 命中 14 处 ✔
    - `tensordg/heat_1d.py:6: * ``"matrix-free"``  -- no global assembly, solve with preconditioned CG`
    - `tensordg/poisson_2d.py:11: * ``"matrix-free"``  -- no global assembly, solve with preconditioned CG`
    - `tensordg/sipg_matrixfree.py:1: """Matrix-free SIPG elliptic solvers using preconditioned CG.`
  - `newton_krylov|newton-krylov|gmres|multigrid|imex|arkode` 期望=absent → 零命中 0 处 ✔
- **参考实现**：
  - `exadg/exadg` (T1, ok) — F10-implicit×51 include/exadg/convection_diffusion/preconditioners/multigrid_preconditioner.cpp；F09-time-integration×61 include/exadg/acoustic_conservation_equations/time_integration/time_int_abm.h
  - `petsc/petsc` (T2, ok) — F10-implicit×31 src/ksp/ksp/impls/gmres/dgmres/dgmres.c；F09-time-integration×8 src/ts/impls/arkimex/arkimex.c
  - `LLNL/sundials` (T2, ok) — F10-implicit×8 bindings/sundials4py/sunnonlinsol/sunnonlinsol_newton.cpp；F09-time-integration×152 bindings/sundials4py/arkode/arkode.cpp
  - `dealii/dealii` (T1, ok) — F10-implicit×70 bundled/boost-1.84.0/include/boost/contract/core/constructor_precondition.hpp；F09-time-integration×65 bundled/boost-1.84.0/include/boost/numeric/odeint/stepper/adams_bashforth.hpp
  - `CEED/libCEED` (T2, ok) — F10-implicit×1 interface/ceed-preconditioning.c
- **验收判据（机器可裁）**：2D 对流扩散（Pe=1e3，p=3）推到稳态：Newton-Krylov + p-多重网格预条件的线性迭代次数从 64x64 到 256x256 增长 <=1.5 倍，达到同一收敛判据所需时间步数 <=显式 SSP-RK3 的 1/10；稳态解与显式长时积分解的 L2 差 <=1e-6。
- **上游要点**：curated：exadg 展示 matrix-free DG + 几何多重网格的组合；sundials/ARKODE 提供 IMEX 与嵌入式误差控制设计；petsc 提供现成的 SNES/PC 栈。本波未读源码。
- **外部锚**：Kronbichler & Kormann, ACM TOMS 45(3) 2019（matrix-free DG 算子求值）；Pazner, SIAM J. Sci. Comput. 42(4) 2020（p-multigrid for DG）；Hindmarsh et al., ACM TOMS 31 (2005)（SUNDIALS）。
- **影响**：粘性主导/低马赫船体绕流的显式时间步受 CFL 与刚性项双重限制；没有隐式-多重网格路线，DG 的工程可用性停在可压显式算例。

### `GAP-DG-007` — 目标导向后验误差估计（DWR）与自适应判据

- **优先级**：P2；**工作量**：M；**功能点**：`F12-error-estimation`
- **现状（TensorDG）**：有可微后端的伴随灵敏度（differentiable.py: adjoint sensitivity of a sensor）与若干算例内的经验误差估计，但没有对偶加权残差（DWR）估计子、没有估计-真误差对照、没有把估计子接到自适应决策上。
- **现状证据（grep @ 281ba5c6cf67）**：
  - `adjoint` 期望=present → 命中 23 处 ✔
    - `tensordg/differentiable.py:235: natural differentiable output functional (adjoint sensitivity of a sensor).`
    - `tensordg/kovasznay.py:350: wall_adjoint: bool = False,`
    - `tensordg/kovasznay.py:458: wall_penalty, wall_adjoint=wall_adjoint,`
  - `dual_weighted|dual-weighted|dwr|a_posteriori|a-posteriori|goal_oriented|goal-oriented` 期望=absent → 零命中 0 处 ✔
- **参考实现**：
  - `McGillCompAero/PHiLiP` (T1, ok) — F12-error-estimation×4 src/functional/adjoint.cpp
  - `firedrakeproject/firedrake` (T2, ok) — F12-error-estimation×24 firedrake/adjoint/__init__.py
  - `dealii/dealii` (T1, ok) — F12-error-estimation×5 include/deal.II/numerics/error_estimator.h
- **验收判据（机器可裁）**：2D Poisson(SIPG) 与 2D Euler 各 2 个算例上，DWR 自适应 5 轮后目标泛函（阻力/线积分）误差 <=同自由度均匀加密解的 50%，效应性指数（估计/真误差）落在 0.8~1.25 区间（4 个算例全部满足）。
- **上游要点**：curated：philip 的 DWR 估计流程、firedrake/pyadjoint 的伴随装配、dealii 的目标泛函装配模式是三条参考路线。本波未读源码。
- **外部锚**：Becker & Rannacher, Acta Numerica 10 (2001)（DWR 双加权残差）；Fidkowski & Darmofal, AIAA J. 49(4) 2011（output-based error estimation and mesh adaptation）。
- **影响**：把有限自由度花在工程真正关心的量（阻力、升力、桨盘推力）而不是全局 L2 上，是 DG 相较二阶 FVM 的卖点能否兑现的关键。

