# DG 开源技术盘点 registry（FEAT-OSS-INVENTORY-DG）

> 程序化生成，勿手改。源：`docs/oss_inventory_dg/registry.json` （`tools/oss_inventory_wave/wave.py build`，生成时间 2026-11-13T01:51:06+00:00）。
> 目标仓 `TensorDG` @ `281ba5c6cf6777ced440b8df756ad4402287d696`（证据来源：live_grep）。

## 1 计数

| 项 | 值 |
|---|---:|
| `pool_entries` | 23 |
| `upstream_repos` | 22 |
| `surveyed_ok` | 20 |
| `survey_failed` | 2 |
| `not_surveyed` | 0 |
| `negative_controls` | 1 |
| `negative_controls_recorded_failed` | 1 |
| `feature_points` | 14 |
| `gaps` | 7 |
| `gaps_evidence_satisfied` | 7 |

## 2 Tier 分级

| Tier | 数量 | 定义 |
|---|---:|---|
| `T1` T1 可吸收 | 13 | absorb now: method-defining DG/FR/spectral-element implementation whose technique maps onto a live TensorDG gap; reading the upstream source is enough to re-implement, licence is permissive or weak-copyleft |
| `T2` T2 参考设计 | 8 | reference design: framework/library that carries the enabling machinery (AMR forests, matrix-free kernels, Krylov/IMEX stacks, MPI abstractions); read for design and cross-check numbers, not scheduled for porting |
| `T3` T3 背景 | 1 | context only: adjacent paradigm (e.g. production FVM) or closed/hostile licence; kept for the maturity baseline and benchmark provenance, never for code reuse |

## 3 证据政策

- **`measured`**：metadata.* and feature_paths.* are fetched in this run from the GitHub REST API (repo metadata + one recursive tree at the default branch); source=github_api marks measured values
- **`curated`**：pool.json tier/focus/why and gaps.json refs/upstream_note are curation (curated_knowledge): they were NOT verified by reading the upstream sources in this run; adoption requires a source-reading pass in wave 2
- **`unresolved`**：entries with status=failed keep null metadata and must not be quoted as current project facts; the negative_control entry exists to prove failures are recorded, not dropped
- **`target_evidence`**：every GAP-DG current-state claim is a grep over the target repository at a pinned commit; patterns, hit counts and file:line samples are stored in registry.json (and re-run by `verify`)

## 4 仓位（普查结果）

| id | repo | tier | status | stars | licence(SPDX) | 可用性判定 | pushed_at | 源码文件数 |
|---|---|---|---|---:|---|---|---|---:|
| `trixi` | `trixi-framework/Trixi.jl` | T1 | ok | 737 | MIT | `port_ok` | 2026-11-12 | 331 |
| `exadg` | `exadg/exadg` | T1 | ok | 194 | GPL-3.0 | `read_only_reimplement` | 2026-11-12 | 662 |
| `philip` | `McGillCompAero/PHiLiP` | T1 | ok | 59 | NOASSERTION | `read_only_until_resolved` | 2026-09-08 | 384 |
| `flexi` | `flexi-framework/flexi` | T1 | ok | 299 | GPL-3.0 | `read_only_reimplement` | 2026-11-01 | 303 |
| `pyfr` | `PyFR/PyFR` | T1 | ok | 585 | NOASSERTION | `read_only_until_resolved` | 2026-11-12 | 253 |
| `nektar` | `nektar/nektar` | T1 | failed | — | — | `no_metadata` | — | — |
| `nek5000` | `Nek5000/Nek5000` | T1 | ok | 493 | NOASSERTION | `read_only_until_resolved` | 2026-10-04 | 378 |
| `nekr` | `Nek5000/nekRS` | T1 | ok | 431 | NOASSERTION | `read_only_until_resolved` | 2026-05-24 | 8284 |
| `dune_fem` | `dune-community/dune-fem` | T1 | failed | — | — | `no_metadata` | — | — |
| `ngsolve` | `NGSolve/ngsolve` | T1 | ok | 594 | LGPL-2.1 | `adapt_at_file_level` | 2026-11-12 | 537 |
| `mfem` | `mfem/mfem` | T1 | ok | 2291 | BSD-3-Clause | `port_ok` | 2026-11-13 | 819 |
| `dealii` | `dealii/dealii` | T1 | ok | 1758 | NOASSERTION | `read_only_until_resolved` | 2026-11-12 | 17539 |
| `gkeyll` | `gkeyllorg/gkeyll` | T1 | ok | 121 | MIT | `port_ok` | 2026-11-12 | 6407 |
| `p4est` | `cburstedde/p4est` | T2 | ok | 321 | GPL-2.0 | `read_only_reimplement` | 2026-11-11 | 109 |
| `amrex` | `AMReX-Codes/amrex` | T2 | ok | 832 | NOASSERTION | `read_only_until_resolved` | 2026-11-12 | 823 |
| `parthenon` | `parthenon-hpc-lab/parthenon` | T2 | ok | 175 | NOASSERTION | `read_only_until_resolved` | 2026-11-13 | 370 |
| `libceed` | `CEED/libCEED` | T2 | ok | 273 | BSD-2-Clause | `port_ok` | 2026-11-12 | 375 |
| `petsc` | `petsc/petsc` | T2 | ok | 550 | NOASSERTION | `read_only_until_resolved` | 2026-11-12 | 2711 |
| `sundials` | `LLNL/sundials` | T2 | ok | 706 | BSD-3-Clause | `port_ok` | 2026-11-13 | 692 |
| `dolfinx` | `FEniCS/dolfinx` | T2 | ok | 1249 | LGPL-3.0 | `adapt_at_file_level` | 2026-11-11 | 207 |
| `firedrake` | `firedrakeproject/firedrake` | T2 | ok | 682 | NOASSERTION | `read_only_until_resolved` | 2026-11-13 | 241 |
| `openfoam` | `OpenFOAM/OpenFOAM-dev` | T3 | ok | 2306 | NOASSERTION | `read_only_until_resolved` | 2026-11-12 | 9311 |
| `negctl` | `tp-agent/dg-negative-control-does-not-exist` | T3 | failed | — | — | `no_metadata` | — | — |

## 5 功能点 × 仓覆盖矩阵

单元格=该仓源码路径中命中此功能点签名路径数（tests/docs/vendor 已剔除）；`—`=该仓未取到树。

| repo | F01-basis | F02-flux | F03-limiting | F04-bound-preserving | F05-entropy-stable | F06-adaptivity | F07-simplex | F08-curved | F09-time-integration | F10-implicit | F11-parallel | F12-error-estimation | F13-multiphysics | F14-dispersion |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| `trixi` | 3 | 4 | 18 | 25 | 9 | 39 | 12 | 1 | 13 | 0 | 13 | 0 | 14 | 3 |
| `exadg` | 1 | 2 | 4 | 0 | 0 | 4 | 1 | 0 | 61 | 51 | 1 | 0 | 218 | 10 |
| `philip` | 5 | 6 | 12 | 10 | 4 | 46 | 4 | 11 | 50 | 3 | 0 | 4 | 16 | 6 |
| `flexi` | 5 | 18 | 2 | 0 | 0 | 0 | 0 | 4 | 0 | 4 | 2 | 0 | 35 | 1 |
| `pyfr` | 0 | 0 | 0 | 0 | 0 | 1 | 1 | 0 | 2 | 16 | 44 | 0 | 1 | 0 |
| `nektar` | — | — | — | — | — | — | — | — | — | — | — | — | — | — |
| `nek5000` | 0 | 0 | 0 | 0 | 0 | 2 | 5 | 0 | 0 | 1 | 4 | 0 | 9 | 0 |
| `nekr` | 2 | 1 | 0 | 0 | 0 | 17 | 0 | 0 | 0 | 87 | 217 | 2 | 27 | 1 |
| `dune_fem` | — | — | — | — | — | — | — | — | — | — | — | — | — | — |
| `ngsolve` | 1 | 0 | 0 | 0 | 0 | 1 | 1 | 0 | 0 | 13 | 54 | 0 | 5 | 5 |
| `mfem` | 6 | 1 | 10 | 0 | 0 | 2 | 5 | 21 | 0 | 6 | 8 | 2 | 31 | 0 |
| `dealii` | 141 | 6 | 182 | 18 | 2 | 343 | 11 | 0 | 65 | 70 | 827 | 5 | 7 | 4 |
| `gkeyll` | 300 | 1146 | 543 | 14 | 0 | 53 | 38 | 0 | 5 | 5 | 21 | 0 | 1367 | 7 |
| `p4est` | 0 | 0 | 0 | 0 | 0 | 49 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| `amrex` | 1 | 26 | 1 | 0 | 0 | 766 | 0 | 0 | 0 | 6 | 47 | 0 | 0 | 0 |
| `parthenon` | 0 | 1 | 0 | 0 | 0 | 12 | 0 | 0 | 4 | 2 | 10 | 0 | 0 | 3 |
| `libceed` | 56 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 1 | 122 | 0 | 0 | 0 |
| `petsc` | 3 | 0 | 0 | 0 | 0 | 47 | 6 | 1 | 8 | 31 | 209 | 2 | 2 | 1 |
| `sundials` | 0 | 0 | 0 | 0 | 0 | 31 | 0 | 0 | 152 | 8 | 47 | 20 | 0 | 0 |
| `dolfinx` | 0 | 0 | 0 | 0 | 0 | 13 | 1 | 0 | 0 | 0 | 4 | 0 | 0 | 0 |
| `firedrake` | 0 | 0 | 3 | 0 | 0 | 1 | 0 | 0 | 0 | 19 | 14 | 24 | 0 | 1 |
| `openfoam` | 46 | 116 | 546 | 2 | 0 | 141 | 21 | 0 | 0 | 49 | 46 | 11 | 101 | 45 |

功能点定义：

- **`F01-basis`** 模态/节点基与数值求积 — 路径签名 `(basis|legendre|quadrature|vandermonde|polynomial|orthonormal)`
- **`F02-flux`** 数值通量（Riemann 求解器） — 路径签名 `(riemann|hllc|rusanov|flux|upwind|godunov)`
- **`F03-limiting`** 限制器 / 激波捕捉 — 路径签名 `(limit|weno|tvd|minmod|tvb|moment|artificial_viscosity)`
- **`F04-bound-preserving`** 保界/正性保持（子单元限界） — 路径签名 `(positiv|bound_preserv|bound-preserv|subcell|zhang|convex|mood)`
- **`F05-entropy-stable`** 熵稳定/分裂型离散 — 路径签名 `(entropy|split_form|split-form|kuscu|esfr|two-point)`
- **`F06-adaptivity`** h/p 自适应与悬挂节点 — 路径签名 `(adapt|amr|p4est|hanging|refine|coarsen|octree)`
- **`F07-simplex`** 单纯形（三角/四面体）单元与网格 — 路径签名 `(triangle|tetra|simplex|gmsh|unstructured)`
- **`F08-curved`** 曲线单元/等参几何 — 路径签名 `(curved|isoparametric|jacobian|nurbs|metric_term)`
- **`F09-time-integration`** 时间推进（SSP-RK / IMEX） — 路径签名 `(time_integ|ssprk|ssp_rk|rk3|imex|arkode|stepper|ode_solver)`
- **`F10-implicit`** 隐式求解（Newton-Krylov/多重网格/预处理） — 路径签名 `(newton|gmres|krylov|multigrid|precondition|implicit)`
- **`F11-parallel`** 并行与加速器后端（MPI/GPU） — 路径签名 `(mpi|kokkos|cuda|hip|gpu|openmp|openacc|device)`
- **`F12-error-estimation`** 后验误差估计/目标导向自适应 — 路径签名 `(adjoint|dual_weighted|dwr|a_posteriori|error_estimat)`
- **`F13-multiphysics`** 多物理（NS/Maxwell/弹性/浅水/等离子体） — 路径签名 `(navier|maxwell|elastic|shallow|swe_|plasma|gyrokinetic|mhd|turbulen)`
- **`F14-dispersion`** 色散/耗散分析与谱性质 — 路径签名 `(dispersion|dissipation|eigen|spectrum|spectral_analysis)`

## 6 取数失败（如实记录，不静默丢弃）

- `nektar` (`nektar/nektar`) — `not_on_github (upstream hosted at https://gitlab.com/nektar/nektar); HTTP 404`
- `dune_fem` (`dune-community/dune-fem`) — `not_on_github (upstream hosted at https://gitlab.dune-project.org/dune-fem/dune-fem); HTTP 404`
- `negctl` (`tp-agent/dg-negative-control-does-not-exist`) — `HTTP 404`
