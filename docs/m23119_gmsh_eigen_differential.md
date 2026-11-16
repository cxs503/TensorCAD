# Gmsh 对拍：Eigen 开关打开后，ICN 家族从 nan/0 变成数字

*锚点:* `GAP-CAD-D1-MESH-QUALITY-LIB`（网格质量指标库 aspect/skew/nonorth/jacobian/水密/流形/取向）
*验收:* 与 meshmagick+Gmsh 指标对拍（≥3 网格、指标差异 < 1e-6 **或差异解释**）；可作 CI 门用
*本单:* `m23119`（platform 席自主工时） · *工件:* `artifacts/mesh_quality_external_diff/gmsh_source_diff.json`
*复现:* `sh tools/mesh_quality_oracles/build_gmsh_source_oracle.sh <gmsh-4.15.2-src> /tmp/gso/bin` 然后
`python tools/mesh_quality_gmsh_source_diff.py --oracle /tmp/gso/bin --require-oracle`

## 1. 起点：Gmsh 一侧「有代码、没数字」

`docs/mesh_quality_external_oracles.md` §4/§6 与 `docs/m22768_gmsh_source_oracle_build.md` §3 记录：
Gmsh 源构建 oracle 能链接、能跑、tet `volume` 对得上，但 `gamma`/`eta` 读回 **0.0**，
`minSICN` 读回 **nan / 0 / DBL_MAX**，并把原因写成一个「一行 dispatch 修复」的假设
（`MElement.h:257` 基类虚函数）。

本单先复核这个假设：重新构建后 `gamma`/`eta` **直接就是对的**（不是 0）——m22768 的
dispatch 假设不成立，真正的原因在下面这一条，而且是**静默**的：

```text
$ /tmp/gso_plain/bin hex_unit.mkm     # HAVE_MESH=1, HAVE_EIGEN undefined
W hex minSICN 0 0                       # ...
$ /tmp/gso_plain/bin tri_warped.mkm
Matrix inversion requires Eigen or LAPACK
W tri minSICN 0 nan
```

Gmsh 的 `minSICN`/`minSIGE` 走 `CondNumBasis` → `src/numeric/fullMatrix.h:670`，逆矩阵是
`#if defined(HAVE_EIGEN)`；没有 Eigen 就打印一句、返回 nan/0，而二进制**退出码 0**。
对拍器会把 nan 当成"数据"，或把 0 当成"Gmsh 说这个单元极差"。这比没有 oracle 更危险。

## 2. 做了两件事

1. **构建脚本（`tools/mesh_quality_oracles/build_gmsh_source_oracle.sh`）**
   - 除 `HAVE_MESH` 外自动探测 Eigen 头文件（`TP_EIGEN_INCLUDE` → 源树 `contrib/eigen` →
     `/usr/include/eigen3`），找到就加 `HAVE_EIGEN`，并在构建行**打印开关状态**
     （`HAVE_MESH=1, HAVE_EIGEN=1 via /usr/include/eigen3`）——nan 型 oracle 不能再无声通过；
   - 缺 `qualityMeasuresJacobian.cpp` 的树直接 `exit 2`（否则静默少一个翻译单元、少一族指标）；
   - 目标文件改为**并行编译**（`xargs -P $TP_JOBS`, 默认 6）：串行 ~70 s（无 Eigen）、>2 min（含 Eigen），
     含 Eigen 后单条命令窗口内跑不完就会被丢掉。实测 **53 s**（6 路）出二进制。
2. **对拍驱动器（`tools/mesh_quality_gmsh_source_diff.py`）** 把判据拆成两张互斥的表：
   - `MAPPING`＝**必须一致**（差异即 FAIL，退出码 1）；
   - `EXPLAINED`＝**测出来的、按定义就不该相等**的配对：记录实测差异 + 解释，不参与退出码。
     两张表互斥由测试强制（一个配对只能有一个身份），解释缺句也由测试拒绝。

## 3. 计算结果（程序化产出，非手抄）

`pairs compared 12, beyond 1e-06: 0; explained 8`（退出码 0）。6 网格 / 4 单元类型。

### 3.1 必须一致（MAPPING）：12 对，0 对超差

| 网格 | 类型 | 指标（库 ↔ Gmsh） | 最大相对差 | 超差元素 |
|---|---|---|---:|---:|
| `tri_warped` | tri | `gamma` / `eta` / `min_sicn`↔`minSICN` | 5.16e-16 / 5.95e-16 / 8.06e-16 | 0/40 |
| `quad_skewed` | quad | `eta` / `area` / `min_sicn`↔`minSICN` | 1.19e-16 / 3.33e-16 / 0 | 0/12 |
| `tet_cube6` | tet | `gamma` / `eta` / `volume` | 0 / 0 / 0 | 0/6 |
| `tet_sliver` | tet | `gamma` / `eta` / `volume` | 1.85e-16 / 2.84e-16 / 1.37e-16 | 0/4 |
| `hex_unit`, `hex_distorted` | hex | （见 3.2：无"必须一致"配对） | — | — |

三角形/四边形 ICN 与 Gmsh `minSICN` **逐元素相同到 1e-16**，是本锚点第一次拿到
Gmsh 的 Jacobian 家族数字（此前一直记账为 NOT RUN）。

### 3.2 差异解释（EXPLAINED）：8 对

| 网格 | 库指标 ↔ Gmsh | 最大相对差 | 为什么不该相等 |
|---|---|---:|---|
| `tet_cube6` | `gmsh_icn` ↔ `minSICN` | **3.55e-01** | 单位相反：实测 Gmsh 对**正四面体**报 **1.0**、对 Gmsh 自己的直角参考四面体报 **0.8165**；库的 `gmsh_icn` 对同样两体报 0.7937 / 1.0。且比值非常数（sliver 上 0.92）——不是换算，是另一个量 |
| `tet_sliver` | `gmsh_icn` ↔ `minSICN` | **9.24e-01** | 同上；`metrics.py` 里"同一量、同一单位、比较即恒等式"的说法被证伪 |
| `tet_cube6` / `tet_sliver` | `min_sicn` ↔ `minSICN` | 2.41e-02 / 9.27e-01 | Gmsh 用**理想（参考）Jacobian + 理想梯度**、并对单元做 Bezier 细分求下界；库取物理角点 Jacobian。正四面体上两者都 =1，离理想后不同 |
| `hex_distorted` | `min_sicn` ↔ `minSICN` | 8.99e-02 | 同上（角点最小 vs 单元内下界），两向偏离 |
| `hex_distorted` | `volume` ↔ `volume` | **1.34e-01** | 两个"体积"：Gmsh 对三线性映射积分 `det J`（`MElement::getVolume` + Gauss）；库用散度定理把 6 面劈成 12 三角形（`_hex_volume`）。平面面片上完全一致（1.33e-15），面片一翘就分叉 |
| `hex_unit` | `volume` ↔ `volume` | 1.33e-15 | 同上，理想立方体 |
| `hex_unit` | `min_sicn` ↔ `minSICN` | 2.22e-16 | 同上，理想立方体两者都 =1 |

Gmsh **自己回答不了**的：hex 的 `gamma`/`eta`（`MHexahedron` 不 override，恒 0.0）；
tri/quad/tet 的 `minSIGE`（nan，只有 hex 有值）。已在 `NOT_COMPARABLE` 里按实测值登记。

## 4. 误差判定（逐项，显式回答）

| 量 | 基准 | 实测 | 是否 <1e-6 | 是否 <3% |
|---|---|---|---:|---:|
| tri `gamma`,`eta`,`min_sicn`（40 元素） | Gmsh 4.15.2（源构建，执行） | ≤8.06e-16 | **是** | 是 |
| quad `eta`,`area`,`min_sicn`（12 元素） | 同上 | ≤3.33e-16 | **是** | 是 |
| tet `gamma`,`eta`,`volume`（10 元素） | 同上 | ≤2.84e-16 | **是** | 是 |
| hex `volume`（平坦面片，8 元素） | 同上 | 1.33e-15 | **是** | 是 |
| hex `volume`（翘曲面片，12 元素） | 同上 | 1.34e-01 | 否 | **否（13.4%）** — 定义不同，见 3.2 |
| tet `min_sicn`/`gmsh_icn` | Gmsh `minSICN` | 2.4e-2 … 9.2e-1 | 否 | **否** — 定义不同，见 3.2 |

结论：验收的「≥3 网格、差异<1e-6 **或差异解释**」在 Gmsh 一侧**成立**（12 对硬指标全 <1e-6，
8 对按定义分离且给出数字与理由）。原先记账的 NOT RUN 现在是 RUN。

## 5. 注意点（踩过的坑）

- **静默全零/nan oracle 比没有 oracle 更坏**：`HAVE_MESH`、`HAVE_EIGEN` 任一缺失，二进制都
  `exit 0`。构建行现在把两个开关的状态打印出来，`tests/test_m23119_gmsh_eigen_differential.py`
  同时断言开关逻辑存在。
- **`libgmsh` 路线不可用**：Debian `libgmsh.so.4.15.2` 需要 ~30 个 OpenCASCADE/FLTK/METIS 依赖，
  本镜像没有；所以走源构建（Gmsh 4.15.2 源码树 + 自造 `GmshConfig.h`）。
- **hex 家族没有 Gmsh 的 gamma/eta**，不是机器问题，是 Gmsh 定义如此（基类返回 0）。
- **不要把 `EXPLAINED` 当通过**：它只记录，不判分；若哪天差异变成 0，测试会因"解释可能是错的"而报警。
- 上一张单的结论（m22768：gamma/eta 归零是 dispatch 问题）**已被实测否定**，本文取代该条。

## 6. 局限与未做完的部分

- Eigen 版本/来源只记录到路径（`/usr/include/eigen3`），未锁版本号；
- `min_sicn`（tet/hex）与 Gmsh `minSICN` 的**公式对齐**未做——本单只把差异量化并给出定性解释；
  若要真正同单位，需要按 Gmsh 的 `getSignedIdealJacobian` + Bezier 下界实现一个"gmsh 兼容"指标，
  或在文档中把库的 `min_sicn` 明确标注为"角点采样版"。
- `metrics.py` 中 `gmsh_icn` 的单位说明**尚未改正**（行为未动，避免影响既有基线）；建议下单一并修文档/公式。
- 本单不改任何指标数值，CI 门仍是既有 `mesh-quality-gate.yml` + `door.admit`；新增的是**外部对拍门**
  （`--require-oracle` 退出码 2 = oracle 未跑，绝不假装通过）。