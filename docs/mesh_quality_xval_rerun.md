# 网格质量指标 × Verdict 对拍：把「已声明」变成「已运行」

锚点 `GAP-CAD-D1-MESH-QUALITY-LIB`（网格质量指标库：aspect/skew/nonorth/jacobian/水密/流形/取向）
验收要求「与 meshmagick + Gmsh 指标对拍（≥3 网格、指标差异 <1e-6 或差异解释）」。
本单（m22268）之前的真实状态是**三处欠账**：

1. 记证工件 `artifacts/mesh_quality_external_diff/oracle_diff.json` **在仓库里不存在**，
   `tests/test_mesh_quality_external_diff.py` 里断言「已录制证据」的 6 个用例长期被
   `pytest.skip` 静默跳过 —— 也就是对拍结论只有散文，没有可机裁的证据；
2. 协议层 2 个用例是**红**的（`compare()` 语义用例写错夹具；oracle 返回长度不足时闸门不失败）；
3. 文档声称「只有 4 处已解释差异」，但那是**在 oracle 未运行**的前提下写的；
   本机 `verdict_oracle`（预编译 loongarch 二进制，`work/_ext/bin/`）其实**可以跑**。

## 1. 做了什么

**(a) 让 Verdict oracle 真正跑起来并落盘证据。**
`work/_ext/bin/verdict_oracle` 在本机可直接执行（`ldd` 无缺失依赖），语料 6 个网格
（tri/quad/tet/hex 四类、含退化元）逐元输出。工件由
`python tools/mesh_quality_external_diff.py --verdict-oracle work/_ext/bin/verdict_oracle`
程序化生成，零手抄。

**(b) 修正了三处把「不同定义」当成「同一指标」的映射。**
上一版 `MAP["verdict"]` 把下列三对当同量比较，因而必然常年报 DIFFER；本单实测
`external/ours` 比值**并非常数**（故不是单位换算，而是定义不同），改为在 `EXPLAINED`
中显式声明并把实测比值区间写进工件：

| 网格/指标 | external/ours 实测区间 | 判定 | 依据 |
|---|---|---|---|
| `tri_warped / aspect_ratio` | 1.1365 … 1.1547 (=2/√3) | 不同量，非单位换算 | Verdict 三角 `aspert_ratio` 的归一化与我们不同 |
| `tet_sliver / skewness` | 1.0 … 1.5704 | 不同量 | Verdict `tet_equiangle_skew`(V_TetMetric.cpp) 用**二面角**、平衡角 70.5287793655093°；我们用**面角**、平衡角 60° |
| `tet_sliver / min_dihedral` | 1.0 … 12.014 | 不同量 | Verdict `tet_minimum_angle` 用**无符号**面法向点积（同一棱取内二面角的补角）；我们取内二面角本身 |

**(c) 把一处真单位差升级为「带校验的换算」。**
四面体 `scaled_jacobian`：Verdict 以**正四面体**为单位（正四面体=1），我们按 Jacobian 列范数
归一（正四面体=1/√2）。实测 `external/ours` 在**每个**单元上都等于 √2（min=max=1.414213562373，
判据 `|max-min| ≤ 1e-12`）。工具现在把 oracle 值除以 √2 再比对（结果 rel ≤1e-6 AGREE），
并**记录实测比值区间**：若哪天这个因子不再恒定，工具会直接把它判失败，而不是被常数掩盖。

**(d) 修闸门语义（CI 可用的前提）。**
oracle 返回的值个数与单元数不一致（`LENGTH_MISMATCH`）此前**不算失败**——即 oracle
被截断/少返回时闸门仍然绿。这是静默丢证据，现改为与 `DIFFER` 同等失败；
`compare()` 那个错夹具的用例也已改正（它断言 `LENGTH_MISMATCH` 的那一对根本没出现在 `ours` 里）。

## 2. 结果（工件 `artifacts/mesh_quality_external_diff/oracle_diff.json`，tolerance 1e-6）

| 网格 | 指标 | max rel diff | 判定 |
|---|---|---|---|
| hex_unit | scaled_jacobian / skewness | 0 / 0 | AGREE |
| hex_distorted | scaled_jacobian / skewness | 3.64e-16 / 6.16e-16 | AGREE |
| tet_cube6 | volume / edge_ratio | 0 / 0 | AGREE |
| tet_cube6 | scaled_jacobian（÷√2 后） | ≤1e-6 | AGREE（此前 2.93e-01） |
| tet_sliver | volume / edge_ratio | 1.37e-16 / 2.67e-16 | AGREE |
| tri_warped | area / edge_ratio / min_angle / max_angle / scaled_jacobian / skewness | 1.34e-16 / 1.56e-16 / 4.23e-16 / 1.70e-16 / 4.85e-16 / 6.30e-16 | AGREE |
| tri_warped | gamma vs 1/radius_ratio | 7.13e-16 | AGREE |
| quad_skewed | edge_ratio / min_angle / max_angle / scaled_jacobian / skewness | 1.66e-16 / 1.69e-16 / 0 / 2.23e-16 / 2.41e-15 | AGREE |

工件 `mismatches: []`（此前同一语料同一 oracle 下为 6 条 DIFFER，其中 3 条是定义混淆、1 条是未做单位换算）。

## 3. 注意点 / 局限（不得粉饰）

- **Gmsh 半边仍是 NOT_BUILT。** 本机只有 `gmsh` 的 **CLI 二进制**（链接 `libgmsh.so.4.15`），
  该 `.so` 的 `DT_NEEDED` 含整套 OpenCASCADE(24 个 `libTK*`)、FLTK、METIS、CGNS、voro++ 等
  本机不存在的库，且 PyPI 无 loongarch64 的 `gmsh` 轮子；`dpkg` 也不在本机。
  本单因此改为**从源码构建**（`gmsh_4.15.2+ds1.orig.tar.xz` 已在盘上，cmake 配置成功：
  `-DENABLE_FLTK=OFF -DENABLE_OCC=OFF -DENABLE_MED=OFF ...`，零外部依赖），构建已按计算调度令
  递交集群（作业 `cj4009`，32 核），交接时仍在排队。**Gmsh 的数值对拍仍未得到结论，
  本文件不对 Gmsh 作任何数值声明**；工件里 `gmsh_oracle.status = NOT_BUILT` 如实保留。
- Verdict 二进制是预编译件，其版本号未内嵌可读字符串，工件记 `Verdict(master)`（取自
  `work/_ext/verdict-master` 源码树），非严格 commit 标识。
- 语料是 6 个小网格（6–40 单元），覆盖理想与退化单元，但**不含**船体实尺度网格；
  `min_dihedral` 的补角差异在 90° 二面角处自然消失，故对接近正交的网格该差异不可见。
- 我们**没有**改动 `skewness` / `min_dihedral` / 三角 `aspect_ratio` 的定义去迎合 Verdict：
  改定义会移动既有阈值的含义（`tensorcad.gate` 的限值迁移依赖既有方言），属于另一单。
  建议后续条目：新增 `skewness_dihedral`（Verdict 定义，θe=70.5287793655093°）与
  `aspect_ratio_radius_ratio`（Verdict 定义），并给出与现有指标的代数关系，供 `MAP` 直接对拍。

## 4. 复现

```
python tools/mesh_quality_external_diff.py \
  --verdict-oracle work/_ext/bin/verdict_oracle \
  --artifact artifacts/mesh_quality_external_diff/oracle_diff.json
python -m pytest tests/test_mesh_quality_external_diff.py -q    # 17 passed（含 6 个证据用例，不再 skip）
```

CI 模式：`--require-oracles` 在任一 oracle 未运行时 exit 2（不允许「没跑=通过」）。

## 5. 误差判定（逐项 vs 基准，显式回答）

| 关键量 | 基准 | 相对误差 | <3%？ |
|---|---|---|---|
| hex scaled_jacobian / skewness | Verdict（运行） | ≤6.2e-16 | 是 |
| tri/quad 面积、边长比、角度、SJ、skew | Verdict（运行） | ≤2.4e-15 | 是 |
| tet volume / edge_ratio | Verdict（运行） | ≤2.7e-16 | 是 |
| tet scaled_jacobian（÷√2 换算后） | Verdict（运行） | ≤1e-6（换算因子恒定到 1e-15） | 是 |
| tet skewness / min_dihedral、tri aspect_ratio | Verdict（运行） | **不适用**：定义不同，比值 1.0–12.0 非常数 | 否（已声明为不同量，不计入误差） |
| Gmsh 各指标 | Gmsh | **无基准**（本机不可运行，构建排队中） | 不可判定 |