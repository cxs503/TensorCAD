# TensorCAD meshio 六格式矩阵往返保真 · 18 组合全绿表（m24628 复核）

**锚点**: GAP-CAD-B2-MESH-FORMAT-MATRIX  
**验收**: 往返保真——体积/拓扑/取向不变（6 格式×3 几何=18 组合全绿表）  
**本单**: m24628（generalist 自找活）  
**环境**: python 3.11.6 / meshio 5.3.5 / numpy 2.4.6  
**判定**: 24 组合（6 格式×4 几何，含验收 18 组合超集） all_pass=**True**

## 1. 全绿表（体积/拓扑/取向三不变量）

| 格式 | 几何 | 拓扑 | 体积(原始) | 体积(往返) | 相对误差 | 体积OK | 拓扑OK | 取向OK | 文件sha256(前12) |
|---|---|---|---|---|---|---|---|---|---|
| stl | cube | triangle | 1.000000000 | 1.000000000 | 0.00e+00 | True | True | True | cde4240d6bce |
| stl | sphere | triangle | 4.158971436 | 4.158971436 | 0.00e+00 | True | True | True | 1e20b0893345 |
| stl | capsule | triangle | 3.639460645 | 3.639460645 | 0.00e+00 | True | True | True | 749eb56d79cf |
| stl | wigley_hull | triangle | 2.771222150 | 2.771222150 | 0.00e+00 | True | True | True | ca9db2276ce2 |
| obj | cube | triangle | 1.000000000 | 1.000000000 | 0.00e+00 | True | True | True | b76fd62395f6 |
| obj | sphere | triangle | 4.158971436 | 4.158971436 | 0.00e+00 | True | True | True | 2231dd51728e |
| obj | capsule | triangle | 3.639460645 | 3.639460645 | 0.00e+00 | True | True | True | 9610cc59145f |
| obj | wigley_hull | triangle | 2.771222150 | 2.771222150 | 0.00e+00 | True | True | True | 730eab23dfac |
| vtk | cube | triangle | 1.000000000 | 1.000000000 | 0.00e+00 | True | True | True | eb791b46fb9e |
| vtk | sphere | triangle | 4.158971436 | 4.158971436 | 0.00e+00 | True | True | True | f7f9e8aa58e8 |
| vtk | capsule | triangle | 3.639460645 | 3.639460645 | 0.00e+00 | True | True | True | 5ef3cde820d2 |
| vtk | wigley_hull | triangle | 2.771222150 | 2.771222150 | 0.00e+00 | True | True | True | faea3dc259e7 |
| gmsh | cube | triangle | 1.000000000 | 1.000000000 | 0.00e+00 | True | True | True | 361975879656 |
| gmsh | sphere | triangle | 4.158971436 | 4.158971436 | 0.00e+00 | True | True | True | d875829b6316 |
| gmsh | capsule | triangle | 3.639460645 | 3.639460645 | 0.00e+00 | True | True | True | b3f636a3de61 |
| gmsh | wigley_hull | triangle | 2.771222150 | 2.771222150 | 0.00e+00 | True | True | True | 1d54ee2c1709 |
| su2 | cube | tetra | 1.000000000 | 1.000000000 | 0.00e+00 | True | True | True | 6d846a494090 |
| su2 | sphere | tetra | 4.158971436 | 4.158971436 | 0.00e+00 | True | True | True | ab8bb57e242a |
| su2 | capsule | tetra | 3.639460645 | 3.639460645 | 0.00e+00 | True | True | True | 367bbc3f49ba |
| su2 | wigley_hull | tetra | 2.771222150 | 2.771222150 | 0.00e+00 | True | True | True | f514c93bbdf6 |
| cgns | cube | tetra | 1.000000000 | 1.000000000 | 0.00e+00 | True | True | True | 62a400ee54ed |
| cgns | sphere | tetra | 4.158971436 | 4.158971436 | 0.00e+00 | True | True | True | 692590b538d8 |
| cgns | capsule | tetra | 3.639460645 | 3.639460645 | 0.00e+00 | True | True | True | 9d78afac5b28 |
| cgns | wigley_hull | tetra | 2.771222150 | 2.771222150 | 0.00e+00 | True | True | True | 1e2e8534d959 |

**判定**：24/24 组合三不变量全绿（all_combinations_pass=True）。验收要求 6×3=18 组合，本表为 6×4=24 组合超集（第 4 几何=真实 Wigley 船体）。

## 2. 几何门（离散体积 vs 解析基准）

| 几何 | 水密 | 体积(测量) | 体积(解析) | 相对误差 | 四面体数 |
|---|---|---|---|---|---|
| cube | True | 1.000000 | 1.000000 | 0.00% | 12 |
| sphere | True | 4.158971 | 4.188790 | 0.71% | 2208 |
| capsule | True | 3.639461 | 3.665191 | 0.70% | 1536 |
| wigley_hull | True | 2.771222 | 2.777778 | 0.24% | 1996 |

解析基准：cube=1.0（精确）；sphere=4π/3；capsule=πr²L+4πr³/3（r=0.5, L=4）；wigley_hull=4LBT/9（n=m=2 闭式，L=10, B=1, T=0.625）。离散体积误差全部 <1%（网格离散化误差，非格式误差）。

## 3. 证据指针

- 机器可读全绿表：`runs/m24628_meshio_matrix18/matrix18.json`
- 24 个往返文件（可第三方复核）：`runs/m24628_meshio_matrix18/files/`
- 图：`artifacts/m24628_meshio_matrix18/`（status_matrix / volume_relative_error / geometry_gate）
- 复现：`python scripts/meshio_matrix18.py --out runs/m24628_meshio_matrix18 --fig artifacts/m24628_meshio_matrix18 && python -m pytest tests/test_meshio_format_matrix.py -q`

