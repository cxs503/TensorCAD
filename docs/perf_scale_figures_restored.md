# TensorCAD 大几何性能基准：缺失图件恢复与回归验证（GAP-CAD-D4-PERF-SCALE）

**锚点**: `GAP-CAD-D4-PERF-SCALE`（大几何性能基准 10⁶–10⁷ 三角级）
**版本**: 2.0.1（图件恢复）· **实测节点**: LoongArch64, 128 核, Python 3.11.6, TensorCAD 0.23.0
**作业**: 本单（m23275）——图件由 `scripts/gen_perf_scale_figures.py` + `scripts/gen_perf_reference_figure.py`
从**受版本控制的工件**（`benchmarks/perf_scale/artifacts/*.json`）重新生成，非手抄。

---

## 1. 缺陷（验证欠账）

`docs/perf_scale_10m.md`（§3.3）与 `docs/perf_scale_reference.md`（§1）引用
`benchmarks/perf_scale/figures/*.png`（12 张 + 1 张参考图），但该目录在干净检出中**不存在**——
图件在合并门禁中被滤除（`artifacts/` 与 `figures/` 均被门禁过滤），而文档与测试仍引用它们。

**直接后果**（`tests/test_perf_benchmark.py` 与 `tests/test_perf_scale_evidence.py` 各 1 个失败）：
- `test_perf_scale_figures_exist`：`AssertionError: missing figures directory: .../benchmarks/perf_scale/figures`
- `test_perf_scale_figures_exist_and_are_valid_png`：同上

这是**验证欠账**：已合并功能（文档 + 工件）引用的图件未随分支落盘，图件有效性测试无法通过。

## 2. 修复

图件生成脚本（`scripts/gen_perf_scale_figures.py`、`scripts/gen_perf_3pt_figures.py`、
`scripts/gen_perf_import_figures.py`、`scripts/gen_perf_single10m_figures.py`、
`scripts/gen_perf_reference_figure.py`）全部从**受版本控制的工件 JSON** 读取数据（零手抄），
因此重新生成即恢复。执行：

```bash
python scripts/gen_perf_scale_figures.py      # 12 张（scaling/point_locate/index_cost/voxelize/fig3pt/import/single10m）
python scripts/gen_perf_reference_figure.py   # 1 张（reference_scaling_bands.png）
```

生成 13 张 PNG 到 `benchmarks/perf_scale/figures/`（受版本控制路径，随分支落盘）。

## 3. 回归验证（验证欠账）

对已合并功能做基线复核（`python -m tensorcad.perf_benchmark --quick --baseline benchmarks/perf_scale/baseline.json`）：

| 回归门 | 实测 vs 基线 | 比值 | 限值 | 判定 |
|---|---:|---:|---:|---|
| regression_load:21x9 | 0.0170s vs 0.0178s | 0.96× | 2.0× | ✅ |
| regression_load:41x13 | 0.0544s vs 0.0521s | 1.04× | 2.0× | ✅ |
| regression_load:81x21 | 0.1724s vs 0.1805s | 0.96× | 2.0× | ✅ |
| regression_load:161x41 | 0.7819s vs 0.7658s | 1.02× | 2.0× | ✅ |
| regression_volume:21x9 | 0.0003s vs 0.0003s | 0.95× | 2.0× | ✅ |
| regression_volume:41x13 | 0.0008s vs 0.0009s | 0.88× | 2.0× | ✅ |
| regression_volume:81x21 | 0.0026s vs 0.0030s | 0.85× | 2.0× | ✅ |
| regression_volume:161x41 | 0.0112s vs 0.0115s | 0.97× | 2.0× | ✅ |
| regression_voxelize:32x0 | 0.1806s vs 0.1812s | 1.00× | 2.0× | ✅ |
| scaling_load | 指数 1.04 | — | 1.50 | ✅ |
| scaling_point_locate | 指数 1.02 | — | 1.50 | ✅ |
| scaling_point_locate_accel | 指数 0.33 | — | 0.80 | ✅ |
| scaling_volume | 指数 1.03 | — | 1.50 | ✅ |

**结论**：全部墙钟回归门在 2× 限值内（比值 0.85–1.04），复杂度指数全部低于限值——
**已合并功能相对基线无回归**。`point_locate` 的 `not_run` 门是**诚实标记**（本单用 200 条射线，
基线记录 1000 条，查询时间跨查询数不可比），不是静默通过。

## 4. 测试

```bash
python -m pytest tests/test_perf_benchmark.py tests/test_perf_scale_artifact.py \
  tests/test_perf_scale_evidence.py tests/test_perf_scale_reference_doc.py -q
# 110 passed（含此前 2 个失败的图件有效性测试）
```

## 5. 误差判定

| 量 | 基准 | 实测 | 相对误差 | <3%? |
|---|---|---|---|---|
| 体积 @25,918 tri | 7.512698 m³（解析） | 7.502488（as-imported） | −0.136% | 是 |
| 体积 @4,003,998 tri | 7.512698 m³ | 7.512474 | −0.003% | 是 |
| 体积 @12,011,994 tri（3 体 fleet） | 22.538095 m³ | 22.537421 | −0.003% | 是 |
| 回归门比值（全部） | 2.0× 限值 | 0.85–1.04× | 全部低于限值 | 是 |

**结论**：关键几何量全部 <3%；图件恢复不改变任何测量（图件从同一工件重新渲染），
回归验证确认已合并功能无退化。
