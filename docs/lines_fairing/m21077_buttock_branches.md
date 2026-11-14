# 型值表三向光顺：第三向（buttock）等值线分支修复 — m21077

锚点：`GAP-CAD-C3-LINES-FAIRING`（型值表三向光顺／曲率能量最小化）。
对象：本仓 `tensorcad/lines_fairing.py`（main 上已有实现，本单不改锚点、不做重写，只做实质缺陷修复 + 覆盖度验证）。

## 1 发现的缺陷（两条，均可复现）

第三向（buttocks，`y = const` 的纵剖线）不是型值表里直接给出的，而是从 `eta(x, z)` 网格**等值线**提取的。
原实现 `buttock_curves` 有两个结构性缺陷：

**D1 把剖面当单调处理**：每个站位的判定只有 `if level < profile[0] or level > profile[-1]: continue`，
随后在 `zeta` 方向二分。一旦某剖面的半宽在最大宽度之上**回落**（tumblehome／舷弧内倾，实船艉体常见），
该 level 的终点值低于 level，于是整条 buttock 被丢弃 —— 而它其实被该剖面**穿过两次**（双分支）。
可复现：`docs/lines_fairing/buttock_branch_evidence.py` 里的 tumblehome 面（`z/T>0.65` 内倾 55%），
level = 0.5·peak 时**旧算法返回 0 条曲线**，新算法返回 **2 条**（`zeta = 0.3588` 与 `0.8668`）。

**D2 level 取自全局峰值的固定比例**（0.2/0.35/0.5/0.65/0.8）：靠近甲板/龙骨的比例在某些船型上只有两三个站位
能越过，而 `buttock_curves` 又静默丢弃少于 3 点的曲线；分族指标因此“看起来跑了、其实没测”。
现改为 `buttock_levels`：在峰值 5%–95% 间扫 61 个候选，保留**最长追踪分支样本数 ≥ max(5, 0.3·N_stations)**
者，再等距取 `count` 条。

同时修掉两处静默失效：
- `_plot` 原先把两条曲线列表 `zip` 起来配 level 标签，光顺前后分支数一旦不同就错配标签（现按 level 逐条遍历）；
- 新增 `buttock_coverage`，把每个 level 的“跨站数 / 分支数 / 每支样本数”写进报告（`fairing_study` 的
  `cases[].buttock_coverage`），第三向不再可能悄悄退化为零。

## 2 方法：追踪而非分段归组

`_trace_branches` 逐站位提取所有穿越点（`_level_crossings`，对分段线性剖面精确反解），
把每个穿越点接到“末端 `zeta` 最近”的开分支上；超过 `_TRACE_JUMP = 0.35`（`zeta` 跨距的比例）视为新分支。
分支在 level set 消失处自然终止（不再用弦连接断点）。

> 为什么不按“落在第几段水线之间”归组：Wigley 21×17 上最顶部 level 中舷落在第 12 段、艏艉落在第 14 段，
> 分段归组会把一条 5 点 buttock 切成 3 + 2（短支被丢弃）。追踪能保持其为一条。

## 3 结果（数字均由 `docs/lines_fairing/buttock_branch_evidence.py` 程序化写出）

三船 21×17 网格，第三向翻转数（曲率符号翻转，`polyline_sign_flips`）：

| 船型 | 旧 level 样本数 | 旧翻转(干净/带噪) | 新 level（峰值比例） | 新样本数 | 新翻转(干净/带噪/光顺后) | 总翻转（光顺前→后） |
|---|---|---|---|---|---|---|
| wigley_iii | [9,11,15,17,17] | 0 / 6 | 0.05, 0.26, 0.47, 0.695, 0.905 | [7,11,15,17,19] | 0 / 6 / **0** | 140 → 13 |
| series60_cb060 | [7,11,13,15,17] | 0 / 0 | 0.05, 0.245, 0.44, 0.635, 0.83 | [7,11,13,17,19] | 0 / 0 / 0 | 180 → 168 |
| kvlcc2_class | [16,17,18,19,19] | 0 / 40 | 0.05, 0.275, 0.5, 0.725, 0.95 | [14,17,18,19,19] | 0 / 36 / **35** | 298 → 293 |

tumblehome 试验船（21×17，`tolerance = 0.032·B/2`）：

| 量 | 旧算法 | 新算法 |
|---|---|---|
| level = 0.5·peak 的 buttock 条数 | 0 | 2（`zeta` 0.3588 / 0.8668，解析双解） |
| 5 个 adaptive level 的翻转（带噪 / 光顺后） | 14 / — | **44 / 26（−41%）** |
| 总翻转（三向） | — | 244 → 121 |
| 偏差 vs 容差 | — | 0.01511 m ≤ 0.01557 m ✔ |

结论：**旧指标在 tumblehome 上完全看不见第三向（level 0.5 时 0 条曲线）**；修复后第三向不仅可见，
且光顺确实把它降了 41%（这是原来根本测不到的改善）。三船算例全部仍满足验收：
`deviation ≤ tolerance`、总翻转非增、确定性可复现（`python3 -m tensorcad.lines_fairing --study` 重跑 results.json 一致）。

## 4 精度判定（honest）

- **等值线反解精度**：对“表格定义的”分段线性曲面，反解在区间内线性插值 —— 与该曲面**逐位精确**
  （`tests/test_lines_fairing.py::test_buttock_contours_are_circular_level_sets_of_a_dome` 断言 `rel=1e-9`）。
- **对连续解析形状的偏差**：tumblehome 面 level 0.5·peak 的下支 `zeta = 0.3588`，解析解 `zeta = 0.3715`，
  差 **0.0127（占吃水 2.0%）** —— 这是 **17 条水线表格本身的离散误差**（`z^(0.7)` 弦线误差 O(h²)），
  不是等值线算法的误差；加密水线该差值按 h² 收敛。
- **翻转计数**：无解析真值可对比（曲率符号翻转数是拓扑量），判定依据是“解析生成的干净表 = 0 翻转”
  （三船解析表新旧两种算法均给 0 ✔）与“带噪表 > 0、光顺后严格下降”。
- **是否 <3%**：曲率能量/翻转数无百分比基准，**不作 <3% 的声称**；几何反解对表格精确、对连续面 2.0%（表格离散，
  非算法）。三船验收条款（偏差在容差内 + 翻转单调非增）**逐项满足**。

## 5 注意点 / 局限

- 顶部 level（0.95·peak）在 kvlcc2 上仍可覆盖 14 点；wigley/series60 顶部被 0.3·N 的覆盖底线压到
  0.905/0.83 peak —— 覆盖度和“贴着甲板”是权衡，现由 `_MIN_COVERAGE=0.3` 显式给出，不再靠运气。
- 近平行相切（level 略低于局部最大半宽）处，穿越 `zeta` 在相邻站位间跳变，追踪会正确切成多条短支；
  短支被 `min_points=3` 丢弃并在 `buttock_coverage` 里可见（不是静默为零）。
- 仍未接入 `geometry_audit` 管线（L3 集成）与 `tensor_benchmark_cases` 判据行（L4）——本单不做，
  留给后续单：建议条目见 finish summary。
- 三船中 kvlcc2/series60 的“解析参照”即生成面本身，外部锚仍缺（与既有报告一致的限制）。

## 6 复现

```bash
cd <worktree>
pip install -e .  # 或 PYTHONPATH=.
python3 -m pytest tests/test_lines_fairing.py -q          # 27 passed
python3 -m tensorcad.lines_fairing --study --stations 21 --waterlines 17 \
        --json docs/lines_fairing/results.json           # exit 0, 3/3 cases meet criterion
PYTHONPATH=. python3 docs/lines_fairing/buttock_branch_evidence.py   # 本文件所有数字 + 两图
```

图：`docs/lines_fairing/buttock_branches.png`（tumblehome 一 level 双穿越 + 旧/新等值线对比）、
`docs/lines_fairing/buttock_coverage.png`（三船覆盖度与翻转数，旧 vs 新）、
`docs/lines_fairing/{wigley_iii,series60_cb060,kvlcc2_class}.png`（光顺前后三族曲线 + 翻转历史）。