# Examples that could not start: a clean-checkout reproducibility gate

- 任务：m16106（selfstart·generalist 席），仓 TensorCAD，分支 `agent/m16106`
- 触发：`audit_report.md`（m14491 质量审计）§4 问题 1「示例不可自举」——审计只抽到 3 个失败脚本
- 可执行证据：`tests/test_examples_bootstrap.py`（69 项，11.13 s 全绿）、
  `artifacts/examples_bootstrap/scan.json`、`artifacts/examples_bootstrap/fig1_clean_checkout_scan.png`
- 结论：干净检出下 **21/48** 个示例入口无法启动（`ModuleNotFoundError`），本单修到 **0/48**，
  并把「不得回归」做成 pytest 机制（静态 + 动态两道检查）

## 1. 案例介绍（问题背景与对象）

TensorCAD 的 README 承诺示例是可执行证据而不是叙述：每个 `examples/<name>/run_*.py` 都要能直接跑通、
并在断言失效时非零退出。这个承诺隐含一个前提——**示例自己负责让 `tensorcad` 可导入**：
仓库根不是已安装包（`pyproject.toml` 未提供安装型入口给示例），而
`python examples/<name>/run_example.py` 启动时 `sys.path[0] == examples/<name>/`，
所以脚本顶部的 `from tensorcad import ...` 会在计算任何数字之前抛 `ModuleNotFoundError`。

对象：仓内全部 48 个示例入口 `examples/*/run_*.py`（不含被示例 import 的辅助模块）。

## 2. 怎么做的（方法 / 判据 / 工具链）

判据先定：**判据 P1** = 在 `cwd=/tmp`、`env -u PYTHONPATH`、未安装 `tensorcad` 的条件下执行
`python3 <script>`，脚本不得以 `ModuleNotFoundError` 失败；**判据 P2** = 该性质由测试守住（不是靠人记）。

测量（零手抄，脚本化）：

```bash
# 静态：哪些脚本完全没有 sys.path 引导
grep -l "sys.path" examples/*/run_*.py
# 动态：干净环境下逐脚本启动，4 s 超时；命中 ModuleNotFoundError 记为 broken
cd /tmp && timeout 4 python3 /abs/path/to/script.py
```

修复：在 21 个失败脚本的模块 docstring（以及 `from __future__` 之后）插入与仓内既有 27 个脚本
完全一致的约定片段：

```python
REPO_ROOT = Path(__file__).resolve().parents[2]
if str(REPO_ROOT) not in sys.path:  # allow running from a checkout without installing
    sys.path.insert(0, str(REPO_ROOT))
```

补丁由 AST 定位插入点（docstring 后、`__future__` 之后），因此不会改变脚本语义，只提前 `sys.path`。

机制（防回归）：`tests/test_examples_bootstrap.py`

1. **静态**（48 参数化用例）：对每个 `examples/*/run_*.py` 用 `ast` 找到第一个 `tensorcad` import 的行号，
   断言其前的源码里出现过 `sys.path` 引导。跑一次 < 1 s，任何新脚本漏写引导立即红。
2. **动态**（21 参数化用例）：把脚本的**模块体**（imports + 定义，`if __name__ == "__main__"` 不执行）
   在子进程里用 `runpy.run_path(..., run_name='tc_example_probe')` 执行，环境剔除 `PYTHONPATH/PYTHONHOME`、
   `cwd` 指到 `tmp_path`。只有真解析得到 `tensorcad` 才能通过——注释和字符串骗不过它。
   选择 runpy 而非整脚本执行，是为了让这道门只考「能否启动」，不承担 30 个数值示例的算力成本。

## 3. 计算结果（数字表，程序化注入）

| 时点 | 干净检出下 `ModuleNotFoundError` 的入口数 | 占比 |
|---|---:|---:|
| 修复前（main @ `2b2499e`） | 21 / 48 | 43.75 % |
| 修复后（本分支） | **0 / 48** | 0.00 % |

| 机制 | 用例数 | 结果 | 墙钟 |
|---|---:|---|---:|
| `tests/test_examples_bootstrap.py` | 69 | 69 passed | 11.13 s |

修复前失败的 21 个入口（`artifacts/examples_bootstrap/scan.json`，机器写入）：
`analytic_hydrostatics/run_analytic_check.py`、`asymmetric_hull/run_example.py`、`bonjean/run_example.py`、
`csg_boolean/run_example.py`、`damage/run_example.py`、`fairness/run_example.py`、`freeboard_seam/run_example.py`、
`geometry_core/run_voxel_export.py`、`gz/run_example.py`、`hydrostatics_verification/run_verification.py`、
`intact_stability/run_example.py`、`nurbs_degree_elevation/run_example.py`、`nurbs_edge_alignment/run_example.py`、
`nurbs_piecewise_elevation/run_example.py`、`nurbs_surface_degree_elevation/run_example.py`、`planing/run_example.py`、
`righting_arm/run_example.py`、`series60_verification/run_verification.py`、`stability/run_example.py`、
`transom_stern/run_example.py`、`wigley_stability/run_example.py`

## 4. 结果图

`artifacts/examples_bootstrap/fig1_clean_checkout_scan.png`：修复前 / 修复后「干净检出下无法启动的入口数」柱状图，
虚线为 48 个入口总数（= 满分线）。数据直接来自本次扫描，无手工填数。

## 5. 注意点（踩过的坑 / 局限 / 适用范围）

1. **修复只是把门推开，门后面还有东西**：去掉 `ModuleNotFoundError` 后，两个示例暴露出**真实的 API 漂移**，
   本单未修（如实登记，留给下一席）：
   - `examples/asymmetric_hull/run_example.py:54`：`analytic_hydrostatics(..., fore_exponent=0.9)`
     → `TypeError: unexpected keyword argument 'fore_exponent'`（`hydrostatics_analytic` 只接受 `breadth_exponent`）。
     疑点：`tests/test_asymmetric_hull.py` 也在传 `fore_exponent`，但该测试文件当前未在本单执行范围内，未判真伪。
   - `examples/series60_verification/run_verification.py:82`：`model.validation.watertight`
     → `AttributeError: 'ValidationReport' object has no attribute 'watertight'`（`validation.py` 只暴露
     `checks/errors/warnings/passed`）。
   两者都说明「示例可启动」≠「示例可跑完」，判据 P1 的边界就在这里，不要把它当成示例已可用。
2. **动态门只覆盖 21 个曾失败的入口**，其余 27 个仅受静态门约束：它们历史上就能启动，但同样可能因
   API 漂移跑不完。全量「跑完」扫描是下一个最小可证单元（见 §6）。
3. **runpy 探针会执行模块体**，故对「模块级就做重计算」的脚本等价于跑一遍；本单 21 个脚本模块体都在
   定义阶段，总耗时 11 s，可接受。若将来有示例在模块级跑算例，应把它降级为 `__main__` 内的调用，
   而不是放宽这道门。
4. 适用范围：只针对「仓库根未安装」这一运行形态。`pip install -e .` 后这些脚本同样可跑，所以这是
   **开发/审计侧**的缺陷，不是终端用户缺陷。

## 6. 误差判定

| 关键量 | 基准 | 修复前 | 修复后 | 是否 < 3 % |
|---|---|---:|---:|---|
| 干净检出下能启动的示例入口占比 | 48/48 = 100 %（README 承诺） | 56.25 %（缺 43.75 %） | **100 %**（缺 0 %） | ✅ 0 % |

- 本任务无物理场量，故不设文献/解析基准；唯一可度量的关键量是「入口可启动率」，其目标值由 README 的
  可执行证据承诺给定（100 %），误差按缺失比例计。
- **诚实边界**：本单只宣称「21 个曾失败入口现在能启动，且该性质被测试守住」；
  **不宣称**全部 48 个示例都能跑完（见 §5.1 的两个真实失败），也不宣称数值示例的物理结论有效。
- 机制本身的有效性有反例保护：修复前扫描命中 21 个 broken，修复后同一扫描命中 0 个——
  同一判据、同一命令，前后对照。

## 7. 给下一席的线索

1. 修 §5.1 的两个 API 漂移（各 1 处调用），然后把 `tools/check_examples.py`（全量跑完 + 超时分级）
   接进 CI，把判据从 P1 升级到 P2「示例跑完且断言成立」。
2. `audit_report.md` §6 剩余的审计整改项仍然开放：`geomquality/pipeline.py:200` 的 `workdir` 死分支
   （`pass` 静默吞参数）、`high_order_geometry.py` 与 `highorder/*` 双实现、`work/` 过程残渣归档。