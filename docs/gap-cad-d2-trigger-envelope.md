# 几何金标门的触发包络：白名单漏掉的 69 个模块

**单号**: m22749（GAP-CAD-D2-GOLDEN-REGRESSION 的验证欠账）
**对象**: `.githooks/pre-commit-geometry-golden.sh`（本地臂）与 `.github/workflows/geometry-golden.yml`（CI 臂）
**结论**: 两臂原来用**手写模块白名单**（11 个）决定"谁算几何改动"。实测四个金标 case 测量时执行
**80 个** `tensorcad` 模块，白名单漏 69 个；其中 `geometry_field.py` 的 +2% 半径改动会让套件判红，
但提交它时钩子返回 0、什么都不打印。已改为按包前缀/通配覆盖，并把"跑过的模块集合"做成每次运行现推导的测试。

## 1. 复现（全部只读命令，除注入副本外不改仓）

| 步骤 | 命令 | 结果 |
|---|---|---|
| 基线门（改前，HEAD `f90fe85`） | `python -m tensorcad.golden --check --json /tmp/g0.json` | 退出 0，4 case，68.0 s |
| 运行时导入迹 | `python -m tensorcad.golden --check --quick` 后读 `sys.modules` | **80** 个 `tensorcad` 模块 |
| 旧白名单覆盖 | 上面的集合减旧白名单 11 项 | **69 个未覆盖**（`loft.py`、`geometry_field.py`、`csg.py`、`nurbs.py`、`analytic_solids.py`…） |
| 旧钩子对漏网模块 | `GOLDEN_STAGED_FILE=<tensorcad/geometry_field.py> sh .githooks/pre-commit-geometry-golden.sh` | 退出 **0**，无输出（门未启动） |
| 材质性：真回归 | 副本 `/tmp/mut_gf`，`SphereField.numpy_callable` 的 `radius = 1.02 * self.radius`，`--check --quick` | 退出 **1**，5 项 FAIL，球体积 `4.42972951` vs 参考 `4.1887902` = **+5.75%**（容差 0.546%） |
| 新钩子同一暂存内容 | 同上命令 | 退出 **1**，`[golden] geometry-affecting change detected: tensorcad/geometry_field.py` |

原始数字：`artifacts/gap-cad-d2-golden/trigger_envelope.json`（程序化写入，未手抄）。

## 2. 改了什么

| 文件 | 改动 |
|---|---|
| `.githooks/pre-commit-geometry-golden.sh` | `WATCHED='^(tensorcad/(hull\|wigley\|…))\\.py)'` → `WATCHED='^tensorcad/'` |
| `.github/workflows/geometry-golden.yml` | `paths` 的 11 行模块清单 → `tensorcad/**`，并补 `docs/data/**`（数据编辑同样能移动金标数字） |
| `tests/test_geometry_golden_gate.py` | 删除"两臂白名单逐名相等"的测试（它把缺陷写成了契约）；新增：整包静态覆盖、**运行时导入迹现推导覆盖**、以及"旧白名单漏掉的模块现在必须触发"的回归探针 |
| `tools/golden/gap_d2_gate_docs.py` | 按新的 `WATCHED` 形式解析触发集合 |
| `docs/gap-cad-d2-geometry-golden-gate.md` | §2 触发路径表与说明同步 |

## 3. 验证

- `python -m pytest tests/test_geometry_golden_gate.py -q` → **14 passed, 11 failed**；11 个失败全部是
  依赖 `artifacts/gap-cad-d2-golden/golden_suite.json` 的既有测试，该目录在本 checkout 里 **0 个受版本
  管理文件**（`git ls-files artifacts/gap-cad-d2-golden | wc -l` = 0），与本次改动无关，改前同样失败。
- 改后 `tensorcad.golden --check`（全量，非 quick）在 HEAD 代码上仍退出 0（几何数字未变：本次只动
  触发集合与测试，不动任何测量/参考代码）。

## 4. 局限（诚实）

- **没跑通第二条臂本身**：dual-derived 套件在 TensorPapers，本 worktree 旁没有该 checkout，
  `suite.py` / `mutants.py` 本轮未执行（`not_run` 已写入证据 JSON）。
- **`docs/gap-cad-d2-geometry-golden-gate.md` 的生成器 `--check` 在本 checkout 跑不起来**：它要读
  `artifacts/gap-cad-d2-golden/golden_suite.json`，该快照在本分支不存在。§2 因此是**按生成器同一套
  解析逻辑手工同步**的（表体与新 `WATCHED` 一致），不是 `--check` 绿过的产物。
- `tensorcad/**` 会多触发：任何包内改动都跑一遍（粗档 ~11 s / 全量 ~68 s）。这是刻意选的方向
  ——多触发一次是几秒，漏触发一次就没有门。
- 触发集合覆盖的是"跑过的模块"。**动态 import / 运行期才加载的路径**（`importlib` 字符串导入）在
  导入迹里会登账，但若某段代码只有换个 case 才加载，本迹只能覆盖现有四个 case 的实际路径。
