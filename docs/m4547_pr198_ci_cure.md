# PR #198 (agent/m4362) — CI 红修复记录（m4547）

## 1. 现象

GitHub 上 PR #198（分支 `agent/m4362`，head `83ecf5e`）的 `tests` workflow 矩阵：

| job | 结论 |
|---|---|
| test (3.10) | cancelled |
| test (3.11) | cancelled |
| test (3.12) | cancelled |
| test (3.13) | failure |

矩阵默认 `fail-fast: true`：一个 job 失败后其余在跑的 job 被取消，所以
「1 failure + 3 cancelled」= **同一确定性故障在 4 个解释器上都会复发**，
不是 3.13 专有缺陷。要定位的是一条与 Python 版本无关、且在 `pytest` 之前执行的步骤。

`tests.yml` 的步骤顺序为：

1. `pip install -e ".[test,dev]"`
2. **Lint**：`ruff check .`
3. Type check：`mypy`
4. Tests：`pytest --cov --cov-report=term-missing`（`fail_under = 90`）
5. CLI / 示例 / 可复现性若干步

## 2. 根因

本地按 CI 环境复现（Python 3.11.6，`ruff 0.15.0`，`mypy` strict，pytest 9.2.0）：

```
$ python3 -m ruff check . --output-format concise
tests/test_merge_cure_m3391.py:78:28: W292 [*] No newline at end of file
tools/make_report_m3859.py:21:34: RUF100 [*] Unused `noqa` directive (unused: `E402`)
tools/make_report_m3859.py:22:21: RUF100 [*] Unused `noqa` directive (unused: `E402`)
tools/make_report_m3859.py:40:15: RUF059 Unpacked variable `rc` is never used
tools/make_report_m3859.py:132:15: UP031 Use format specifiers instead of percent format
tools/make_report_m3859.py:163:68: W292 [*] No newline at end of file
tools/verify_merge_m3391.py:30:8: F401 [*] `math` imported but unused
tools/verify_merge_m3391.py:197:38: B905 `zip()` without an explicit `strict=` parameter
tools/verify_merge_m3391.py:371:21: W292 [*] No newline at end of file
tools/verify_merge_m3859.py:124:9: F841 Local variable `lost` is assigned to but never used
tools/verify_merge_m3859.py:178:21: W292 [*] No newline at end of file
Found 11 errors.
```

`pyproject.toml` 里 dev 依赖是 `ruff>=0.6`（无上限），CI 装的是最新 ruff（本地 0.15.0）。
新 ruff 的点亮规则（`RUF059`、`UP031`、`B905`、`RUF100`、`W292`）命中了
m3391/m3859 遗留的工具脚本与一个测试文件，**与 PR #198 的新内容（`tools/wetted_uncertainty*.py`、
`tests/test_wetted_uncertainty.py`、`data/wetted-uncertainty-2026-11-04/`）无关**——
同样的 11 条在 `main` 上也存在，属"共享祖先上的历史欠账"，PR 只是被迫继承。

判据链（逐条排除）：

- `mypy`：`Success: no issues found in 34 source files` — 不是类型检查步骤。
- `pytest`（遮掉 torch，复刻 CI 不装 `neural` extra 的环境；`test_cli_neural.py`
  与 `test_neural_pipeline.py` 在缺 torch 时整体 `importorskip` 跳过）：
  `907 passed, 12 skipped in 79.26s` — 测试步骤本身不是红的来源。
- PR 新增的 `tests/test_wetted_uncertainty.py`：`11 passed` — 交付物自身健康。
- 因此红点收敛到第 2 步 `ruff check .`。

## 3. 修复

只做行为等价的 lint 清理（4 个文件、9 增 11 删），不触碰 PR 的领域代码与数据：

| 文件 | 位置 | 处理 |
|---|---|---|
| `tests/test_merge_cure_m3391.py` | EOF | 补文件末尾换行（W292） |
| `tools/make_report_m3859.py` | 21–22 | 删掉已失效的 `# noqa: E402`（RUF100） |
| `tools/make_report_m3859.py` | 40 | `name, rc, summary` → `name, _rc, summary`（RUF059） |
| `tools/make_report_m3859.py` | 132 | `"..." % box["limit_deg"]` → f-string（UP031） |
| `tools/make_report_m3859.py` | EOF | 补换行（W292） |
| `tools/verify_merge_m3391.py` | 30 | 删未使用的 `import math`（F401） |
| `tools/verify_merge_m3391.py` | 197 | `zip(curve.gz, reference, strict=True)`（B905） |
| `tools/verify_merge_m3391.py` | EOF | 补换行（W292） |
| `tools/verify_merge_m3859.py` | 124 | 删死变量 `lost: list[str] = []`（F841，实际统计走 `losses`） |
| `tools/verify_merge_m3859.py` | EOF | 补换行（W292） |

`zip(..., strict=True)` 是有意选的方向：`curve.gz` 与解析解 `reference` 逐角对应，
长度不一致说明数据构造出错，宁可显式报错也不静默截断。

## 4. 验证

```
$ python3 -m ruff check . --output-format concise
All checks passed!

$ python3 -m mypy
Success: no issues found in 34 source files

$ python3 -m pytest -q -n 16 --ignore=tests/test_cli_neural.py --ignore=tests/test_neural_pipeline.py
907 passed, 12 skipped, 19 warnings in 79.26s

$ python3 -m pytest -q tests/test_merge_cure_m3391.py tests/test_wetted_uncertainty.py
12 passed, 2 skipped
```

`ruff` 由 11 errors → 0；`mypy` 与测试通过数不变，说明改动是纯 lint 层面的等价整理。

## 5. 局限与注意点

- **覆盖率门禁未在本机复跑完整值**：本机 `coverage` 缺 C tracer（`coverage.tracer` 无法加载），
  带 `--cov` 的整套跑要 10 分钟以上，超出本地单次命令窗口。PR 未改动 `tensorcad/` 下任何一行，
  覆盖率分子分母均不变，故 `fail_under = 90` 的结果与 main 一致；此点如实标注为"未直接实测"。
- **ruff 未固定上限**是这次红测的结构性原因。彻底解法是在 dev extra 里给 ruff 加上限
  （如 `ruff>=0.6,<0.16`）或升级时统一跑一次清理；本轮按"修代码不藏问题"处理，
  建议下一轮由平台侧决定是否加 pin。
- 修复落在 `tools/`、`tests/` 的历史文件上，属于"让 PR 的 CI 绿所必需"的最小改动，
  不改变任何船舶/几何算法行为。