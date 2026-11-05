# PR #75 CI 治愈（m4087）：三个 mission 工具的 lint 红 → 绿，报告正文零字节改动

**仓** TensorCAD（`agent/m3778`，PR #75） · **单** m4087（修复·第 1 轮） · **红项** `test (3.11)=failure`，其余矩阵项 cancelled
**产物** 本文件 + `tools/m3778_publish_report.py`、`tools/m3778_pr41_figures.py`、`tools/verify_merge_pr41.py` 的 lint 治愈

## ① 案例介绍

PR #75 是对 `main` 提出的 8 文件 / +954 行的交付（`git diff --stat origin/main...HEAD`，合并基 `2ff80b2`）：

| 文件 | 行数 | 作用 |
|---|---|---|
| `tools/verify_merge_pr41.py` | +216 | PR #41 合并治愈的机读证据生成器（22 个保内容探针） |
| `tools/m3778_pr41_figures.py` | +167 | 两张证据图（公共面保留矩阵 + 合并后树重跑 GZ 物理） |
| `tools/m3778_publish_report.py` | +175 | 平台报告 payload 生成器（数字从 artifacts 注入，不手抄） |
| `docs/m3778_pr41_merge_cure.md` | +160 | 治愈报告 |
| `artifacts/m3778/*` | 4 文件 | 证据 JSON + 两张 PNG |

三个脚本是 PR #75 的全部可执行代码。`git merge-tree --write-tree origin/main HEAD` 无冲突（exit 0，无 conflict 行），
因此「CI 绿」的唯一障碍就是它们自身的静态检查。

## ② 怎么做的（复现 → 定位 → 最小治愈 → 等价性证明）

1. **工作树**：从 `agent/m3778` 切出（`git switch agent/m3778`，head `273e772`），
   与本机 Python 3.11.6 组合，正好是红掉的那个矩阵项。
2. **复现**：`python -m ruff check .` → **41 errors**，全部落在上述三个脚本，零落点在 `tensorcad/`、`tests/`、`examples/`。
   CI `tests.yml` 的步骤顺序是 Lint → Type check → pytest → 各 example 步骤，第一个步骤就红，
   所以 3.10/3.12/3.13 只是 fail-fast 被 cancel，**不是**四个各自独立的故障。
3. **错误构成**（41 = 26 E501 + 6 UP031 + 3 W292 + 2 C416 + 2 F401 + 1 B905 + 1 RUF100）：
   E501 全部是中文报告正文的长行；其余是 `%` 格式化、缺尾换行、无用 import、`zip` 缺 `strict=`、多余 `noqa`。
4. **治愈策略**：不动一个字的内容，只动源码形状。
   - 长正文行按**隐式字符串拼接**切分（`"…"` + 换行 + `"…"`，拼接不插入任何字符）；切点只允许落在
     `{…}` 插值之外与转义序列之外，且按 **显示宽度**（东亚宽字符计 2）判定，而不是 `len()`。
   - 因为 `return f"""…"""` 后面跟续行会变成**两条语句**（见 ⑤ 第一坑），整块包进 `return ( … )`。
   - `%` → f-string / `{}`，`zip(..., strict=True)`，`json.dumps(dict(gz))`，删无用 import 与多余 `noqa`，补尾换行。
5. **等价性证明**（这一步是治愈的验收判据）：把治愈前的三个脚本副本留在 `/tmp`，用同一份合成 evidence 调
   `markdown(ev, gz)`，比较治愈前后渲染出的报告正文：

   | | 字符数 | sha256 |
   |---|---|---|
   | 治愈前 | 3794 | `3ab01d1656d8d3ace9fd4b0f959d5bbd3f5283ee027edeb10253b98d2b4e9b47` |
   | 治愈后 | 3794 | `3ab01d1656d8d3ace9fd4b0f959d5bbd3f5283ee027edeb10253b98d2b4e9b47` |

   字节级相同 ⇒ 已发布的报告正文（平台 mission 3778 的报告）不受本次改动影响。

## ③ 计算结果（程序化，命令即证据）

| 判据 | 治愈前 | 治愈后 | 说明 |
|---|---|---|---|
| `ruff check .` | 41 errors | **All checks passed!** | CI 红掉的第一步，机器判定 |
| `mypy`（strict，32 source files） | success | **Success: no issues found** | CI 第二步 |
| 测试（12 文件） | — | **331 passed / 17 warnings / 20.8 s** | `-n 16` |
| 测试（14 文件） | — | **329 passed, 10 skipped / 24.6 s** | `-n 16` |
| 合计观测 | — | **660 passed, 10 skipped, 0 failed** | 见 ⑤「未测部分」 |
| `tools/verify_merge_pr41.py` 实跑 | — | 正常输出 22 个探针 `[PASS]/[FAIL]` | 格式化改动未破坏工具 |
| `git merge-tree --write-tree origin/main HEAD` | — | exit 0，0 处冲突 | PR 可常规合并 |
| 治愈后最长源码行 | 357 | **110**（= 项目 line-length） | 全仓 `ruff` 通过 |

## ④ 结果图

**本任务不产图，理由**：改动是静态检查层面的源码形状 + 一个 byte-for-byte 的正文等价性证明，
没有任何场量、曲线或几何量可画。唯一的量化证据是 ③ 的表格（lint 计数、测试计数、sha256）。
强行画「错误数柱状图」不属于场量云图，也不属于与解析解/基线的比较曲线，故按规矩声明无图。

## ⑤ 注意点（踩过的坑 / 局限 / 适用范围）

- **坑 1（险些静默改坏报告）**：`return f"""A""" ` 换行 `f"""B"""` 在缩进相同的函数体里是**两条语句**，
  第二条在 `return` 之后永不执行 —— 语法合法、ruff 全绿、`markdown()` 却只返回 155 字符（原本 3794）。
  是等价性比对（sha256）把它抓出来的；没有那个比对，这会是一次「CI 绿 + 报告变空」的假治愈。
  正解：整块包进 `return ( … )`。
- **坑 2**：ruff 的 E501 按**显示宽度**计数（实测 `"# " + "中"*60` → 122 > 110），
  按 `len()` 切分会「切完还红」。切点选择必须用 `unicodedata.east_asian_width`。
- **坑 3**：`{{`/`}}`（f-string 里的字面花括号）不能当插值括号计数，否则 `artifacts/{{…}}` 段落找不到合法切点。
- **坑 4**：`ruff check --fix` 会把没有占位符的续行 `f"""` 降级成 `"""`（F541，17 处自动修复）；
  对正文无影响（已验证），但人工复核时容易被误读为「字符串被改了」。
- **局限 1（如实记）**：本沙箱单命令 120 s 上限 + 机器 load average ≈1900（128 核被其他 mission 占满），
  **全量套件与 `--cov`（`fail_under = 90` 覆盖率门）未能本地跑完**：带 `--cov` 的全量跑到 ~50%、不带 cov 跑到 ~96% 均被
  120 s 预算收割。已跑完的 660 项 0 失败，未跑部分为 `test_benchmark / test_cli / test_design / test_extract /
  test_neural_* / test_suboff / test_suboff_appendages / test_mesh_conditioning`（重算型，本机单文件即超预算）。
  本次改动只碰 `tools/`，而覆盖率口径 `source = ["tensorcad"]`，故覆盖率数值**按构造不变**；
  但「覆盖率仍 ≥90」这一条我**没有本地实测**，交给 CI 复核，不粉饰。
- **局限 2**：`tools/m3778_pr41_figures.py` 的 `gz_panel` 在本分支无法执行（`ModuleNotFoundError: tensorcad.heel`）——
  该脚本是为 `agent/m3213` 的合并树写的，PR #41 的 `heel` 模块从未进入本分支。同理 `verify_merge_pr41.py`
  在本分支只 14/22 探针 PASS（缺的都是 heel 侧内容）。二者都是**既有内容缺口**，不是本次 lint 治愈引入的，
  也不在 CI 执行路径上（CI 不跑 `tools/`），故不在本单范围，原样保留并在报告里留证。
- **适用范围**：本单只做 lint 治愈，不动任何几何/数值/证据逻辑；PR #75 的物理内容（公共面并集、GZ 曲线）未触碰。

## ⑥ 误差判定

| 关键量 | 基准 | 实测 | 误差 | 是否 <3% |
|---|---|---|---|---|
| lint 错误数 | 0（CI 判据） | 0 | 41 → 0，精确归零 | **是**（严格 0） |
| 报告正文 sha256 | 治愈前渲染文本 | 相等（3794/3794 字符） | 0 字符 | **是**（严格 0） |
| 已跑测试 | 0 failed | 660 passed / 0 failed | 退化率 0 | **是**（严格 0） |
| 合并冲突 | 0 | `merge-tree` exit 0 | 0 | **是**（严格 0） |
| mypy | 无错误 | strict 32 文件无错误 | — | **是** |
| 覆盖率门（≥90%） | — | **未实测** | — | **无基准，不判定**（原因见 ⑤ 局限 1） |

**结论**：CI 红项（Lint）已本地复现并治愈，副作用经字节级比对证明为零；覆盖率与全量套件交由 CI 复核。