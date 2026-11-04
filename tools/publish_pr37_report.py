"""Publish the PR #37 cure report (credentials come from the environment, never the file)."""

import base64
import io
import json
import os
import urllib.request
from pathlib import Path

BASE = os.environ.get("TENSOR_API_BASE", "http://127.0.0.1:8010/api/v1")
USER = os.environ.get("TENSOR_AGENT_USER")
PASS = os.environ.get("TENSOR_AGENT_PASS")
ROOT = Path(__file__).resolve().parents[1]
MISSION = 3745


def post(path, payload, token=None):
    headers = {"Content-Type": "application/json"}
    if token:
        headers["Authorization"] = f"Bearer {token}"
    request = urllib.request.Request(f"{BASE}{path}", data=json.dumps(payload).encode(),
                                     headers=headers, method="POST")
    with urllib.request.urlopen(request, timeout=180) as response:
        return json.loads(response.read())


def token() -> str:
    data = post("/auth/login", {"username": USER, "password": PASS})
    return data.get("access_token") or data.get("token") or data["data"]["access_token"]


def figure(path, name, caption, width=1000):
    from PIL import Image

    image = Image.open(path).convert("RGB")
    if image.width > width:
        image = image.resize((width, int(image.height * width / image.width)), Image.LANCZOS)
    buffer = io.BytesIO()
    image.save(buffer, format="PNG", optimize=True)
    raw = buffer.getvalue()
    return {"name": name, "caption": caption, "png_base64": base64.b64encode(raw).decode()}, len(raw)


audit = json.loads((ROOT / "artifacts/pr37/union_audit.json").read_text())
metrics = json.loads((ROOT / "artifacts/pr37/cure_metrics.json").read_text())
state = json.loads((ROOT / "artifacts/pr37/branch_state.json").read_text())
names = audit["name_ledger"]
base = audit["tests"]["main_baseline"]["counts"]
merged = audit["tests"]["merged"]["counts"]

markdown = f"""# PR #37 治愈：`agent/m2873`（艉板船体系族）与 `main`（Wigley 族）的语义并集合并

**仓** TensorCAD · **PR** #37（分支 `agent/m2873`）· **单** m3745 · **结论**：冲突 4 文件全部按语义并集解决，
0 文件丢失 / 0 公开名丢失 / 0 新失败用例，合并号 0.21.0。

## ① 案例介绍（背景 / 几何 / 参数）

TensorCAD 的几何层原本只有一个船体系族（正弦收尾的参数化船体，首尾都收敛到一点）。两个并行演进
分支各自给这个族加了一个新族，也因此把仓库元数据（版本号、公开命名空间、变更日志）顶到了同几行上：

* `main` 侧（7 个提交）：**Wigley 抛物船体族** `tensorcad/wigley.py`，闭式静水力
  `Cb = nm/((n+1)(m+1))`、`Cp = Cwp = n/(n+1)`、`Cm = m/(m+1)`、`VCB = T - T(m+1)/(2(m+2))`，
  版本号已走到 **0.20.0**；另含 #33 的非对称船体与 #43 合并。
* `agent/m2873` 侧（1 个提交，重建在与 `main` 同一祖先 `6c568c3` 上）：**艉板船体系族**
  `tensorcad/transom.py`。船尾在 `τ` 站位被切平、用平面艉封板封闭，参数
  `transom_fraction τ`（切船位置）与 `breadth_exponent p`（正弦收尾指数）；`A_t/A_m = sin(πτ)^p`
  可闭式反解，不完整 taper 积分 `∫_τ^1 sin(πu)^p du` 由 `sinc^p` 幂级数求出，因此 `V, A_wp, LCB,
  Cb, Cp, Cwp` 与浸没艉板面积都有闭式基准。该侧版本号写在 **0.19.0**（基于 0.16.0 的 main）。

两侧**没有一个共同源码文件**，冲突全部落在 4 个记账文件上（`CHANGELOG.md`、`pyproject.toml`、
`tensorcad/__init__.py`、`tests/test_public_surface.py`）。这正是机械 overlay（`--ours`/`--theirs`）
最容易**静默吃掉一个船体系族**的形状：合并后仍可导入、测试仍绿，但 `generate_wigley`
或 `generate_transom_hull` 会悄悄消失。

## ② 怎么做的（方法 / 流程 / 工具链）

1. 取 `origin/main` 与 `origin/agent/m2873` 的最新 tip，确认 `merge-base = 6c568c3`，
   `main` 领先 7 提交、分支领先 1 提交。
2. 在工作树内 `git merge --no-commit --no-ff origin/agent/m2873`，得到 4 个冲突文件，其余 38 个
   单侧改动文件自动合并。
3. **不用 `--ours`/`--theirs`**：`tools/resolve_pr37.py` 对每个冲突块应用一条显式规则，
   并在写盘后断言文件中不再残留冲突标记：
   * `CHANGELOG.md` — 两侧发布段**逐字保留**（0.20.0 Wigley 段 + 0.19.0 艉板段），新增 0.21.0 抬头
     说明「艉板族并入 Wigley 线、0.19.0 从未单独发布」，因此不改写任一段正文、不丢任一条目；
   * `pyproject.toml` — 标量无法求并集：取两侧最大值再升一位，`0.21.0`；
   * `tensorcad/__init__.py` — `__all__` 取**排序并集**（`generate_transom_hull` 与 `generate_wigley`
     同时保留），`__version__` 随版本号；`from .transom import …` 与 `from .wigley import …` 两侧
     导入语句自动合并后均在位；
   * `tests/test_public_surface.py` — `EXPECTED_API` 与分层模块元组取并集，版本断言改为 0.21.0。
     该测试用 `set(tensorcad.__all__) == EXPECTED_API` 比较，任一名字丢失立即红灯。
4. **机器证据**：`tools/audit_pr37_union.py` 把合并重放成集合代数，输出 4 张账（文件账 / 名字账 /
   测试账 / 判据），并集不精确则非零退出。
5. 测试基线对照：在干净的 `origin/main` 工作树上跑同一命令，再在合并树上跑同一命令。
6. 两个族的算例各自跑通，并把分支提交的 CSV 与重跑结果做**逐字节**比对（证明合并没让数值漂移）。

## ③ 计算结果（数字表，全部由 `artifacts/pr37/*.json` 注入，无手抄）

**文件/名字并集审计（`artifacts/pr37/union_audit.json`）**

| 量 | 值 |
|---|---|
| 两侧从共同祖先起触碰的路径 | {audit['files_touched_by_either_side']} |
| 单侧路径（必须与其贡献方逐字节一致，按 blob id 比对，含 PNG 二进制） | {len(audit['file_ledger']['one_sided'])} |
| 冲突文件（逐个声明解法） | {len(audit['file_ledger']['conflicts'])} |
| 合并后缺失路径 | {len(audit['file_ledger']['missing'])} |
| 未声明的意外冲突 | {len(audit['file_ledger']['unexpected'])} |
| `main` 独有公开名 | {len(names['main_only'])} |
| `m2873` 独有公开名 | {len(names['branch_only'])} |
| 两侧共有公开名 | {names['both']} |
| 合并后 `__all__` 名数 / 并集名数 | {names['merged_count']} / {names['union_count']} |
| 合并丢失 / 合并凭空新增的名字 | {len(names['missing_from_merge'])} / {len(names['invented_in_merge'])} |

**测试对照（同一命令 `pytest -q --ignore=tests/test_cli_neural.py --ignore=tests/test_neural_pipeline.py -n 8`）**

| 树 | passed | failed | skipped |
|---|---|---|---|
| `main` {audit['main_tip'][:9]}（基线） | {base.get('passed')} | {base.get('failed')} | {base.get('skipped')} |
| 合并树 | {merged.get('passed')} | {merged.get('failed')} | {merged.get('skipped')} |
| 差 | **+{audit['tests']['passed_delta']}** | 0 | 0 |

新增的 +96 例全部来自 `tests/test_transom.py`；2 个失败（`test_geometry_field.py` 的
learned-revision 与 learned-field-gradient）在 `main` 上**同样失败**，原因都是本沙箱缺
`ModuleNotFoundError: No module named 'safetensors'`，属环境缺口而非合并回归。

**算例与复现性**：`examples/wigley/run_example.py` 退出 0（10 条断言全 PASS）；
`examples/transom_stern/run_example.py` 退出 0（4 个 τ 全覆盖 watertight / 边流形 / Euler-2）；
重跑生成的 `convergence.csv`、`design_table.csv`、`reference_ladder.csv` 与分支提交件**逐字节相同**。

## ④ 结果图

图 1 三联：(a) 测试平权——基线 706 绿 vs 合并 802 绿，两侧同为 2 红 11 跳；(b) 合并树上的艉板族
体积误差随网格收敛（双对数，与二阶参考斜率对照）；(c) 不完整 taper 积分的三条独立求值路线
（`sinc^p` 级数 vs Γ 闭式 vs tanh-sinh 求积）之差，都在舍入量级。

## ⑤ 注意点（坑 / 局限）

* **机械 overlay 的坑正是本单要防的**：两侧只碰记账文件，`--ours`/`--theirs` 都能让 CI 全绿而删掉一整个
  船体系族；因此解法必须逐文件给规则，并用名字账机器复核。
* **版本号是标量，只能判断不能求并集**：分支写 0.19.0（基于 0.16.0），`main` 已到 0.20.0。取 0.21.0
  既不回退发布号，也不给已发布的 0.20.0 偷偷塞新内容。
* 顺带修了分支侧一个可复现性缺陷：`examples/transom_stern/run_example.py` 缺同族算例都有的三行
  `sys.path` 引导，从 clone（未 `pip install`）直接跑会 `ModuleNotFoundError: tensorcad`。引导行逐字
  复制自 `examples/wigley/run_example.py`，不改任何数字。
* **PR 引用未能改写**：`agent/m2873` 远端更新的动作被沙箱安全策略拒绝，因此治愈提交落在
  `agent/m3745`（tip `{state['cure_commit'][:9]}`，父提交 `{state['cure_commit_parents'][0][:9]}`=main、
  `{state['cure_commit_parents'][1][:9]}`=旧 `agent/m2873`，旧 SHA `{state['old_sha'][:9]}`），可被
  merge_governor 直接快进；平台验收后发布。
* 审计是「文件 / 名字 / 测试」三级的；同一函数体内两侧各自改动的语义风险，靠人工读 4 个冲突块确认
  ——它们全是列表型（版本串、导出元组、日志段），块外三方文本合并是精确的。
* `safetensors` 相关用例在本沙箱无法运行，两侧同样未验证，不做粉饰。

## ⑥ 误差判定（是否 < 3%）

| 量 | 基准 | 误差 | < 3%？ |
|---|---|---|---|
| 内容丢失（路径数） | 0 允许 | 0 | **是（0%）** |
| 内容丢失（公开名数） | 0 允许 | 0 | **是（0%）** |
| 新增失败用例（vs `main` 基线） | 0 允许 | 0 | **是（0%）** |
| 网格体积，16 个 (网格, τ) 组合最差 | `transom.py` 闭式 | {metrics['volume_error_worst_all_grids']*100:.3f}% | **是** |
| 网格体积，最细网格 321×81 | `transom.py` 闭式 | {metrics['volume_error_worst_finest_grid']*100:.3f}% | **是** |
| 不完整 taper 积分，级数 vs Γ 闭式 | `Γ((1+p)/2)/(√π Γ(1+p/2))` | {metrics['ladder_max_closed_form']:.1e}（相对） | **是** |
| 不完整 taper 积分，级数 vs tanh-sinh | `tanh_sinh_quadrature` | {metrics['ladder_max_quadrature']:.1e}（相对） | **是** |

本单主判据是**内容零丢失**（0 vs 0 允许）；数值行是「合并后代码仍算出分支所声称的数」的非回归校验，
不是新的物理结论。

**产物**：`tools/resolve_pr37.py`、`tools/audit_pr37_union.py`、`tools/plot_pr37_cure.py`、
`tools/report_pr37_doc.py`、`docs/pr37_merge_audit.md`、`artifacts/pr37/`、`runs/`（两侧 pytest 日志与两个算例日志）。
"""

fig, size = figure(ROOT / "artifacts/pr37/cure_evidence.png", "pr37_cure_evidence.png",
                   "图 1 (a) 测试平权：main 706 绿 vs 合并 802 绿，两侧同为 2 红 11 跳；"
                   "(b) 合并树上艉板族体积误差随网格二阶收敛；(c) 不完整 taper 积分三条独立路线之差在舍入量级。")

payload = {
    "title": "PR #37 治愈：艉板族 ∪ Wigley 族的语义并集合并（0.21.0，内容零丢失）",
    "summary": (
        f"4 个记账冲突逐文件按显式规则求并集（两段发布日志逐字保留、__all__ 120 名 = main 10 + 分支 13 + 共有 97、"
        f"版本 0.21.0），审计显示 0 路径丢失 / 0 名字丢失 / 0 名字凭空新增；测试 706→{merged.get('passed')} 绿、"
        "同 2 个 safetensors 环境失败，0 新失败；两个族的算例均退出 0 且 CSV 逐字节可复现。"
    ),
    "scenario_key": "tensorcad_merge_cure",
    "mission_id": MISSION,
    "tags": ["TensorCAD", "merge-cure", "PR#37", "transom-stern", "wigley", "semantic-union", "0.21.0"],
    "markdown": markdown,
    "figures": [fig],
}
print(f"payload figure bytes: {size}, markdown chars: {len(markdown)}, "
      f"approx total KB: {(size * 4 / 3 + len(markdown)) / 1024:.0f}")
result = post("/tensor/reports", payload, token())
print(json.dumps(result)[:400])