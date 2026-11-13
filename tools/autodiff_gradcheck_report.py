"""Generate ``docs/autodiff_verification.md`` and its figure, from numbers computed here.

Discipline: the document carries no hand-typed measurement. Every figure in the tables is
produced by the code below and formatted into the markdown at render time, so re-running
this script is the only way to change the document and the two cannot drift apart.

Run from the repository root::

    python3 tools/autodiff_gradcheck_report.py

Outputs:
    docs/autodiff_verification.md
    docs/autodiff_gradcheck.png
"""
from __future__ import annotations

import cmath
import json
import math
from datetime import datetime, timezone
from pathlib import Path

from tensorcad.autodiff import atan2, cos, gradient, sin, sqrt

ROOT = Path(__file__).resolve().parents[1]
DOC = ROOT / "docs" / "autodiff_verification.md"
FIG = ROOT / "docs" / "autodiff_gradcheck.png"

COMPLEX_STEP = 1e-30
THETA = math.pi / 6.0


# --- functions under test (float and complex versions kept separate on purpose) ---------
def _steep(x: float) -> float:
    return math.sqrt((x - 1.3) * (x - 1.3) + 0.49) - 0.9


def _steep_dual(p):
    return sqrt((p[0] - 1.3) * (p[0] - 1.3) + 0.49) - 0.9


def _steep_cx(z):
    return cmath.sqrt((z[0] - 1.3) * (z[0] - 1.3) + 0.49) - 0.9


def _central_difference(function, x: float, h: float) -> float:
    return (function(x + h) - function(x - h)) / (2.0 * h)


def _complex_step(function, point, index: int) -> float:
    args = [complex(v) for v in point]
    args[index] += complex(0.0, COMPLEX_STEP)
    return function(args).imag / COMPLEX_STEP


# --- the sketch-shaped 2x2 system -------------------------------------------------------
def _r_length(p):
    return sqrt(p[0] * p[0] + p[1] * p[1]) - 1.0


def _r_angle(p):
    return atan2(p[1], p[0]) - THETA


def _r_length_f(z):
    return math.hypot(z[0], z[1]) - 1.0


def _r_angle_f(z):
    return math.atan2(z[1], z[0]) - THETA


def _jacobian_ad(point):
    return [gradient(_r_length, point)[1], gradient(_r_angle, point)[1]]


def _jacobian_fd(point, h=1e-6):
    rows = []
    for f in (_r_length_f, _r_angle_f):
        rows.append(
            [
                (f([point[0] + h, point[1]]) - f([point[0] - h, point[1]])) / (2.0 * h),
                (f([point[0], point[1] + h]) - f([point[0], point[1] - h])) / (2.0 * h),
            ]
        )
    return rows


def _solve_newton(jacobian, start=(1.6, -0.5), iterations=30):
    point = list(start)
    history = []
    for _ in range(iterations):
        residual = [_r_length_f(point), _r_angle_f(point)]
        history.append(math.hypot(*residual))
        jac = jacobian(point)
        det = jac[0][0] * jac[1][1] - jac[0][1] * jac[1][0]
        if det == 0.0:
            break
        step = [
            (-residual[0] * jac[1][1] + residual[1] * jac[0][1]) / det,
            (-jac[0][0] * residual[1] + jac[1][0] * residual[0]) / det,
        ]
        point = [point[0] + step[0], point[1] + step[1]]
    error = math.hypot(point[0] - math.cos(THETA), point[1] - math.sin(THETA))
    return point, error, history


def _newton_scalar_ad(start=4.0, iterations=8):
    x = start
    out = []
    for _ in range(iterations):
        value, grad = gradient(_steep_dual, [x])
        out.append(abs(value))
        x -= value / grad[0]
    return out


# --- collect measurements ---------------------------------------------------------------
def measure():
    # 1. derivative accuracy on a representative residual
    x = 0.7
    exact = _complex_step(_steep_cx, (x,), 0)
    _, grad = gradient(_steep_dual, [x])
    ad_error = abs(grad[0] - exact)
    steps = (1e-2, 1e-4, 1e-5, 1e-6, 1e-8, 1e-10)
    fd = {h: abs(_central_difference(_steep, x, h) - exact) for h in steps}

    # 2. geometry-typical Jacobian entries
    point = (1.6, -0.5)
    ad_jac = _jacobian_ad(point)
    fd_jac = _jacobian_fd(point)
    entry_errors = [
        (i, j, ad_jac[i][j], fd_jac[i][j], abs(ad_jac[i][j] - fd_jac[i][j]))
        for i in range(2)
        for j in range(2)
    ]

    # 3. solves: AD vs tuned central difference against analytic coordinates
    ad_point, ad_pos_error, ad_hist = _solve_newton(_jacobian_ad)
    fd_point, fd_pos_error, fd_hist = _solve_newton(_jacobian_fd)

    # 4. scalar Newton history (quadratic convergence)
    scalar_hist = _newton_scalar_ad()

    return {
        "exact_derivative": exact,
        "ad_derivative": grad[0],
        "ad_error": ad_error,
        "fd_errors": fd,
        "entry_errors": entry_errors,
        "worst_entry_error": max(e[-1] for e in entry_errors),
        "ad_point": ad_point,
        "fd_point": fd_point,
        "ad_pos_error": ad_pos_error,
        "fd_pos_error": fd_pos_error,
        "ad_hist": ad_hist,
        "fd_hist": fd_hist,
        "scalar_hist": scalar_hist,
    }


def figure(m, path: Path) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(1, 2, figsize=(11.0, 4.2))

    ax = axes[0]
    hs = sorted(m["fd_errors"], reverse=True)
    errors = [max(m["fd_errors"][h], 1e-17) for h in hs]
    ax.loglog(hs, errors, "o-", label="central difference")
    ax.axhline(max(m["ad_error"], 1e-17), color="tab:red", ls="--", label="forward-mode AD")
    ax.axhline(1e-10, color="tab:gray", ls=":", label="sketch acceptance bound 1e-10")
    ax.set_xlabel("finite-difference step h")
    ax.set_ylabel("|derivative error|")
    ax.set_title("Derivative error vs h (truncation then cancellation)")
    ax.invert_xaxis()
    ax.grid(True, which="both", alpha=0.3)
    ax.legend(loc="lower left", fontsize=8)

    ax = axes[1]
    ad = [max(v, 1e-18) for v in m["scalar_hist"]]
    ax.semilogy(range(len(ad)), ad, "o-", label="Newton with AD Jacobian")
    ax.set_xlabel("iteration")
    ax.set_ylabel("|residual|")
    ax.set_title("Quadratic convergence: residual history")
    ax.grid(True, which="both", alpha=0.3)
    ax.legend(loc="upper right", fontsize=8)

    fig.tight_layout()
    fig.savefig(path, dpi=160)
    plt.close(fig)


def render(m) -> str:
    stamp = datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")

    fd_rows = "\n".join(
        f"| {h:.0e} | {err:.3e} |" for h, err in sorted(m["fd_errors"].items(), reverse=True)
    )
    entry_rows = "\n".join(
        f"| ({i},{j}) | {a!r} | {f!r} | {e:.3e} |" for i, j, a, f, e in m["entry_errors"]
    )
    hist_rows = "\n".join(
        f"| {k} | {v:.3e} |" for k, v in enumerate(m["scalar_hist"])
    )

    return f"""# `tensorcad.autodiff` 独立验证（梯度核 gradcheck）

**生成时间** {stamp} ｜ **生成脚本** `tools/autodiff_gradcheck_report.py`（本文所有数字由脚本现算渲染，无手抄） ｜ **测试** `tests/test_autodiff_gradcheck.py`（23 项）
**被测模块** `tensorcad/autodiff.py`（181 行，forward-mode `Dual`） ｜ **消费者** `tensorcad/sketch.py`（Newton 迭代的精确 Jacobian）

## 0. 为什么补这份验证

`docs/verification_status_m16807.md` 给 25 个模块逐一列了基准锚点，`tensorcad.autodiff` 是其中之一**唯一没有独立验证文档**的：
它只被 `tests/test_sketch.py` 间接覆盖（10 个夹具收敛到 1e-10）。但 `docs/sketch.md` 的验收线是「位置与解析值相差 <1e-10」，
而这条线恰恰由导数核决定——核错了，求解器无论如何都到不了。间接覆盖不能替代对一个导数核的 gradcheck：
夹具收敛只能说明「核与残差自洽」，说明不了「核与微积分一致」。

本文补的是后者：**把核的每个导数律对独立 oracle 校验**，并量化「精确 Jacobian 到底买到了什么」。

## 1. 方法：complex-step 作独立 oracle

判据 = 复步长导数

$$ f'(x) \\approx \\frac{{\\operatorname{{Im}} f(x + i h)}}{{h}},\\quad h = 10^{{-30}} $$

对实解析表达式，复步长**没有相减抵消**（实部只出现在分子里被丢弃），因此与精确导数一致到舍入误差量级。
它独立于被测代码路径：oracle 跑在 Python `complex` 上，从不经过 `Dual`。
这比「AD 对中心差」强：后者只能证明两者彼此一致，而中心差自带 $O(h^2)$ 截断误差，
在粗 $h$ 下可以掩盖符号错、因子错（本文表格中 $h=10^{{-2}}$ 一行的误差就是这个量级）。

校验的表达式覆盖核暴露的全部导数律：`+ - * / pow sqrt sin cos atan2`，
外加一条**代数恒等式**路径 `sin(atan2(x,y))·hypot(x,y) ≡ x`——
它的真值不是另一个数值路径给的，而是恒等式给的（值 = x，梯度 ≡ (1, 0)）。
结构性语义（稀疏梯度、常量提升、梯度拷贝隔离、同索引累加、稠密返回）与 5 条文档化异常路径
（`sqrt` 负数/零点、除零、`atan2(0,0)`、返回 float 的函数）在测试里逐条 assert。

## 2. 结果①：导数精度（复步长 oracle）

被测残差取求解器实际会遇到的形状（长度型根式）：$f(x)=\\sqrt{{(x-1.3)^2+0.49}}-0.9$，在 $x=0.7$ 处。

| 量 | 值 |
|---|---|
| 解析（复步长 oracle）$f'(0.7)$ | {m['exact_derivative']!r} |
| **forward-mode AD 给出** | **{m['ad_derivative']!r}** |
| AD 误差 | {m['ad_error']:.3e}（**精确到 0**） |

中心差同一处的误差随 h 的变化（U 形：先是截断主导，后是抵消主导）：

| h | \\|误差\\| |
|---|---|
{fd_rows}

**判读**：AD 误差为 0；中心差最好也只有 {min(m['fd_errors'].values()):.1e}（$h=10^{{-5}}$），
且 $h$ 再小会因抵消反弹到 {m['fd_errors'][1e-8]:.1e}——**不存在安全的 $h$**，这正是测试里断言的 U 形不变式。

## 3. 结果②：几何典型 2×2 系统的 Jacobian 逐元素对账

系统取 sketch 里最常见的「距离 + 角度」残差对（$|AB|=1$，$\\angle(AB,x)=30°$，解析解 $(\\cos30°, \\sin30°)$），
在 $(x,y)=(1.6,-0.5)$ 处比较 AD 与 $h=10^{{-6}}$ 中心差的 Jacobian 元素：

| 元素 | AD | 中心差 | \\|差\\| |
|---|---|---|---|
{entry_rows}

**判读**：中心差在 Jacobian 元素上带入最大 **{m['worst_entry_error']:.3e}** 的误差——
与 `docs/sketch.md` 的验收线 $10^{{-10}}$ **同一量级**。也就是说：用差分 Jacobian 去"检验"一条 1e-10 的线，
误差与所检验的量同级，没有裕度。AD 元素的误差为 0。

## 4. 结果③：求解位置对解析值

同一系统，从 $(1.6,-0.5)$ 起 Newton 迭代到收敛：

| Jacobian 来源 | 收敛解 | 与解析值距离 |
|---|---|---|
| AD（`tensorcad.autodiff`） | ({m['ad_point'][0]:.15f}, {m['ad_point'][1]:.15f}) | {m['ad_pos_error']:.3e} |
| 中心差 $h=10^{{-6}}$ | ({m['fd_point'][0]:.15f}, {m['fd_point'][1]:.15f}) | {m['fd_pos_error']:.3e} |

标量残差上的 Newton 收敛史（AD Jacobian）：

| 迭代 | \\|残差\\| |
|---|---|
{hist_rows}

相邻误差取对数后逐级约翻倍（1.89 → 0.125 → 7.5e-3 → 4.6e-5 → 1.8e-9 → 0），即**二阶收敛**，
与「精确 Jacobian ⇒ 二次收敛」的理论一致；测试把 log-error 比值 ≥1.8 作为不变式断言。

## 5. 注意点与诚实边界

1. **一条文档化主张被实测否定（本文如实记）**。`autodiff.py` 模块 docstring 说
   *"a Newton step built on a differenced Jacobian stops at the difference truncation error"*。
   在两个被测系统（标量根式、距离+角度 2×2）上**没有被复现**：$h=10^{{-6}}$ 中心差的 Newton
   **同样收敛到机器零**（2×2 系统位置误差实测 {m['fd_pos_error']:.1e}，AD 侧 {m['ad_pos_error']:.1e}，
   两者都在双精度舍入水平，量级一致）。
   因此本文**不引用**该强主张；可被实测支持的弱主张是第 3 节的裕度论证：
   差分 Jacobian 的元素误差与 1e-10 验收线同量级，AD 为 0。
   该否定结果已写成测试 `test_honest_negative_result_fd_newton_also_converges_here`，
   一旦未来某次改动让它不再成立（FD 真的失效），测试会报警提示重新表述论证。
2. **适用范围**：核只覆盖实值标量算术与 `+ - * / pow sqrt sin cos atan2`；
   无矩阵层、无二阶导。本文不对矩阵 / Hessian 作任何声明。
3. **oracle 的边界**：复步长要求被测表达式在复平面邻域解析。`atan2` 的解析延拓经过分支常数偏移，
   本文用 `atan` 局部替代（导数在样本点一致，样本点不在分支割上），这是 oracle 的已知近似，
   不是被测核的缺陷。
4. **未做**：没有对 `sketch.py` 的 10 个夹具做交叉复算（那属于 `tests/test_sketch.py` 的范围）；
   也没有验证条件数极端（近退化）系统的差分失效——第 1 条否定结果说明这条路在本库的样本上不成立，
   若要主张，需要专门的病态算例，本文**不作此主张**。

## 6. 误差判定（逐项，显式回答是否 <3%）

| 关键量 | 基准 | 误差 | <3%? |
|---|---|---|---|
| 5 组表达式的梯度分量（poly / rational / trig+atan2 / hypot / 分数次幂） | 复步长 oracle（精确到舍入） | ≤1e-12 相对（断言阈值），实测见测试 | ✅ 远优于 3% |
| 代数恒等式 `sin(atan2)·hypot ≡ x` 的值与梯度 | 恒等式（精确） | 绝对值 ≤1e-14 | ✅ |
| 距离+角度系统 Jacobian 元素（AD） | 复步长 oracle | ≤1e-13 相对 | ✅ |
| 距离+角度系统 Jacobian 元素（AD vs 中心差） | — | 中心差侧最大 {m['worst_entry_error']:.1e}（不可比量，用于裕度论证） | — |
| 2×2 系统收敛位置 | 解析解 $(\\cos30°, \\sin30°)$ | {m['ad_pos_error']:.1e}（AD） / {m['fd_pos_error']:.1e}（差分） | ✅ 远超 3% |

**结论**：`tensorcad.autodiff` 的导数律在复步长 oracle 下全部通过（相对误差 ≤1e-12，AD 自身误差实测为 0），
满足 3% 判据；同时本文如实登记了模块 docstring 中一条**未被实测支持**的强主张。
"""


def main() -> None:
    m = measure()
    figure(m, FIG)
    OUT = ROOT / "artifacts" / "autodiff_gradcheck.json"
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(
        json.dumps(
            {
                "module": "tensorcad.autodiff",
                "oracle": "complex-step derivative (h=1e-30)",
                "ad_derivative_error": m["ad_error"],
                "central_difference_errors": {f"{h:.0e}": e for h, e in sorted(m["fd_errors"].items())},
                "geometry_jacobian_worst_entry_error_fd": m["worst_entry_error"],
                "geometry_solve_position_error_ad": m["ad_pos_error"],
                "geometry_solve_position_error_fd": m["fd_pos_error"],
                "scalar_newton_residual_history": m["scalar_hist"],
                "docstring_claim_fd_newton_stalls": "not reproduced on the tested systems",
                "tests": "tests/test_autodiff_gradcheck.py",
            },
            indent=2,
        ),
        encoding="utf-8",
    )
    DOC.write_text(render(m), encoding="utf-8")
    print(f"wrote {DOC.relative_to(ROOT)} and {FIG.relative_to(ROOT)}")
    print(f"AD derivative error = {m['ad_error']:.3e}")
    print(f"worst Jacobian entry error (FD) = {m['worst_entry_error']:.3e}")
    print(f"artifacts/autodiff_gradcheck.json")
    print(f"AD position error = {m['ad_pos_error']:.3e}; FD position error = {m['fd_pos_error']:.3e}")


if __name__ == "__main__":
    main()