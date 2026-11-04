"""Publish the m2873 report (credential-free script: reads env for API base and login)."""
import base64
import io
import json
import os
import urllib.request

BASE = os.environ.get("TENSOR_API_BASE", "http://127.0.0.1:8010/api/v1")
USER = os.environ.get("TENSOR_AGENT_USER")
PASS = os.environ.get("TENSOR_AGENT_PASS")


def post(path, payload, token=None):
    headers = {"Content-Type": "application/json"}
    if token:
        headers["Authorization"] = f"Bearer {token}"
    request = urllib.request.Request(
        f"{BASE}{path}", data=json.dumps(payload).encode(), headers=headers, method="POST"
    )
    with urllib.request.urlopen(request, timeout=180) as response:
        return json.loads(response.read())


def login():
    data = post("/auth/login", {"username": USER, "password": PASS})
    return data.get("access_token") or data.get("token") or data["data"]["access_token"]


def figure(path, name, caption, width=900):
    from PIL import Image

    image = Image.open(path).convert("RGB")
    if image.width > width:
        image = image.resize((width, int(image.height * width / image.width)), Image.LANCZOS)
    buffer = io.BytesIO()
    image.save(buffer, format="PNG", optimize=True)
    raw = buffer.getvalue()
    encoded = base64.b64encode(raw).decode()
    return {"name": name, "caption": caption, "png_base64": encoded}, len(raw)


MARKDOWN = """# 艉板船体系族：把「刀尖」船尾换成真实艉封板，并保住闭式静水力基准

**仓** TensorCAD（分支 `agent/m2873`，v0.14.0） · **单** tp-mod-4 自找活（models 席）
**产物** `tensorcad/transom.py`、`tests/test_transom.py`（96 例）、`examples/transom_stern/run_example.py`、
`docs/transom_stern.md`、`outputs/transom_stern/`

## ① 案例介绍（物理背景 / 几何 / 参数）

参数化船体 `b(x,ζ) = (B/2)·sin(πx/L)^0.7·√ζ` 首尾两端都收敛到一点：闭合性满分、可验证性最好，
但**不是一条船尾**。排水型船在 Fn≳0.25 普遍带**浸没艉板**（外壳被切平、尾封板直接浸水），
这块板先于一切改变湿表面积、压力场与兴波，是下游 BEM/阻力码看到的第一个几何特征。
缺口账 `backlog_tensorcad_section.md` 将其列为 **P0「非对称船体 + 艉板」**；本单完成其中**艉板**一半。

| 符号 | 字段 | 含义 |
|---|---|---|
| `L` | `length` | **总长 LOA，艉板到首柱** |
| `B` / `T` | `beam` / `draft` | 水线最大型宽 / 吃水 |
| `τ` | `transom_fraction` | 切船位置（占母体比例），`0 ≤ τ < 0.5` |
| `p` | `breadth_exponent` | 正弦收尾指数，`0.25 ≤ p ≤ 3.0` |
| `nx`,`nz` | — | 站数 / 垂向层数（默认 41×13） |

参考算例 `L=10 m, B=2 m, T=0.8 m, p=0.7, τ=0.2, ρ=1025`（母体长 `L_parent = LOA/(1−τ) = 12.5 m`）。
坐标同 `hull.py`：`x` 由**艉板所在的 x=0** 指首，`z=0` 设计水线，龙骨 `z=−T`。

## ② 怎么做的（方法 / 流程 / 关键设置）

1. **拓扑**：每站一个闭合剖面环（右舷水线→龙骨→左舷水线，龙骨点共用），站间放样成四边形，水线盖条封闭；
   `τ>0` 时用**平面扇形**封住艉剖面。剖面在龙骨以上为凸（抛物线龙骨剖面 + 平水线），扇形无需新增顶点、
   无需三角化库。`τ=0` 跳过艉板即精确退回尖尾船体（面表完全一致）。
2. **保住闭式基准（本单技术核心）**：截断把完整积分 `I=∫₀¹sin(πu)^p du` 变成不完全积分
   `J(τ)=∫_τ¹sin(πu)^p du`，一般 `p` 无初等闭式。做法是拆项 `sin(πs)^p=(πs)^p·(sinc(πs))^p`：
   前项是唯一奇点、逐项积分精确，后项对 `s` 是偶函数即 `s²` 的幂级数（阶乘收敛），于是
   `∫₀^a s^m sin(πs)^p ds = π^p Σ_k c_k a^{m+p+2k+1}/(m+p+2k+1)`（24 项）。基准因此**仍与网格无关**。
3. **三重独立校验**：`τ=0` 必须复现 Γ 闭式 `Γ((1+p)/2)/(√πΓ(1+p/2))`；首矩必须满足
   `M(0)/J(0)=1/2`（尖尾 LCB 在舯）；`τ>0` 必须与仓内已有的双指数积分 `tanh_sinh_quadrature` 一致。
4. **设计侧反解**：`A_t/A_m = sin(πτ)^p` → `τ = asin(r^(1/p))/π`，设计环可直接给"10% 浸没艉板"。
5. 工具链：标准库核心 + pytest 9.2 / matplotlib 3.10。
   `run_example.py` 程序化写出全部 CSV/PNG，本报告数字零手抄。

## ③ 计算结果

**基准自校验**（`reference_ladder.csv`）：级数 vs Γ 闭式最大相对差 **2.2e-16**；
级数 vs tanh-sinh **1.1e-15**。

**网格积分 vs 闭式基准**（`convergence.csv`，排水体积相对误差）：

| τ | 41×13 | 81×25 | 161×41 | 321×81 |
|---:|---:|---:|---:|---:|
| 0.00 | −0.880 % | −0.308 % | −0.136 % | −0.048 % |
| 0.05 | −0.813 % | −0.286 % | −0.129 % | −0.045 % |
| 0.20 | −0.785 % | −0.278 % | −0.126 % | −0.045 % |
| 0.35 | −0.770 % | −0.274 % | −0.125 % | −0.044 % |

水线面积 / LCB 更快：41×13 最差 −0.175 % / −0.065 %；321×81 为 −0.005 % / −0.002 %。
艉板浸没面积（平面扇形 vs 解析剖面）：−0.71 %(nz=13) / −0.12 %(nz=41) / −0.043 %(nz=81)，单调。

**设计空间**（`design_table.csv` 节选）：

| p | A_t/A_m | τ | Cb | Cp = Cwp | LCB/LOA |
|---:|---:|---:|---:|---:|---:|
| 0.7 | 0.05 | 0.0044 | 0.4715 | 0.7073 | 0.4979 |
| 0.7 | 0.15 | 0.0212 | 0.4784 | 0.7177 | 0.4905 |
| 1.5 | 0.10 | 0.0691 | 0.3965 | 0.5947 | 0.4653 |
| 3.0 | 0.15 | 0.1783 | 0.3386 | 0.5079 | 0.3989 |

拓扑不变量（τ ∈ {0, 0.01, 0.05, 0.2, 0.35, 0.49}）：水密 True、边流形 True、`V−E+F=2`、外向、无零面积面。

## ④ 结果图

图 1 `convergence.png`：体积与 LCB 误差随站距的双对数收敛（四档 τ），全部朝基准单调下坠——
自refinement只能显示数字在动，有独立基准才画得出这张图。图 2 `geometry.png`：左=纵剖线 + 艉板位置（红），
右=横剖线（含 τ=0.2 的艉封板剖面）。图 3 `design_space.png`：`Cb` 随 τ 走向（p=0.7/1.5/3.0）：
`p` 是 `Cb` 的主旋钮（0.31–0.48），`τ` 只在 τ≲0.2 前有效，之后趋平。

## ⑤ 注意点（坑 / 局限 / 适用范围）

- **`Cp ≡ Cwp`**：族内中横剖面恒定、水线面恰为剖面分布×2/3，两系数恒等；**不要用它研究 Cp/Cwp 解耦**
  （真实船靠平行中体与首尾不对称）。
- **τ≥0.5 被构造性拒绝**：切到舯前会削掉最大剖面，`Cm=2/3` 及全部基准随之失效。
- **`p=0.7` 时"真实"艉板极度靠后**：要 5–15 % 浸没艉板面积比，τ 只有 0.005–0.02（10 m 船上切掉 5–20 cm）。
  这是正弦收尾本身太瘦的物理事实，也是这类参数化族长期被记为"无艉板"的原因；要显眼的艉板需抬 `p`。
- 艉板垂直、单站位；无倾斜艉板 / 轴包 / 球鼻艏；**只做几何 + 静水力，不做阻力**。
- 收敛阶 ≈1.4–1.5（非 2）：`√ζ` 龙骨剖面斜率无穷，与尖尾船体同源。默认 41×13 的 −0.88 % 适合看趋势，
  不适合作 3 % 验收。
- 湿表面积可从网格取（含艉板）但本族**无闭式**；`hydrostatics_analytic` 只覆盖尖尾船体。

## ⑥ 误差判定（关键量 vs 基准，显式回答是否 <3%）

| 关键量 | 基准 | 41×13 | 161×41 | 321×81 | <3% ? |
|---|---|---:|---:|---:|---|
| 排水体积（τ=0…0.35） | Γ/级数闭式 | −0.88 % | −0.14 % | −0.05 % | **是**（四档网格全部） |
| 水线面积 | 同上级数 | −0.18 % | −0.02 % | −0.005 % | **是** |
| LCB（自艉板） | 同上级数（含首矩） | −0.065 % | −0.007 % | −0.002 % | **是** |
| 艉板浸没面积 | 解析剖面 `(2/3)BT sin(πτ)^p` | −0.71 % | −0.12 % | −0.043 % | **是** |
| 基准自身 | 级数 vs Γ / vs tanh-sinh 两独立法 | 2.2e-16 / 1.1e-15 | — | — | **是**（远超 3%） |

结论：**所有关键量在所有测试网格上均 <3 %**，161×41 起全部 ≤0.14 %。无 honest-FAIL 项。
未做：与实验/文献的阻力对标（本单仅几何 + 静水力，见 ⑤）。

## 证据指针

commit `8a1abf2`（分支 `agent/m2873`）· `python3 -m pytest tests/test_transom.py -q` → **96 passed** ·
`tests/test_public_surface.py` → 6 passed（公共面登记 13 个新符号）· 文档 `docs/transom_stern.md` ·
数据与图 `outputs/transom_stern/`。
"""

FIGURES = (
    ("outputs/transom_stern/convergence.png", "convergence.png",
     "体积与 LCB 相对误差随站距收敛（τ=0/0.05/0.2/0.35），基准为闭式解"),
    ("outputs/transom_stern/geometry.png", "geometry.png",
     "τ=0.2 船体纵剖线与艉板（红）及横剖线（含艉封板剖面）"),
    ("outputs/transom_stern/design_space.png", "design_space.png",
     "块系数 Cb 随 τ 的走向：p=0.7/1.5/3.0"),
)


def main():
    token = login()
    figures, sizes = [], []
    for path, name, caption in FIGURES:
        entry, size = figure(path, name, caption)
        figures.append(entry)
        sizes.append((name, size))
    payload = {
        "title": "艉板船体系族：真实艉封板 + 保住闭式静水力基准（TensorCAD 0.14.0）",
        "summary": (
            "TensorCAD 参数化船体首尾刀尖、缺真实艉板（backlog P0）。新增 tensorcad/transom.py：τ 站截断 + "
            "平面扇形封艉，τ=0 精确退回尖尾船体；不完全正弦积分用 sinc^p 幂级数算到机器精度，并由 Γ 闭式"
            "(2.2e-16) 与 tanh-sinh (1.1e-15) 双重独立校验；96 例测试、4τ×4 网格收敛（体积 0.88%→0.045%）、"
            "设计反解 τ=asin(r^(1/p))/π、三张图。"
        ),
        "scenario_key": "tensorcad_transom_stern_family",
        "mission_id": 2873,
        "tags": ["TensorCAD", "hull-form", "transom-stern", "hydrostatics", "verification",
                 "closed-form-reference", "m2873"],
        "markdown": MARKDOWN,
        "figures": figures,
    }
    print("payload bytes:", len(json.dumps(payload).encode()), "figures:", sizes)
    result = post("/tensor/reports", payload, token)
    print("published:", json.dumps(result)[:500])


if __name__ == "__main__":
    main()
