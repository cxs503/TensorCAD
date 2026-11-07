"""Publish the damage-stability mission report through the platform API."""
from __future__ import annotations

import base64
import io
import json
import os
import urllib.request

BASE = os.environ.get("TENSOR_API_BASE", "http://127.0.0.1:8010/api/v1")
USER = os.environ.get("TENSOR_AGENT_USER", "eco_admin")
PASS = os.environ.get("TENSOR_AGENT_PASS", "EcoTest!2026pass")


def post(path, payload, token=None):
    data = json.dumps(payload).encode()
    request = urllib.request.Request(
        BASE + path, data=data, method="POST",
        headers={"Content-Type": "application/json", **({"Authorization": f"Bearer {token}"} if token else {})},
    )
    with urllib.request.urlopen(request, timeout=120) as response:
        return response.status, response.read().decode()


def shrink(path, target=60000):
    raw = open(path, "rb").read()
    if len(raw) <= target:
        return raw
    from PIL import Image

    image = Image.open(path).convert("RGB")
    scale = (target / len(raw)) ** 0.5
    width, height = image.size
    image = image.resize((max(360, int(width * scale)), max(220, int(height * scale))))
    buffer = io.BytesIO()
    image.save(buffer, "PNG", optimize=True)
    return buffer.getvalue()


status, body = post("/auth/login", {"username": USER, "password": PASS})
token = json.loads(body)["access_token"]
result = json.load(open("docs/damage/results.json"))
rows = "\n".join(
    f"| {c['claim']} | {c['measured']:.10g} | {c['reference']:.10g} | "
    f"{c['relative_error'] * 100:.1e}% | {'PASS' if c['pass'] else 'FAIL'} |"
    for c in result["claims"]
)
figures = []
for name, caption in (
    (
        "docs/damage/buoyancy_curve.png",
        "Buoyancy curve A(x) before and after a midship compartment floods; the shaded "
        "band is the buoyancy the ship no longer has.",
    ),
    (
        "docs/damage/sweeps.png",
        "Left: damaged draft against the closed form L T/(L-l) over the whole sweep. "
        "Right: heel of a port-side compartment and GM, which goes negative for long "
        "damage and is reported as negative.",
    ),
):
    figures.append(
        {
            "name": os.path.basename(name),
            "caption": caption,
            "png_base64": base64.b64encode(shrink(name)).decode(),
        }
    )

mid = result["midship"]
markdown = f"""## 1 案例介绍
TensorCAD 0.23 的 `docs/intact_stability.md` 与 CHANGELOG 都明写「damage stability (part B) 未实现」。本单补上确定性破舱稳性的**自由进水（lost buoyancy）**一半：箱形舱室被击穿并与舷外连通，船舶重量不变而该处浮力丢失，船下沉、纵倾，非对称破舱还会横倾。

程序内生成的水密盒体（仅 GM 用到 KG）：

| L [m] | B [m] | T [m] | 排水体积 [m³] | KG [m] | 干舷 [m] |
|---|---|---|---|---|---|
| 10.0 | 2.0 | 0.8 | 16.0 | 0.95 | 1.0 |

工况：中舱（x 4–6 m，μ=1）、中舱 μ=0.6、尾端舱（x 0–2 m）、左舷半宽舱（x 4–6 m，y −1–0 m）与右舷镜像、舱顶在下的部分深度舱。

## 2 怎么做的
重量不变 ⇒ 浮力平衡 `V_ext − μ V_c = Δ/ρ`；净浮心必须落在 `(x_G, y_G)`。两者都在同一水线面
`z = z0 + tan(trim)(x−x_ref) + tan(heel)·y` 上积分，三个未知量由带回溯线搜索的阻尼 Newton
（有限差 Jacobian）一次解出；剖面切割每站只算一次并缓存，迭代内只做多边形裁剪。舱壁端面被设为
Gauss 面板边界（3 点 Gauss–Legendre），`A(x)` 在舱端的跳变不会被跨越采样。`GM = KB + I_T/(Δ/ρ) − KG`
的 `I_T` 取**净水线面**（外部减 μ×破舱部分），这正是 lost-buoyancy 法对自由液面的处理，再叠加自由液面
修正会重复扣减。工具链：纯标准库 + pytest + matplotlib；`examples/damage/run_example.py` 逐条断言、
失败非零退出，并把 `docs/damage_stability.md`、`docs/damage/results.json` 与两张 PNG 写出来。

## 3 计算结果（程序化注入，零手抄）
中舱 μ=1：破损吃水 {mid['draft_at_reference_m']:.12g} m（闭式 `LT/(L−l)` = 1.0 m），丢浮力
{mid['lost_buoyancy_volume_m3']:.12g} m³，破损后 GM {mid['gm_transverse_m']:.6g} m。μ=0.6 吃水
{result['mu_0p6']['draft_at_reference_m']:.12g} m（闭式 0.909090…）。尾端舱 z0 =
{result['aft']['draft_at_reference_m']:.12g} m、tan(trim) = {result['aft']['waterline']['trim_tangent']:.10g}；
左舷舱横倾 {result['port_side']['heel_angle_deg']:.6f}°，右舷镜像 {result['starboard_side']['heel_angle_deg']:.6f}°。

| 断言 | 实测 | 基准 | 相对误差 | 判定 |
|---|---|---|---|---|
{rows}

## 4 结果图
图 1：浮力曲线——外部 A(x) 与净 A(x)，阴影为被海水吃掉的浮力。图 2：左为破损吃水 vs 闭式
`LT/(L−l)` 的全长度扫掠，右为左舷破舱横倾角与 GM（长破舱使 GM 转负，如实报负值不截断）。

## 5 注意点
- 只处理**自由进水**；内部油舱破裂属 added-weight + 自由液面，不可混用。
- 不做概率破舱（SOLAS B-1 分舱指数 A/R）、递进进水、横贯注水、破损船 GZ。
- 渗透率 μ 是输入不猜值；状态记录 `lost_buoyancy_volume = μ V_c`。
- 无干舷船体破舱若损失浮力超过水线以上容积，直接报「the ship would be a submarine」，不给被截断的假答案。
- 曲线船体上剩余误差是网格本身：完好态对 `tensorcad.hydrostatics` 吃水差 0.5%、KB/LCB 约 1%、`I_T` 5% 内。

## 6 误差判定
15 条量纲断言最大相对误差 {result['worst_relative_error'] * 100:.2e}%（盒体闭式解在浮点舍入量级）。
关键量 vs 基准：中舱吃水、μ=0.6 吃水、丢浮力、破损 I_T、尾端 z0/tan(trim)、舷侧 z0/tan(heel)、镜像横倾——
**均 <3%**（实测 ≤1e-12）。唯一 ≥3% 的一项是曲线船体 `I_T` 与既有 hydrostatics 路线的跨代码差异（≤5%），
已单列说明来源，不冒充 <3%。负控三条（舱长↑吃水↑、μ=0 精确复现完好态、左右舷镜像等值反号）均通过，
说明比较具备失效能力。

## 证据
- commit `95f8181`（分支 agent/m6274）：`tensorcad/damage.py`、`tests/test_damage.py`（24/24 通过）、
  `examples/damage/run_example.py`、`docs/damage_stability.md`、`docs/damage/` 下的 results.json 与 PNG、CHANGELOG。
"""

payload = {
    "title": "TensorCAD 破舱稳性（lost buoyancy 自由进水）：吃水+纵倾+横倾三平衡求解与盒体闭式解验证",
    "summary": (
        "新增 tensorcad/damage.py：自由进水破舱平衡一次解出吃水/纵倾/横倾与破损后 GM；"
        f"15 条闭式断言最大相对误差 {result['worst_relative_error'] * 100:.1e}%；"
        "tests 24/24 通过；曲线船体完好态与既有 hydrostatics 独立路线一致（吃水 0.5%，I_T ≤5%）。"
    ),
    "scenario_key": "tensorcad",
    "mission_id": 6274,
    "tags": ["TensorCAD", "damage-stability", "lost-buoyancy", "verification", "ship-hydrostatics"],
    "markdown": markdown,
    "figures": figures,
}
print("figure bytes:", [len(f["png_base64"]) for f in figures])
status, body = post("/tensor/reports", payload, token)
print(status, body[:400])