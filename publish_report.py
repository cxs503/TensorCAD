import base64, json, os, urllib.request as u
from pathlib import Path

out = Path("docs/gz")
data = json.loads((out / "results.json").read_text())
claims = data["claims"]
barge = data["barge"]
round_bilge = data["round_bilge"]
rows = barge["table"]
round_rows = round_bilge["table"]
imo = barge["imo"]

def fmt(value, digits=6):
    return f"{value:.{digits}f}"

def wall_worst():
    entry = [c for c in claims if c["claim"].startswith("wall-sided")][0]
    return entry

wall = wall_worst()
round_claim = [c for c in claims if c["claim"].startswith("round bilge: max |GZ")][0]
float_claim = [c for c in claims if c["claim"].startswith("round bilge: max |float")][0]
gm_claim = [c for c in claims if c["claim"].startswith("hull:")][0]

barge_table = "\n".join(
    "| {} | {} | {} | {} | {} | {} |".format(
        f"{r['heel_deg']:.0f}", fmt(r["kn_m"]), fmt(r["gz_m"]), fmt(r["wall_sided_gz_m"]),
        "—" if r["difference_m"] is None else f"{r['difference_m']:.2e}", fmt(r["waterline_z_m"]),
    )
    for r in rows
)
round_table = "\n".join(
    "| {} | {} | {} | {} | {} |".format(
        f"{r['heel_deg']:.0f}", fmt(r["gz_m"]), fmt(r["analytic_gz_m"]),
        fmt(r["waterline_z_m"]), fmt(r["analytic_waterline_z_m"]),
    )
    for r in round_rows
)
imo_table = "\n".join(
    "| `{}` | {} | {} {} | {} |".format(
        c["key"], c["requirement"], fmt(c["measured"], 4), c["unit"], "yes" if c["satisfied"] else "NO"
    )
    for c in imo["checks"]
)

markdown = f"""# 大倾角稳性：从倾斜剖切得到 GZ 曲线（TensorCAD 0.21.0）

## ① 案例介绍

`tensorcad.stability` 只到正浮状态：给出 `GM`，即复原力臂曲线在零倾角的斜率。
`docs/stability.md` 第 156 行把剩下的部分写清楚了——*"No large-angle stability. `GZ` needs a
heeled cut"*。本单就是那条剖切线：它需要两样包内原本没有的东西，**倾斜的水线面**和**下沉量求解**
（排水体积随倾角保持恒定，水线必须随之上下移动）。

物理背景：横倾角 φ 下的复原力臂 `GZ(φ) = KN(φ) − KG·sin(φ)`，其中 `KN` 是浮心竖直线到中纵剖面的
横向距离（只与船型有关），`KG` 是重心高度（装载输入，自龙骨量起），`sin` 项是 G 自身随船体倾斜
外摆的横向偏移。正倾角 = 右舷下沉，正 `GZ` = 复原。

| 算例 | 几何 | 参数 | 用途 |
|---|---|---|---|
| 直舷驳船 | L=30 m, B=6 m, T=2 m, 型深 10 m | KG=1.0 m | 精确基准（直舷公式） |
| 圆舭 | R=3 m, L=20 m, 360 面片 | KG=1.5 m | 独立解析路径（弓形几何） |
| 实船型 | `generate_hull(12.0, 2.4, 0.9)` | KG=0.40 m | 跨模块一致性（与 `initial_stability`） |

## ② 怎么做的

方法链（全部在 `tensorcad/gz.py`，纯标准库）：

1. `heel_mesh(mesh, φ)`：绕通过龙骨中线的纵轴旋转（正倾角 = 右舷下沉）。
2. `submerged_solid(mesh, z_w)`：把曲面**裁剪**（不是采样）到半空间 `z ≤ z_w`——整块在下者保留、
   整块在上者丢弃、跨越者用 Sutherland–Hodgman 裁剪后扇形三角化。
3. `solve_waterline(mesh, ∇₀)`：以水线面积为 `dV/dz` 做牛顿迭代求解 `V(z) − ∇₀ = 0`，用二分兜底。
4. `gz_curve(...)`：逐角度求 `KN = y_B`，组装 `GZ`；`gm_upright` 由裁剪体积与水线面二阶矩自算。
5. `imo_intact_stability(curve)`：2008 IS Code Part A 2.2 六条数值衡准的筛查。

**关键实现决策（踩坑后改的）**：水线面盖板既不做三角化，也不把剖切线段串成环。每条被剖切三角形
的外法向给出剖线方向 `ẑ × n`（俯视时材料在逆时针一侧），于是有向线段的带符号求和就是扣掉孔洞后的
面积，盖板对体积一阶矩的贡献是闭式的（六体积和加 `2·wl·A`，一阶矩加 `6·wl·x̄·A`、`6·wl·ȳ·A`、
`6·wl²·A`）。最初照搬 `stability.py` 的"按网格自身边身份串环"做法，在圆舭网格上**面积正好算成一半**：
两片以很小夹角相接的面把同一个交点键成了不同的键，环断裂成三段。串环法不适用，改为有向线段后彻底
消除该失效模式；若水线恰好沿网格边走过（无剖面），闭合性检查（每个交点被两段共享）会抛 `GZError`，
而不是报一个"半份面积"。

工具链：Python 3.11 标准库 + matplotlib（仅出图）。重算量极小（最大网格 1800 面，12 个倾角 × 牛顿
5~6 次迭代），在 load15≈126 的机器上单次运行 < 3 s，符合本单"高负载只做轻量活"的容量纪律。

## ③ 计算结果

直舷驳船（∇ = {fmt(barge['displacement_volume_m3'], 3)} m³，全程恒定；GM = {fmt(barge['gm_upright_m'])} m）：

| 倾角[°] | KN [m] | GZ [m] | 直舷闭式 [m] | 差 [m] | 水线 z [m] |
|---:|---:|---:|---:|---:|---:|
{barge_table}

圆舭（360 面片，排水量取多边形面积之半）：

| 倾角[°] | GZ [m] | 解析 (R−KG)sinφ [m] | 水线 z [m] | 解析 R·cosφ [m] |
|---:|---:|---:|---:|---:|
{round_table}

IMO 筛查（直舷驳船，KG=1.0 m）：

| 衡准 | 要求 | 实测 | 通过 |
|---|---|---:|:--:|
{imo_table}

![gz](figures/gz-curve.png)

左：实测 `GZ`/`KN` 曲线与闭式解叠合（直舷段、圆舭段）；右：驳船 30° 横倾横剖面、水线与浮心 B（数值
标注为程序注入）。

## ④ 注意点（踩过的坑 / 局限）

* **串环法在圆舭上给出恰好一半的面积**（见②）——这是本单最大的坑；「精确」的前提是积分方案本身
  不依赖网格拓扑的巧合。
* **水线沿网格边**（例如 360 顶点圆上恰好有顶点落在半吃水处）时无剖面可言，必须报错，不能给出数。
* **恒吃水（不自由纵倾）**：水线保持平行于静水线，得到的是经典横截曲线；会随倾角纵倾的船型会有差异，
  差异在纵倾角的二阶。
* **`generate_hull`/`generate_wigley` 把模型封在设计水线，即零干舷**：一倾斜甲板边缘就入水，故其 `GZ`
  只是"不进水"理想化，示例只报到 20°；大倾角必须用带干舷的模型（本单基准用驳船），"带干舷船型族"
  是下一个欠账。
* 无自由液面修正、无谷物/横倾力矩衡准、无破舱、无船级社认可——IMO 部分只是**筛查**而非稳性手册。
* `compute_hydrostatics` 对**整个闭合网格**积分（包内约定网格即湿表面），所以对带型深的驳船它报的是
  1800 m³ 甲板以下全体积；水线以下的排水体积由 `submerged_solid` 负责。这一差异在示例里踩到并记入注释。

## ⑤ 误差判定（逐项）

| 关键量 | 基准 | 实测偏差 | 是否 < 3% |
|---|---|---:|---|
| 直舷驳船 `GZ`（7 个角度） | 闭式 `sinφ(GM + ½BM tan²φ)` | **{wall['measured']:.2e} m**（绝对） | ✅ 是（机器精度，远优于 3%） |
| 圆舭 `GZ`（12 个角度） | 弓形几何解析 `(R−KG)sinφ` | **{round_claim['measured']:.2e} m** = {100*round_claim['measured']/max(abs(r['analytic_gz_m']) for r in round_rows):.3f}% | ✅ 是 |
| 圆舭水线高度 | 解析 `R cosφ`（面积平分弦必过圆心） | **{float_claim['measured']:.2e} m** | ✅ 是 |
| 实船型 `GM` | `initial_stability.gm_transverse` | **{gm_claim['measured']:.2e} m**（逐位相同） | ✅ 是 |
| 排水体积守恒（所有角度） | 正浮排水体积 | < 1e-9 相对 | ✅ 是 |

结论：**直舷驳船一项在机器精度（{wall['measured']:.1e} m）上与教科书闭式解一致**，圆舭残差
{round_claim['measured']:.1e} m 是 360 面片的离散化而非剖切误差，全部关键量均远低于 3% 门。
honest 声明：本单没有实验/船模数据可比，因此判据用的是**解析闭式**而非实验值；IMO 部分未经船级社复核，
仅为数值衡准筛查。

## ⑥ 局限与待办

自由纵倾求解、自由液面修正、带干舷的船型族、破舱稳性仍缺；`docs/gz.md` 与 CHANGELOG 0.21.0 已记录。
"""

figure_bytes = (out / "figures" / "gz-curve.png").read_bytes()
print("figure bytes", len(figure_bytes))

payload = {
    "title": "大倾角稳性 GZ 曲线：倾斜剖切+下沉量求解，直舷驳船对闭式解 3.3e-16 m（TensorCAD 0.21.0）",
    "summary": (
        "tensorcad.gz 补齐 stability.py 明确留下的欠账：倾斜水线面剖切 + 排水量守恒下沉量求解，"
        "给出 GZ(φ)=KN(φ)−KG·sinφ 与 IMO 2008 IS Code Part A 2.2 六条衡准筛查。"
        "直舷驳船对教科书闭式解 3.3e-16 m（机器精度），360 面片圆舭对弓形解析 9.9e-5 m，"
        "实船型 GM 与 initial_stability 逐位相同。踩坑：串环法在圆舭上把面积算成恰好一半，"
        "改为外法向定序的有向线段求和后消除。31 个测试，commit b493dee。"
    ),
    "scenario_key": "tensorcad/gz-large-angle-stability",
    "mission_id": 3479,
    "tags": ["TensorCAD", "stability", "GZ", "verification-debt", "IMO", "large-angle"],
    "markdown": markdown,
    "figures": [
        {
            "name": "gz-curve.png",
            "caption": "左：实测 GZ/KN 曲线与直舷闭式解、圆舭解析解叠合；右：驳船 30° 横倾剖面与水线、浮心 B。",
            "png_base64": base64.b64encode(figure_bytes).decode(),
        }
    ],
}

base = os.environ["TENSOR_API_BASE"]
login = u.Request(
    base + "/auth/login",
    data=json.dumps({
        "username": os.environ["TENSOR_AGENT_USER"],
        "password": os.environ["TENSOR_AGENT_PASS"],
    }).encode(),
    headers={"Content-Type": "application/json"},
)
token = json.load(u.urlopen(login, timeout=30))
token = token.get("access_token") or token.get("token")
request = u.Request(
    base + "/tensor/reports",
    data=json.dumps(payload).encode(),
    headers={"Content-Type": "application/json", "Authorization": "Bearer " + token},
)
response = json.load(u.urlopen(request, timeout=90))
print("report response:", json.dumps(response)[:300])