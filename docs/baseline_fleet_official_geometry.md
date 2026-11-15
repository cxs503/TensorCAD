# 基线船队官方几何路线：KCS 通过，JBC/KVLCC2 拒收并给出根因（GAP-CAD-T1-BASELINE-FLEET）

**锚点**：`GAP-CAD-T1-BASELINE-FLEET` — 6 船（KCS/KVLCC2/ONRT/JBC/DTC/Series60）从 NMRI T2015 /
SIMMAN2020 / Gothenburg2022 官方源取真实几何，转水密网格，逐船 provenance + 验证；验收＝
**6 船全部水密 + 流形；体积 vs 公开排水量误差 <1%；每船 manifest（源 URL/license/sha/验证结论）**。

**本单增量**：把"官方 IGES → 水密网格 → 判据"这条此前 `PENDING` 的管线**真正跑起来并落成可复算的判据机器**
（`tensorcad/fleet/official_geometry.py`），不是新增一份代理几何。结论：**6 船里 KCS 一条真通过
（体积 +0.348%，水密、流形、χ=2）；JBC/KVLCC2 被拒收，且拒收根因由数据给出，不再是一句话的举白旗**。
本单没有把锚点做完，这是有意为之的 honest FAIL：剩下三条船的阻塞点是具体的、可命名的、下一步可攻的。

---

## ① 案例介绍

| 船 | Lpp (m) | B (m) | T (m) | Cb | 公开排水体积 (m³) | 官方几何档案 |
|---|---:|---:|---:|---:|---:|---|
| KCS  | 230.0 | 32.2 | 10.8 | 0.6505 | 52030.0 | NMRI T2015 `KCS_hull_SVA.zip` → `KCS_hull_SVA.igs` |
| JBC  | 280.0 | 45.0 | 16.5 | 0.8580 | 178369.9 | NMRI T2015 `JBC_IGES.zip` → `jbc.igs` |
| ONRT | 154.0 (LWL) | 18.78 | 5.494 | 0.535 | 8500.1（页面 Cb 闭合；页面只给 8507 t 质量） | NMRI T2015 `onrt_wo_gap.zip` |
| KVLCC2 | 320.0 | 58.0 | 20.8 | 0.8098 | 312622.0 | SIMMAN 2014 `KVLCC2-Hull.zip` → `kvlcc2.igs` |
| DTC | 355.0 | 51.0 | 14.5 | 0.661 | 173467（闭合） | UD `dtc_igs_ordner.zip` → `DTC_hull.igs`（1345 patch） |
| Series60 | 100.0 | 13.333 | 5.333 | 0.60 | 闭合 Cb·L·B·T | **无官方 IGES（方法型系列，只有型值表）** |

物理背景：几何是整个船舶总体性能链的地基——阻力/耐波/稳性算例的 Cb/L/B/T 若不先验证，
下游 15% 的精度门根本无法判定是算法错还是几何错。所以本单把"几何先行验证（Cb/L/B）"做进判据。

## ② 怎么做的（方法 / 工具链）

新模块 `tensorcad/fleet/official_geometry.py`（`python -m tensorcad.fleet.official_geometry`）：

1. **按船登记 provenance**：源 URL、档案内成员名、publisher、license、发布页 URL + 页面 sha256、
   particulars 读法（`page-parsed` 还是 `literature`）、体积基准口径 `volume_basis`；
2. **字节级钉扎**：读盘后现算 sha256 与登记值比对，不一致直接判 `DIGEST_MISMATCH`（源变了不许静默通过）；
3. **按需解档**：`materialise()` 从本地档案缓存里按成员名解出 hull 文件；缓存缺失就如实报 `NO_OFFICIAL_GEOMETRY`，
   不下载、不编几何；
4. **建网格**：复用已有 `fleet_baseline.build_ship()`（IGES type-128 曲面 → 半宽场 → 水线下闭合实体），
   单位从 global 段读（不假定）；
5. **判据**（每条都带着它所依据的数）：水密、边流形、0 边界边、0 零面积面、χ=2（亏格 0）、
   |体积误差|≤1%、|Cb 误差|≤1%；trim 普查里只要有 patch 被 trim 边界绑定，就判 `TRIM_UNHANDLED`
   并说明"体积误差是漏 trim 造成的过覆盖，不是网格噪声"。

**关键纪律（踩坑后立的）**：JBC/KVLCC2 都声明 unit flag 2 = 毫米，但**抽样出来的坐标不是毫米**
——所以"读 global 段"还不够，必须把声明单位与几何自身尺度对账。这一步在本单被写成**判据**而不是注释。

## ③ 计算结果（程序化注入，未手抄）

**KCS —— VERIFIED**（原始字节 sha256 `6f5556b2…af17db917d`，与登记值一致）

| 量 | 计算值 | 公开值 | 相对误差 |
|---|---:|---:|---:|
| 排水体积 V (m³) | 52211.28 | 52030.0 | **+0.348%** |
| 方形系数 Cb | 0.65277 | 0.6505 | +0.348% |
| 中横剖面系数 Cm | 0.98541 | 0.9849 | +0.052% |
| LCB (%Lpp) | −1.509 | −1.48 | 0.029 %Lpp 绝对差 |
| 水线面系数 Cwp | 0.82174 | —（页面未给） | — |
| 船长（x 向 extent）(m) | 239.0 | 230.0 (Lpp) | +3.9%（含两舷外悬，非 Lpp） |
| 型宽 (m) / 吃水 (m) | 32.241 / 10.800 | 32.2 / 10.8 | +0.13% / 0.00% |

拓扑：三角形 52320、顶点 26162、**水密 = True、边流形 = True、边界边 0、零面积面 0、χ = 2**。
半宽场积分体积 52222.55 m³ 与网格体积差 −0.0216%（离散误差，可忽略）。用时 8.3 s（单核）。

**JBC —— WITHHELD（声明单位与坐标不符，实测）**
抽样（45 个 type-128 patch，全未 trim）：raw x ∈ [−3856.19, 3649.01]、y ∈ [−562.5, 562.5]、
z ∈ [0, 630.39]。声明 mm（×0.001）→ 船只有 7.5 m 长，与 280 m 差 37 倍。
按三轴比值反解：raw 比 7505:1125:630 = 6.67:1.0:0.56，公开 280:45:16.5 = 6.22:1.0:0.37；
用 B 定标得 **1 raw unit ≈ 0.04 m（＝1/25）**，则 L≈300 m（Lpp+外悬）、B=45.0 m、
z 向 25.2 m 正对上公开型深 25.0 m。**结论：文件声明 mm，实际是 25×m 的非标准尺度**，
按声明单位算出来的体积没有意义，故拒收而不是硬算一个数。

**KVLCC2 —— WITHHELD（两个独立缺陷，实测）**
`19/19` 个 type-128 曲面全部被 type-141 trim 边界绑定（304 条 type-126 曲线，`t126_trim_audit` 已定
trim 环在 BPTR 参数空间）。把 basis 面整体积分必然过覆盖，实测抽样：y ∈ [−8.16, 31.80]（**不对称**，
跨中过覆盖 31.8 m ≫ 半宽 29 m）、z ∈ [−29.18, 0.03]（含水上部分）、x 向 extent 333.6 m。
另外声明 mm、实际坐标是**米**（333.6 ≈ Lpp 320 + 外悬）。两条都要修，本单未修。

**字节级复核**：`KVLCC2-Hull.zip` 本单重新抓取，sha256 `dc032647f00a68d354ede2f4b1ee3d68cd048bc3d34d221c2a8c9aac3a6dd548`
与 `tensorcad/fleet/data/acquisition_plan.json` 里 first_retrieval 的 pin **逐字节相同**（ETag `"5cb08c74-148738"`、
Last-Modified `Fri, 12 Apr 2019 13:02:44 GMT`）。这是 acquisition plan 的第一次可复算验证。

## ④ 结果图

报告随附 `fleet_official_geometry.png`（matplotlib 生成，程序化注入）：左＝各船体积相对误差与 ±1% 验收带；
右＝各船 type-128 patch 总数与其中被 trim 边界绑定的数量（KVLCC2 19/19，JBC 0/45，KCS 0/55）。

## ⑤ 注意点（坑 / 局限 / 适用范围）

1. **声明单位不可信**：JBC（声明 mm、实为 1/25 m）与 KVLCC2（声明 mm、实为 m）各自把船做成 1/37 与
   1/1000。任何"读 global 段就好"的说法都被这两条否证。下一步必须做**单位对账**（用 Lpp/B/T 三轴共识
   反解尺度，且只允许离散倍率，不许连续拟合，否则体积门就变成自证）。
2. **trim 不是可选装饰**：KVLCC2 是 trim-heavy 档案，整面积分会产生**非对称**过覆盖（实测跨中 31.8 m），
   这类错误在自洽判据（水密/χ）里**完全看不出来**——正是"锚点验收必须对公开排水量"的原因。
3. **ONRT/DTC 本单未跑**：ONRT 已解出（132 patch、未 trim、含附体与甲板室，`z∈[0,draft]` 采样只取裸船体
   但附体会污染半宽场，需先剥离）；DTC 已解出 `DTC_hull.igs`（1345 patch、24 个 type-143 包装），
   单船构建超出本单剩余算力预算（本机 load15≈79，禁起重进程）。**这两条是 NOT_RUN，不是 PASS。**
4. **KCS 的 +3.9% 船长是特征不是误差**：文件含两舷外悬，x-extent 是总长量级而非 Lpp，勿据此判几何错。
5. **字节不随仓分发**：publisher 未声明 licence，仓里只放 sha256 与解档规则，跑之前需本地档案缓存
   （`TENSORCAD_HULL_CACHE`，默认 `/root/.cache/hull_sources`）。

## ⑥ 误差判定（关键量 vs 基准，逐项）

| 关键量 | 基准 | 误差 | <1%（锚点） | <3%（平台门） |
|---|---|---:|---|---|
| KCS 排水体积 | 52030 m³（NMRI T2015 页面，`page-parsed`） | **+0.348%** | **PASS** | PASS |
| KCS Cb | 0.6505（同页） | +0.348% | PASS | PASS |
| JBC 排水体积 | 178369.9 m³（同页） | **未评估（拒收：声明单位与坐标不符，实测 1/25 m）** | N/A | N/A |
| KVLCC2 排水体积 | 312622 m³（文献） | **未评估（拒收：19/19 trim 未处理 + 单位不符）** | N/A | N/A |
| ONRT / DTC / Series60 | — | **NOT_RUN** | N/A | N/A |

**关键量是否 <3%？** KCS 是（+0.348%，且同时落在 <1% 锚点验收内，体积与 Cb 两个独立量一致）。
JBC/KVLCC2 是**诚实拒收**：没有算出可信体积，就不给数字，也不拿代理几何的体积冒充官方几何的体积。
所以本锚点**仍未完成**（1/6 真通过）。

## ⑦ 下一步（可单干的次序）

1. **单位对账判据**：`official_geometry.reconcile_scale()`——用 Lpp/B/T 三轴比值找离散倍率（×1000、×100、
   ×25、×1…），把选中的倍率、三轴残差写进 manifest；对不上就拒收（JBC 会是 1/25，KVLCC2 会是 ×1）。
2. **trim 环裁剪**：用已有 `IgesFile.trimmed_surfaces()` + `Boundary` + BPTR 复合曲线，把每个 patch 的
   参数空间 trim 多边形建出来，做点在多边形内测试再采样；验收＝KVLCC2 体积进 312622 ±1%，且
   y 向抽样恢复对称。
3. **ONRT 附体剥离 / DTC 1345 patch 建网**（单船预算够时递交集群，别压本机）。
4. Series60 定性为"方法型系列无官方 IGES"，在 manifest 里写成 `NO_OFFICIAL_GEOMETRY` 而不是待办——
   锚点验收对它可以改成"型值表 vs Todd(1963) 数字化一致性"，需 owner 裁决口径。

## 复现

```bash
# 档案进本地缓存（本仓不分发 publisher 字节）
export TENSORCAD_HULL_CACHE=/root/.cache/hull_sources/extracted
python -m tensorcad.fleet.official_geometry --out artifacts/baseline_fleet/official --only KCS
# → artifacts/baseline_fleet/official/official_geometry.json
pytest tests/test_fleet_official_geometry.py -q
```