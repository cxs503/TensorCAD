# 几何模型库资产化：持久店 · 登记簿 · 增长指标 · 消费面

锚点 `GAP-GEOM-ASSET-LIB`（第 1 轮）。仓库 TensorCAD，分支 `agent/m20465`。

## ① 问题设定与基准出处

**问题。** 平台把 KCS/KVLCC2/JBC/ONRT/KP505 的 IGES 字节收进了持久店
`$TENSORCAD_ASSET_ROOT/geometry/hulls/<hull>/`（默认 `/root/tp-agent/assets`，每船一份
`MANIFEST.json`），但仓库里的消费者仍然指着**任务工作树**：

```
"iges_path": "/root/tp-agent/work/m18488/artifacts/baseline_fleet/raw/KCS/KCS_hull_SVA.igs"
"published_source": "artifacts/baseline_fleet/raw/pages/t2015_kcs.html"
"evidence.KCS.iges.path": "artifacts/baseline_fleet/raw/KCS/KCS_hull_SVA.igs"
```

工作树随任务关闭而删除，这些指针全部失效，而**失效是静默的**：`benchmarks/verified/baseline_fleet/result.json`
依旧写着 KCS `VERIFIED`，但没有任何一条数字还能被重跑出来。这就是本锚点要解的缺口：
字节有了持久家，**名字却没有**——没有可校验的登记簿、没有增长指标、没有消费面。

**基准（ref）出处。** 物理基准不是自造的，沿用该锚点既有的公开工作坊数据：

| 量 | 值 | 出处 |
|---|---|---|
| KCS 排水体积 ∇ | 52030 m³ | NMRI T2015 KCS 页面 Table（页面 SHA-256 `68d9c246…`，摘要串存于 `fleet_specs.json` evidence）|
| KCS Cb / LCB | 0.6505 / −1.48 %Lpp | 同上 |
| KVLCC2 ∇ | 312622 m³ | SIMMAN 2020 Ship Data Summary（页面 SHA-256 `c1a668f6…`）|
| JBC ∇ | 178369.9 m³ | NMRI T2015 JBC 页面（页面 SHA-256 `fc2ce905…`）|
| ONRT ∇ | 8500.78 m³ | NMRI T2015 ONRT 页面（页面 SHA-256 `5e134273…`）|
| 验收容差 | ∇ 相对误差 ≤ 1 % | `result.json.acceptance.volume_rel_error_percent_max` |

平台精度标准（>15 % 不可接受）之上，本锚点用更严的 1 % 与"是否 <3 %"两条线同时报。

## ② 方法与实现（格式/出处 + 代码指针）

### 2.1 四层结构

| 层 | 文件 | 作用 |
|---|---|---|
| 字节层 | `tensorcad/assetlib/store.py` | 解析店根（`$TENSORCAD_ASSET_ROOT`）、SHA-256、读每船 `MANIFEST.json` 并**对其声明做复核**、按内容哈希定位（改名不丢）、带期望哈希的原子导入 |
| 登记簿 | `tensorcad/assetlib/registry.py` | `assets/geometry_registry.json`：13 条资产行（asset_id/uri/sha256/bytes/hull/status/duplicate_of）+ 4 条船舶绑定 + 5 船块 + 快照；`verify()/resolve()/growth()/coverage` |
| 指针卫生 | `tensorcad/assetlib/specs.py` | 扫描/修复 `fleet_specs.json` 的 43 个指针，三态 `resolved / dead / withheld`，另 `cited`=已指向持久 URL |
| 消费面 | `tensorcad/assetlib/consumer.py` + `tensorcad/fleet_baseline.py::specs_from_json` | `resolve_pointer(value, sha256)`：`asset://` → 存在且哈希相符的路径 → **按内容哈希在店里找** → 原样返回。`geometry_path(ship)` / `open_geometry(ship)` 供下游直接用 |

### 2.2 关键格式选择

1. **指针写成 `asset://geometry/hulls/KCS/KCS_hull_SVA.igs`**，不是绝对路径：仓库要能跨机移植，
   店根由环境变量给。绝对路径只出现在运行工件里（证明这次跑的是哪块盘）。
2. **按内容修复（repair by content）**：`iges_path` 只在 `iges_sha256` 于店中命中时改写，
   绝不按文件名猜。8 条几何指针全部命中（KCS/KVLCC2/JBC/ONRT 各 2 条：`ships[*]` + `evidence.*`）。
3. **字节互相重复时保留规范副本**：`KCS_hull_SVA.igs` 与 `KCS_hull_SVA/KCS_hull_SVA.igs`
   同哈希，后者登记为 `duplicate`，`resolve()` 归一化到规范副本。
4. **字节完整性与物理验收分离**：`resolve()` 默认只要哈希相符就给路径；`require_verified=True`
   才要求店 MANIFEST 声明 `status=verified`。理由见第⑥节（KVLCC2 字节干净但排水量未修好）。
5. **重复实现规避**：复用店内既有的 `MANIFEST.json`（协调者资产）与
   `tensorcad/fleet/{manifest,sources}.py` 的 sha/来源范式，`fleet_baseline.specs_from_json`
   只加一个解析钩子（+ 一个可选 `asset_registry` 参数），不改几何/水动力链路。

### 2.3 页面与压缩包的“诚实修复”

31 条页面指针与 4 条压缩包指针**在店里没有字节**（发布方未授予再分发权，仓库不存字节）。
修复不是删掉它们，而是升级为**持久引用**：路径 → 发布方 URL（真值来源：
`tools/baseline_fleet/fetch_sources.py::PAGES/TARGETS`），同时保留 `sha256`（可核对是否同一页面）、
`pointer_status: {path: "cache-lost"}` 与原路径 `legacy_path`。于是"丢了"是一等状态，不是静默。

## ③ 参数表

| 参数 | 值 | 说明 |
|---|---|---|
| 店根 | `/root/tp-agent/assets` | `TENSORCAD_ASSET_ROOT` 可覆盖 |
| 布局 | `geometry/hulls/<HULL>/…` + `MANIFEST.json` | 5 船 |
| 登记资产 | 13 文件 / 20 165 115 B（哈希吻合 13/13：`ok` 2 + `present_unverified` 11）| 重复字节计一次 |
| 船体 | KCS, KVLCC2, JBC, ONRT, KP505 | 5 |
| 船舶绑定 | 4（fleet_specs 的 4 条） | coverage 1.00 |
| 网格分辨率 | dx = 0.5 / 0.25 / 0.125 m，dz = 0.1 m，ring = 240，采样 8 /m | `fleet_baseline.DEFAULT_*` |
| 验收容差 | ∇ 相对误差 ≤ 1 % | `result.json` |
| 基准状态 | 15.470 °C? 不适用（本锚点无粘流） | — |
| 环境 | python3.11, 本机 1 核串行；单分辨率 KCS 7.1–7.5 s | 计算量远低于集群门槛 |

## ④ 对照表（ref vs got + 误差，逐条判据）

### 4.1 指针修复（判据 ①：死指针全清）

| 指标 | 修复前（`git show HEAD:…fleet_specs.json`） | 修复后 | 判据 |
|---|---|---|---|
| 指针总数 | 43 | 43 | — |
| 解析成功 | 0 | 8（全部 `asset://`） | ✅ |
| 死指针 | 8 | **0** | ✅ |
| 无米之炊（页面/压缩包） | 35 | 0（→ 35 条 `cited` 持久 URL + `cache-lost` 标记） | ✅ |
| 几何指针可解析率 | 0.00 | **1.00** | ✅ |
| 船舶覆盖（绑定/船舶） | 0/4 | **4/4** | ✅ |

### 4.2 登记簿与字节一致性（判据 ②：登记簿可校验）

| 判据 | 结果 |
|---|---|
| 全部 13 条资产重新哈希 vs 登记值 | 13/13 哈希全部相符（`ok` 2 + `present_unverified` 11；后者仅表示店 MANIFEST 未声明物理已验证，非哈希问题）✅ |
| `fleet_specs` 记录 sha 与店字节逐位一致 | 4/4 船 ✅（KCS `6f5556b2…`、KVLCC2 `ed06c287…`、JBC `2de08d29…`、ONRT `ff061ab8…`）|
| 登记簿 vs 现店重建（防漂移） | asset_ids 与总字节完全一致 ✅ |
| 篡改检测（同长度改 1 字节） | `sha_mismatch` ✅（`tests/test_assetlib.py::test_tampered_byte_is_detected_with_the_reason`）|
| 长度变化 | `size_mismatch` ✅；文件删除 → `missing` ✅ |
| 无 sha 的死指针 | 拒绝改写，落 `unresolved` ✅ |

### 4.3 物理（判据 ③：店里的字节能用，不只是“在”）

消费者 `specs_from_json` 解析后，四条 `iges_path` 全部落在店内（`Path.is_relative_to(store.root)` 有测试守护）。

| 船 | 指针 | 读到 sha | ∇ ref [m³] | ∇ got [m³] (dx=0.5) | 误差 | Cb ref/got | LCB ref/got [%Lpp] | 判定 |
|---|---|---|---|---|---|---|---|---|
| KCS | `asset://…/KCS/KCS_hull_SVA.igs` | `6f5556b2…` ✅ | 52030 | 52214.08 | **+0.354 %** | 0.6505/0.65280 (+0.354 %) | −1.48/−1.506 | VERIFIED ✅（<1 %，<3 %）|
| KVLCC2 | `asset://…/KVLCC2/kvlcc2.igs` | `ed06c287…` ✅ | 312622 | — | — | — | — | WITHHELD（`ValueError: field produced fewer than two usable waterline slices`）❌（遗留缺陷）|
| JBC | `asset://…/JBC/jbc.igs` | `2de08d29…` ✅ | 178369.9 | 4.44 | −99.998 % | 0.858/0.00002 | — | WITHHELD（几何伸出 7.499 m vs Lpp 280 m，模型尺度未换算）❌ |
| ONRT | `asset://…/ONRT/…mergedR.igs` | `ff061ab8…` ✅ | 8500.78 | 0.01 | −100.0 % | 0.535/0.000 | — | WITHHELD（伸出 1.000 m vs Lpp 154 m，模型尺度 1:48.94）❌ |

结论：**1/4 船以店字节通过 1 % 判据（KCS +0.354 %）**；KVLCC2/JBC/ONRT 的失败原因与
`result.json`（m18488）记录**同类**——水位切片、模型尺度换算——属既有缺陷，不在本锚点修复范围，
但**现在它们是可复现地失败**，而不是"文件找不到"。

## ⑤ 收敛性 / 参数敏感性（两分辨率→三分辨率）

KCS 在店字节上做网格加密（同一命令，三个 dx）：

| dx [m] | ∇ got [m³] | 误差 vs 52030 | 网格 vs 场积分体积差 | Cb | LCB [%Lpp] | 用时 |
|---|---|---|---|---|---|---|
| 0.500 | 52214.08 | +0.354 % | −0.019 % | 0.65280 | −1.506 | 7.3 s |
| 0.250 | 52134.65 | +0.201 % | −0.021 % | 0.65181 | −1.502 | 7.1 s |
| 0.125 | 52075.28 | +0.087 % | −0.068 % | 0.65106 | −1.493 | 7.5 s |

误差随 dx 减半按 ≈0.57→0.43 倍衰减（近一阶），外推剩余误差 <0.05 %，同时 `mesh` 与
`field` 两条积分路径相差 <0.07 %（离散自洽）。即 **KCS 的 +0.354 % 不是巧合，是收敛到 ref 的序列**。

参数敏感性第二点（库增长侧，两组状态）：见 `artifacts/assetlib/growth_report.json` ——
可解析几何指针比例 0.00 → 1.00；店管字节 0 → 20 165 115 B；覆盖 0 → 1.00。

## ⑥ 失败模式与工件

**失败的（honest-FAIL，留第 2 轮）**

1. **KVLCC2 水位切片**：`field produced fewer than two usable waterline slices`（dx=0.5/0.25 均失败）。
   已试：更换解析来源（店字节，sha 一致）、两分辨率重跑——均同错。缺的方法：该 IGES 的
   水线附近面片可能只覆盖部分船艏艉，需要按 z 分层的切片回退策略（或改半边宽场重建的内插）。
2. **JBC / ONRT 尺度换算**：厂方字节是模型尺度（JBC 1:37.34 由 7.499 m→280 m 推得；ONRT 1:48.9355
   由 spec 记录），`fleet_specs` 的 `geometry_scale` 对 JBC 为 `null`（"normalised to published L_PP"
   未真正实现）。缺的方法：把 `build_ship` 的规范化步骤与 spec 的 `geometry_scale/scale_basis` 打通。
3. **页面/压缩包字节**：35 条指针只能升级为 URL + `cache-lost`，因为发布方未授权再分发、且本环境无出网。
4. **归档漂移**：`result.json`（m18488）里 KCS 记 +0.304 %，本次同默认参数得 +0.354 %；JBC 记
   145769 m³（−18.28 %），本次 4.44 m³。二者之间 `tensorcad/iges.py` 被 m20163 改过（`git log 35f9cf5..HEAD -- tensorcad/iges.py`）。
   即**归档基准已不能复现**，本次数字是当前代码的真值。建议单开一单：把归档 result.json 重新生成并加
   "代码版本→数字"绑定（本次已把 sha 与读数写进工件，可直接复用）。

**工件（worktree 内）**

| 路径 | 内容 |
|---|---|
| `assets/geometry_registry.json` | 登记簿（13 资产 + 4 绑定 + 5 船块 + 快照） |
| `assets/geometry_registry_history.jsonl` | 增长快照历史（每行一次构建） |
| `benchmarks/verified/baseline_fleet/fleet_specs.json` | 43 条指针全修：8 `asset://` + 35 URL |
| `artifacts/assetlib/pointer_scan_before.json` / `_after.json` | 修复前后指针三态逐条 |
| `artifacts/assetlib/pointer_repair.json` | 43 条改动日志（每条含 sha 与新旧值） |
| `artifacts/assetlib/registry_build.json` | 登记簿构建快照 + 问题清单（空） |
| `artifacts/assetlib/physics_from_store.json` | 每船：解析路径、实读 sha、两/三分辨率 ∇/Cb/LCB 与误差 |
| `artifacts/assetlib/growth_report.json` | 增长指标（before/after/delta + 快照历史） |
| `artifacts/assetlib/assetlib_growth.png` | 图：(a) 死指针清零 (b) 各船店管字节 (c) KCS 误差-分辨率收敛 |
| `tensorcad/assetlib/*.py`, `tools/assetlib/*.py`, `tests/test_assetlib.py` | 实现、运维面、14 条测试 |

**repro_cmd（一键复现）**

```bash
export TENSORCAD_ASSET_ROOT=/root/tp-agent/assets
python -m tensorcad.assetlib build-registry            # 重建登记簿（失败退出码 1 并列出问题）
python -m tensorcad.assetlib verify                    # 重新哈希 13 条资产
python -m tensorcad.assetlib scan-specs                # 指针三态（期望 dead=0 withheld=0）
python -m tensorcad.assetlib resolve KCS               # 船名 → 店内路径（哈希校验）
python tools/assetlib/verify_from_store.py --ships KCS --dx 0.5 0.25 0.125
python tools/assetlib/growth_report.py
python -m pytest tests/test_assetlib.py -q             # 14 passed
```

**误差判定（显式回答）**：KCS 排水体积 +0.354 %（dx=0.5）/ +0.201 %（0.25）/ +0.087 %（0.125），
**<3 % 成立、<1 % 成立**；Cb 同误差；LCB 绝对偏差 0.026 %Lpp（ref −1.48）。KVLCC2/JBC/ONRT
**无有效 got 值，不能声称 <3 %**，如实记 WITHHELD，原因见上。

---

## ⑦ 第 3 轮补正（m23063）：按名字可达性 · hull-primary · 诚实的 hull_coverage

### 7.1 缺口是什么（本轮实测，不是推测）

第 1 轮把「字节有持久家、名字有可校验登记簿」做完了，登记簿里 **6 个船块**（DTC/JBC/KCS/
KP505/KVLCC2/ONRT）齐备。但第 1 轮的 `resolve()` 只会走 **船舶绑定**（`bindings`），而
`bindings` 只由 `fleet_specs.json` 的 4 条 ships 生成。于是：

```
$ python -m tensorcad.assetlib resolve DTC
{"name": "DTC", "resolved": false, "reason": "not_registered",
 "detail": "ship 'DTC' has no geometry binding in the registry"}
```

`DTC` 与 `KP505` 的字节明明在店里、DTC 的 `MANIFEST.json` 甚至声明 `status=verified` 且有
m21785 的实测量记录，**消费面却取不到**。更糟的是快照指标把这件事盖住了：

| 指标（第 1 轮） | 值 | 为什么掩盖了缺口 |
|---|---|---|
| `coverage` | **1.0** | 分母是「被 fleet_specs 绑定的船」=4，只统计 4 条，天然满分 |
| 实际可达船数 | **4 / 6** | DTC、KP505 存在但按名字不可达 |

即「库存完备感」是分母造成的假象——这正是本库存在的理由（静默失效要变成可检出失效）。

### 7.2 修法：两条规则 + 拒绝猜测

`AssetRegistry.hull_primary(hull)`（`tensorcad/assetlib/registry.py`）按固定优先级取
**船体几何**：

| 顺序 | `via` / `rule` | 判据 | 例 |
|---|---|---|---|
| 1 | `spec_binding` | fleet_specs 绑定的那条字节（实测基准指向谁就继续指向谁）| KCS→`KCS_hull_SVA.igs` |
| 2 | `hull_primary` / `unique_geometry` | 该船非重复 geometry 资产**恰好 1 条**，无可比 | KP505→`KP505_modelScale_250mm.igs` |
| 3 | `hull_primary` / `hull_form_name` | 唯一一条资产文件名主干**就是** `<hull>` 或 `<hull>_hull` | DTC→`DTC_hull.igs` |
| 4 | 拒绝 | 其余一切情况抛 `RegistryError` 并**列出全部候选** | 两条 `omega_*.igs` → 报错不选 |

第 4 条是关键：`DTC_bilgekeel.igs`、`DTC_rudder_skeg.igs` 的主干**也**以 `DTC` 开头，前缀
匹配会把艏/舭龙骨当船体发出去。因此只接受精确主干集 `{hull, hull_hull}`；不确定就报错，
把候选摊开给调用者——本库的全部价值就在「不按文件名猜」。

`resolve_explained(name)` 返回 `{asset_id, path, sha256, bytes, verified, via, rule}`，
让调用者看得见这次解析是**实测绑定**还是**规则推断**（`resolve()` 行为不变，返回 path，
内部改为委托）。

### 7.3 指标补正：`hull_coverage`

`snapshot_dict()` 新增 `hulls` / `hulls_resolvable` / `hull_coverage`（分母=登记簿船块数，
分子=按上述规则真的能取到字节的船数），`growth()` 的 `delta` 同步带 `hull_coverage`。
旧 `coverage` 保留（语义未变：fleet-spec 船的解析率），两个数并列就不再有「1.0 假象」。

现状（本轮重建后写入 `assets/geometry_registry.json` + 一行 history）：

| 指标 | 值 |
|---|---|
| `spec_ships` / `spec_bindings` / `coverage` | 4 / 4 / **1.0** |
| `hulls` / `hulls_resolvable` / `hull_coverage` | 6 / 6 / **1.0** |

### 7.4 验证（判据 ④：可达 = 哈希一致 = 就是被实测的那块字节）

```
$ python -m tensorcad.assetlib resolve DTC
{"name": "DTC", "resolved": true, "asset_id": "geometry/hulls/DTC/DTC_hull.igs",
 "path": ".../assets/geometry/hulls/DTC/DTC_hull.igs",
 "sha256": "53c3be90385c7abae6590276a7225edd1b78edb0e5e4c4a446428296e0edc274",
 "verified": true, "via": "hull_primary", "rule": "hull_form_name"}
```

该 sha256 与实测量记录 `benchmarks/verified/fleet_refs/dtc_hull_measured.json`
（`source.hull_sha256`）**逐位相同**——即解析到的就是被测量过的那块字节，
不是「同名但另一份文件」。对应实测量：∇ 173458.99 m³ vs 173467.0 m³ = **−0.0046 %**，
Cb −0.039 %，B +0.004 %。

| 判据 | 结果 |
|---|---|
| DTC/KP505 按名字可达 | ✅（`geometry_path("DTC")` / `resolve KP505` 均返回店内路径） |
| 可达字节 == 登记簿 sha256 | ✅（`store.verify_asset` 每次解析都重算哈希） |
| 可达字节 == 实测量那块 | ✅ DTC `53c3be90…` 双向一致 |
| 歧义必须报错而非猜 | ✅ 合成 2 候选船 `OMEGA` → `RegistryError` 并列出 `omega_a.igs, omega_b.igs` |
| 物理误差 < 3 % | ✅ DTC ∇ **−0.0046 %**（既有记录，本轮只复核字节同一性；本轮**未新做** CFD/几何计算） |

**本轮不声称的部分**：KP505 仍是 `status=pending`（`verified:false`，`require_verified=True`
时正确抛 `unverified`），它可达 ≠ 它已验证；JBC/KVLCC2/ONRT 的 WITHHELD 结论不变。

### 7.5 新增测试与复现

`tests/test_assetlib_hull_primary.py`（13 条）：`hull_form_name` 选中船体而非附体、
`unique_geometry` 覆盖单文件船、歧义拒绝并点名候选、未知船仍按名拒绝、绑定优先级高于
一切命名规则、`resolve_explained` 的 `via/rule/sha` 契约、`coverage` vs `hull_coverage`
的诚实性（4 船 1 歧义 → 0.75）、`growth` 带 `hull_coverage` delta；另有 4 条真店集成
（六船全部可达、DTC 解析到 `DTC_hull.igs` 而非 `DTC_bilgekeel.igs`、DTC sha 与实测量一致、
CLI 输出带 `via/rule`）。

可达性图与数据（程序化生成，零手抄）：
`python tools/assetlib/resolution_report.py --out artifacts/assetlib_hull_primary.png --json artifacts/assetlib_hull_primary.json`
→ 左图「6 船可达性 before 4/6 → after 6/6（标注每条用的规则）」，右图「DTC 被解析到的那块字节
的相对误差 vs 3 % 平台线 / 1 % 验收线」。数值：`hull_coverage_before = 0.666667`（新分母套旧规则）
→ `hull_coverage = 1.0`。

```bash
export TENSORCAD_ASSET_ROOT=/root/tp-agent/assets
python -m tensorcad.assetlib resolve DTC              # via=hull_primary rule=hull_form_name
python -m tensorcad.assetlib resolve KP505            # rule=unique_geometry, verified=false
python -m tensorcad.assetlib growth                   # delta.hull_coverage
python -m pytest tests/test_assetlib_hull_primary.py tests/test_assetlib.py -q   # 13 + 14 passed
python -m pytest tests/test_fleet_baseline.py tests/test_fleet_gate.py -q        # 36 passed
```

**踩过的坑**：① 前缀匹配会让 `DTC_*` 附体全部"像船体"，必须用精确主干集；
② `coverage` 这种「分母只含已绑定项」的指标会给出满分假象，新指标必须以**全部船块**
为分母；③ 修 `resolve()` 时若直接在里面加分支，会把「按 id 解析」和「按名推断」混成
一条路径，`via/rule` 就无处记录——所以拆出 `resolve_explained()`，`resolve()` 只做兼容委托。
