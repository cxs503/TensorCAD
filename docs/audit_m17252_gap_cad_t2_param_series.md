# 交差审计 m17252 — GAP-CAD-T2-PARAM-SERIES（TensorCAD Series 60 家族：方尾 + 球鼻艏）

审计人：board-audit（m17308）· 对象：`tensorcad/series60.py` @ TensorCAD `848e4b6`
（分支提交 6784903 / 1b94526 / 848e4b6）· 审计工作树：`/root/tp-agent/work/m17308`
被审工作树：`/root/tp-agent/work/m17252`（`git status` 干净，HEAD=848e4b6）

## 0. 结论

**级别 L2**（独立功能实现：真代码 输入→计算→输出；含真实算例、误差分析、可复现工件）。
**L3 未达成**：本单声明的接入口 `tensorcad/hull.py` 一行未改，`tensorcad/__init__.py` 未导出
`series60`，`tensorcad/cli.py` 只认 `generate_hull`，平台 `/generate` 仍走盒子路径——新内核
没有进入任何生产入口，只在仓内 examples/tools/geomquality 基线里被调用。
**L4 未达成**：唯一的"公开值对照"是**闭合校验**（公开表即输入），不是外部基准验证；
外部锚（Todd 原始型值表）为 WITHHELD（见 F1/F4）。
三件套（run_artifacts / ref-vs-got / repro_cmd）**真实存在且可复现**（本文 §2 亲测）。

## 1. 取证：三件套逐项核验

| 项 | 声明 | 实测 | 裁定 |
|---|---|---|---|
| run_artifacts | `runs/m17252_series60_family/*` | 5 个 json/csv + 3 张 png + gate_power.log，均为真实数值（非占位） | ✅ 真 |
| ref vs got | 五档 Cb × 五系数 vs `docs/data/series60_published_parents.json` | 逐行核对（本文 §2.1）：最大偏差 +0.38%（Cm @Cb=0.80），<3% 预置门 | ✅ 数值真，**但性质是闭合**（F1） |
| repro_cmd | `python3 tools/series60_family/audit.py --out runs/...` | 在隔离拷贝 `/tmp/aud` 重跑 25.7 s，与提交版工件 **460 个数值字段 0 差异** | ✅ 可执行、可复现 |

基准裁判：`SELECT case_name,metric,tolerance FROM tensor_benchmark_cases WHERE
anchor_slug='GAP-CAD-T2-PARAM-SERIES'` → **0 行**（全表 41 行/38 个 anchor 均无此 slug）。
即本锚点没有在台面上登记过判据与容差；3% 判据是写死在 `tools/series60_family/audit.py:63`
（`TOLERANCE = 0.03`）里的**代码内预置**，非台账预注册（F3）。

## 2. 我自己跑出来的数（零手抄，程序化注入）

### 2.1 复现与自洽

* `python3 -m pytest tests/test_series60.py tests/test_series60_fidelity.py
  tests/test_series60_bulb_transom.py -q` → **94 passed**（32.1 s，与 m17252 声明一致）。
* 重跑 `tools/series60_family/audit.py` → 与提交工件数值完全一致（460 字段，0 差异）；
  重跑输出：worst closed-form vs published `+0.38 %`、transom closure `2.2e-16`、
  bulb slab 差 `4.4e-16`、bulb mesh 体积差 `+0.96 %`。
* 门禁突变：提交件 `gate_power_series60.json` 记录 14 个突变体全杀；我抽验第一轮曾幸存的
  3 个（S60-9 / S60-12 / S60-13）→ **3/3 killed**，且各由具名测试杀死
  （`test_fully_submerged_volume`、`test_piercing_the_waterline_keeps_only_the_submerged_part`、
  `test_wetted_area_against_a_brute_force_summation` 等），"3 个测试洞已补"属实。

### 2.2 第三路径独立复核（本审计新增，作者没有这条路径）

作者的两条路径（闭式积分 / 二次参数化）共用同一套代数，无法发现"矩定律本身写错"。
我另写 `runs/m17308_audit/independent_recheck.py`：只从 `series60_section_law()` 取
解析出的 `n_aft/n_fore/cm_mid/gain/transom_ratio`，用 **scipy 自适应积分**（rtol 5e-12，
换元 `t=|2ξ−1|` 完全不用）重积分文档声明的三个矩 + 球鼻体积积分：

| 复核项 | 最坏误差 | 判据 | 结论 |
|---|---|---|---|
| 方尾扫掠 τ=0…0.30：重积分 Cwp/LCB/Cb vs 请求值 | 1.98e-14 | 1e-9 | ✅ 矩定律成立 |
| 五档 Cb：重积分 Cwp/LCB/Cb vs 公开请求值 | 1.69e-14 | 1e-9 | ✅ 成立 |
| 五档 Cb：重积分 VCB/T vs 工件列 | 6.0e-3 | 1e-2（网格 vs 闭式） | ✅ 一致（网格首阶亏损） |
| 球鼻 27 组：独立闭积分 vs 工件闭式体积 | 4.4e-16 | 1e-12 | ✅ 成立 |

工件：`runs/m17308_audit/independent_recheck.json`（含每行明细）。

## 3. 发现（按严重度）

### F1【中·声明问题，证据确凿】`role` 字段把"闭合"标成"独立"

`runs/m17252_series60_family/five_cb_gradient.json` 每行 `error_vs_published.*.role`
把 `waterplane_coefficient/lcb_over_length/block_coefficient` 标为
`input target (closure check)`（正确），却把 `prismatic_coefficient` / `midship_coefficient`
标为 `shape-derived (independent of the target)`。事实是：
`tensorcad/series60.py:42-51` 的 `_PARENT_FORMS` **硬编码了公开 Cp 阶梯**（0.625…0.805），
`Cm = Cb/Cp`、`Cp = Cb/Cm` 都是该阶梯的算术结果——**没有任何一个系数独立于公开表**。
所以 +0.38% 这个"对公开值的最大偏差"实际测量的是公开表自身的 `Cb = Cp·Cm` 互不自洽
（0.8/0.805=0.9938 vs 表里 0.990），不是实现精度。报告正文与 `docs/series60_bulb_transom.md`
§3.1 的文字是诚实的（明说"Cp/Cm 由公开 Cp 阶梯派生""不把闭合当独立验证"），
**但机器可读字段与文字互相矛盾**，任何按 JSON 自动判级的门读到的都是"独立验证通过"。
门槛证据：同一文件里唯一的真独立对照是**负对照**（未耦合家族 Cwp 偏 5–11%，能触发判据），
说明管线有咬合力，但它验证的是 `cwp` 目标耦合，不是 Cp/Cm。
**修复建议**：`role` 改为 `input ladder (closure check)`；报告 §6 表首行改为"闭合校验（非外部验证）"。

### F2【低-中·工件一致性】`gate_power.log` 与 `gate_power_series60.json` 结论相反

`gate_power.log`（mtime 03:13）写 `gate power: 11/14 killed`，幸存 S60-9/S60-12/S60-13
→ **78.6% < 80% 验收线**；`gate_power_series60.json`（03:17）为 14/14 全杀；
二者同在提交 `848e4b6` 里。只读 log 的读者会判本单**不满足**"gate-power ≥80%"验收。
我抽验 3 个原幸存者已全部被杀（§2.1），所以**能力为真、工件过期**：
log 是第一轮的中途快照，第二轮未覆盖。**修复**：重跑后覆盖 log，或在 log 头写"第一轮快照"。

### F3【低·判据登记】锚点无台面判据

`tensor_benchmark_cases` 无 `GAP-CAD-T2-PARAM-SERIES` 行；3%/2% 判据虽在代码内**先于**结果
声明（`audit.py:63-66`，符合预置纪律），但未进台账，因此门禁无法机器裁定容差。
**修复**：本单 verdict 已把两条容差写进 evidence，建议后续由台账补齐。

### F4【低·外部锚】external anchor = WITHHELD

`docs/data/series60_published_parents.json` 的 provenance 自陈：值来自仓内继承的文献链条
（`docs/series60_verification.md` §2 等），**无 URL/license/新抓取 sha**，只有该 JSON 自身的
sha256（只能防篡改，不能防转录错）。更关键的是**船型本身没有对外部型值表**：
`tests/test_series60_stations_published_comparison.py` 文件头即写"可与 Todd 1951 外部表比较"
——即比的是网格 vs 自身解析律。所以：系数闭合 ✅、形状对公开型值为 WITHHELD。
（不是造假：作者在 provenance 与 §5 里如实写了局限；这是覆盖边界而非诚信问题。）

### F5【提示·工具】`audit.py --out` 只接受 ROOT 内相对路径

`--out /abs/path` 会崩：`ValueError: '/tmp/aud_out/...' is not in the subpath of '/tmp/aud'`
（`audit.py:791` 的 `relative_to(ROOT)`）。我亲测复现；不影响仓内默认用法。

## 4. 六段报告核验（tensor_reports id=16260，4522 字符，3 图）

① 案例介绍（含几何/参数表/provenance）✅ ② 方法 ✅ ③ 数字表 ✅（与工件一致，我逐行核过）
④ 图 3 张（五档梯度/方尾/球鼻）✅ ⑤ 注意点（球鼻是独立实体需 boolean union、τ 是宽度比、
网格首阶亏损、不验证阻力/稳性）✅ ⑥ 误差判定（含"<3%"逐项回答 + 负对照）✅
**但**⑥ 的表首行把闭合当"published comparison"，与 F1 同源——文字层面已注明，不算虚标。

## 5. 未做 / 不能做的

* 未复核公开母型数值的文献原始出处（沙箱禁外网，无 Todd 原表可对）——F4 的 WITHHELD 即此。
* 未验证带球鼻/方尾船型的水动力（阻力/稳性）——本单也未声称，`docs §5` 已声明。
* 未跑全仓回归（205 passed 声明只抽验 3 个 series60 测试模块的 94 例 + 3 个突变体）。
* 未检查 ship-performance-platform 侧是否消费该几何（跨仓，超出本单）。

## 6. 处置

`request changes`（非阻塞返修）：F1 role 标签、F2 过期 log 两处属工件/声明缺陷，改一行即可；
能力本身（方尾 + 球鼻闭式、三路径自洽、14/14 突变、可复现）为真，**判 L2**，
不进 L3（未接 hull.py/API/平台）也不进 L4（外部基准为闭合）。