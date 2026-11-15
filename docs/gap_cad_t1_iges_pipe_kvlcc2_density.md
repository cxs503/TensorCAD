# KVLCC2 出版方 IGES：采样密度跟随实测单位因子（锚点 GAP-CAD-T1-IGES-PIPE，m21995）

## 1. 案例
输入 = KRISO/SIMMAN 出版方原样交付的 `KVLCC2/kvlcc2.igs`（6 270 129 B，sha256 `ed06c2…e2f9`），
D 段 `{102:38, 126:304, 128:19, 142:19, 144:19}` —— 19 张裁剪面（144），全部基曲面被裁剪，船体只存在于 trim 环内。
出版基准（SIMMAN 2020）：L_pp 320.0 m、B 58.0 m、T 20.8 m、∇ 312 622 m³、C_B 0.8098；判据 = 锚点验收 3 %。

## 2. 缺陷与根因
`fleet_baseline.build_ship` 用**声明的**单位读法算采样间距：
`spacing = (1/samples_per_metre) / (unit_scale(flag)·model_scale)`。KVLCC2 全局段标志 2（毫米）压着米制坐标；
单位门 `measured_unit_factor()` 已实测出矛盾并采纳米制读法，但间距仍按声明读法：0.125 m / 1e-3 = **125 文件单位**。
后果可量测：每张裁剪面只剩 `[4,4]` 网格、16 点中 4..9 点落在环内（`domain_area_fraction` 0.25–0.56），
半宽场只覆盖约四分之一船体，∇ = **10 875.6 m³**，误差 **−96.52 %**，判决 WITHHELD。
（`docs/baseline_fleet_trim_intake.md` 记的 −0.13 % 在本树不可复现；同一输入同一判据实测 −96.52 %，本单以实测为准。）
根因一句话：密度是物理量，必须除以**实际起作用**的单位因子；该因子由测量决定，不由文件头决定——门采纳了它，密度没跟随。

## 3. 修复
`tensorcad/fleet_baseline.py::build_ship`：单位门给出 `factor` 后，若与采样所用的 `metres_per_unit` 不一致，
按采纳因子重算 `spacing` 并重采一次（并记录 `unit_gate.density_resampled_for_adopted_factor`）。
新增开关 `density_follows_adopted_factor: bool = True`，`False` 精确复现修前行为 → before/after 是一条命令的实测而非断言。
米制档案（KCS 等）该分支不触发，行为不变。

## 4. 结果（`runs/kvlcc2_density/evidence.json`，程序化注入）
`python scripts/kvlcc2_density_evidence.py`

| 量 | 修前 | 修后 | 出版基准 | 修后误差 |
|---|---|---|---|---|
| 采样间距（文件单位） | 125.0 | **0.125** | — | — |
| 每面网格 | `[4,4]` | `[400,257] …` | — | — |
| 环内采样点 | 4 502 | **417 632** | — | — |
| 裁剪面采用 | 19/19 | 19/19 | 19 | — |
| ∇ | 10 875.6 m³ | **312 568.0 m³** | 312 622 m³ | **−0.017 %** |
| C_B | 0.0282 | **0.8097** | 0.8098 | **−0.017 %** |
| 水密 / 边流形 | ✓ / ✓ | ✓ / ✓ | — | — |
| 判决 | WITHHELD | **VERIFIED** | — | — |

dz 细化阶梯（0.1 → 0.05 m）：333 314.7 → 312 568.0 m³（−6.22 %），记录在 evidence 中不隐去。

## 5. 注意点
- 只在"声明单位被几何反驳"的档案上暴露（KVLCC2）；JBC 同类矛盾、KCS/ONRT 不触发，故此前未显形。
- `tests/test_iges_pipeline_gate.py::test_kvlcc2_fleet_baseline_status` 原断言修前行为（WITHHELD），
  且局部 `KVLCC2_SPEC` 把 T 抄成 KCS 的 10.8 m、C_B 抄成 0.8495；已改出版值并断言 VERIFIED。
- 本机 shell 120 s 墙钟上限：`test_iges_pipe_gate.py` 中走 `import_cad`（trim 三角化路线）的 KVLCC2 用例在此上限内
  跑不完（`docs/gap_cad_t1_iges_pipe_real_files.md` §5 已记该路线耗时问题），本单未据此判其通过或失败，记 not_run。

## 6. 误差判定
| 关键量 | 基准 | 实测 | 误差 | < 3 % ? |
|---|---|---|---|---|
| ∇ | 312 622 m³（SIMMAN 2020） | 312 568.0 m³ | −0.017 % | **是** |
| C_B | 0.8098 | 0.8097 | −0.017 % | **是** |
| L_pp（包围盒跨度） | 320.0 m | 322.0 m | +0.63 % | **是** |

## 7. 文件
`tensorcad/fleet_baseline.py` · `tests/test_iges_pipeline_gate.py` · `scripts/kvlcc2_density_evidence.py` ·
`runs/kvlcc2_density/{evidence.json,figure.png}`
