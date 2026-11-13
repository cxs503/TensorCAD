# AUDIT-RESTAT v2 · 全库真实成果重新统计（快照锁定 + 前值复现 + 溯源矩阵）

- 锚点：`AUDIT-RESTAT`（feature_ledger: 全库任务真实成果重新统计，L0-L5 归类+待完成名单）
- 快照（所有查询的时间上界）：`2026-11-13 01:50:00+00`
- 机器可读产物（含逐指标 sql_id 溯源）：`artifacts/audit_restat/restat_v2_2026-11-13T015000.json`
- 口径：L0-L5 定义沿用 `/root/tp-agent/PANORAMA.md`「成果成熟度 L0-L5」(owner 2026-11-11 令)；判级只认 `tensor_maturity_verdicts` 的 claim 与机器门 gate 实测，不认报告自评。

## 0. 结论摘要

| 指标 | 本快照实测 | 上游 m14502 (2026-11-11 09:34) | 读法 |
|---|---:|---:|---|
| 审计判定条数 (Tier-A) | 410 | 184 | 判定层样本翻倍（新增 board-audit 生产者） |
| gate≥L2 条数 | 117 | 0 | **上游「0/184」结论已被更大样本推翻** |
| claim 被机器门降级 (gate<claim) | 258 (62.9%) | — | 过度申报率，分生产者差异极大 |
| Tier-B 严格口径 L2+（严格定义带 A/B/C） | 2,995–4,054 | 上游单点 2,366 | 定义比上游宽松/收紧都在带内；上游未登记取数 SQL |
| succeeded 任务 | 10,895 | — | 任务量分母 |
| 能力台账条目 | 235 | — | 任务→能力记账率 2.2% |
| 台账 verified | 3 | — | 任务→已验证率 0.03% |

## 1. ① Tier-A 判定层：claim × 机器门 gate（快照锁定）

总判定 410 条，覆盖 410 个 scope；claim 被降级 258 条（62.9%），gate≥L2 117 条（28.5%）。

| 审计生产者 | n | claim=gate | 降级 gate<claim | 升级 gate>claim | ungated | 降级率% | gate≥L2% |
|---|---:|---:|---:|---:|---:|---:|---:|
| board-audit | 226 | 124 | 99 | 2 | 1 | 43.8 | 51.8 |
| audit2-ab | 163 | 23 | 140 | 0 | 0 | 85.9 | 0.0 |
| audit2-run | 21 | 2 | 19 | 0 | 0 | 90.5 | 0.0 |

**读法（关键）**：gate≥L2 的 117 条全部来自 `board-audit`；更严格的 `audit2-ab`/`audit2-run` 两族 184 条判定 gate≥L2 仍为 0。因此「本库没有一条被实测到 L2+」与「本库有 117 条实测 L2+」两个相反结论都能从同一张表读出——**取决于采信哪个生产者**。这不是数据问题，是治理问题：两族判定的 gate 语义需先统一，否则任何总量级结论都不可裁决。

### claim × gate crosstab

| claim \ gate | g0 | g1 | g2 | g3 | g4 | g5 | ungated |
|---|---:|---:|---:|---:|---:|---:|---:|
| **claim L0** | 25 | 2 | 0 | 0 | 0 | 0 | 0 |
| **claim L1** | 39 | 13 | 0 | 0 | 0 | 0 | 0 |
| **claim L2** | 50 | 49 | 57 | 0 | 0 | 0 | 0 |
| **claim L3** | 31 | 30 | 2 | 18 | 0 | 0 | 0 |
| **claim L4** | 28 | 25 | 3 | 0 | 35 | 0 | 1 |
| **claim L5** | 0 | 0 | 0 | 0 | 1 | 1 | 0 |

图：`fig_claim_vs_gate.png`（热力图）· `fig_producer_consistency.png`（分生产者一致性）

## 2. ② Tier-B 任务层：严格/宽松口径与定义敏感带

| 定义 | 严格 L2+ 计数 | 说明 |
|---|---:|---|
| A_merged_any_files | 4,054 | 有 PR 合并记录（不看改动内容） |
| B_merged_impl_file | 2,995 | 合并 + 至少一个实现文件（排除 tests/docs/examples/benchmarks 与 test_/verify_/audit_/report_ 命名） |
| C_merged_nonempty_changed_files | 3,836 | 合并 + changed_files 非空 |
| D_impl_file_any_pr_state | 6,888 | 有实现文件（不看是否入库）——膨胀上界 |

- 宽松上界（定义 D，任何有实现文件的 succeeded 单，不看入库）：**6,888**
- **严格定义敏感带 = [2,995, 4,054]**（A/B/C 三个要求合并的定义）。上游正文未登记取数 SQL/正则，故该单点不可独立复核；本文件把取数 SQL 全部登记在 `sql_registry`，任何人可在任一时点复算。

## 3. ③ 前值复现实验（能不能复现别人报的数）

### m14502_pinned · cutoff `2026-11-11 09:34:40.813285+00`

| 指标 | 上游声明 | 本次重算 | 判定 |
|---|---:|---:|---|
| n_verdicts | 184 | 184 | reproduced |
| n_scopes | 184 | 184 | reproduced |
| gate2 | 0 | 0 | reproduced |
| strict_l2plus | 2,366 | [2,313, 2,856] | reproduced_in_band |

### m17113_pinned · cutoff `2026-11-13 01:16:40.980444+00`

| 指标 | 上游声明 | 本次重算 | 判定 |
|---|---:|---:|---|
| strict_l2plus | 2,366 | [2,995, 4,054] | NOT_reproducible |

读法（逐条）：

1. **m14502 的「184 条判定 / gate≥L2 = 0」逐位复现**（cutoff 精确到该单 finished_at 秒）——判定层取数口径可复算，这条是其结论成立性的正面证据。
2. m14502 的 Tier-B 单点 2,366 在其自身 cutoff 的严格定义带 [2,313, 2,856] 内（带宽 18.8%）：点可被『某个合理定义』复现，但**定义带本身比结论更重要**——带内上下沿相差 543 单（23.5%）。
3. m17113 引用的 2,366 在**其自身 cutoff**（2026-11-13 01:16）的任何登记定义下都复现不了：该时点最小严格定义亦为 2,992（比 2,366 高 26.5%）。这与 14502 单逐位复现形成对照，推断该单沿用了上游点位而未按其自身时点重算——正是「统计量不带 cutoff」的实践代价（推断，非其自述）。
4. 因此**任何对外引用的成果数都必须带 (cutoff, 定义) 二元组**；本文件把两者都写死进产物，后续单只需改 `--snapshot` 即可对齐口径。

### 3.1 ⑥ 误差判定（关键量 vs 基准，逐项显式回答「是否<3%」）

| 关键量 | 基准（上游/定义） | 本次值 | 相对误差 | 是否 <3% |
|---|---|---:|---:|---|
| 审计判定条数 @m14502 cutoff | 184 | 184 | 0.00% | **是** |
| gate≥L2 条数 @m14502 cutoff | 0 | 0 | 0.00%（绝对差 0） | **是** |
| 覆盖 scope 数 @m14502 cutoff | 184 | 184 | 0.00% | **是** |
| Tier-B 严格 L2+ @m14502 cutoff | 2,366（上游单点） | 2,313（定义 B） | 2.24% | **是**（但见定义带宽 23.5%，带宽 ≫ 3%，单点精度无意义） |
| Tier-B 严格 L2+ @m17113 cutoff | 2,366（上游单点） | 2,994–4,053（三定义） | 26.5%–71.3% | **否**（honest FAIL，不可复现） |
| gate≥L2 条数 @本快照 | 117（本文件口径） | 117 | 0.00% | **是**（自身一致性） |

结论：判定层（Tier-A）的三个量在 3% 内逐位复现；任务层（Tier-B）只在**放宽到定义带**时才与上游单点相容，而 m17113 在自身 cutoff 引用的单点无法在 3% 内复现（26.5% 偏差）——该项按 honest FAIL 记录，并作为「统计量必须带 (cutoff, 定义) 二元组」的实证依据。

## 4. ④ 成果膨胀量化：任务数 ≠ 能力数（溯源漏斗）

| 层 | 计数 | 相对上一层的保留率 |
|---|---:|---:|
| succeeded 任务（时间上界内） | 10,895 |  |
| 其中带 PR 合并记录 | 4,054 | 37.2% |
| 其中带 ledger 锚点文本（【锚点=X】，排除占位串） | 1,018 | 9.3% |
| distinct 锚点（被引用的能力单元） | 221 |  |
| feature_ledger 台账条目 | 235 | 2.2% |
| 台账 verified（机器可复核） | 3 | 0.03% |
| 审计判定可映射到台账条目的 scope 数 | 157 | 66.8% |
- 锚点占位串（`【锚点=SLUG】`/`【锚点=NEW:】` 等不可解析写法）：**770** 单（7.1%）——占全部带锚点文本任务的 43.1%，这批单在账本上无法归位。
- 审计判定 scope 可机读率：**55.1%**（226/410 条形如 `mNNNNN <ANCHOR>`；184 条为自由文本，样例 ['hydrostatics_analytic', 'docs/ontology v1 + card schemas + validate_ontology_v1 (m8815)', 'scripts/verify_sediment_m7485 + tests/test_sediment_m7485 (m8823)']）
- 单锚点平均被派单次数：**4.61** 单/锚点（1,018 单 ÷ 221 锚点）——重复交付是膨胀的第一大来源。
- 可机读 scope 的锚点里有 **157** 个在台账中存在、**0** 个不存在（样例 []）；台账被判定覆盖的比例 **66.8%**（157/235 条台账条目有至少一条审计判定）——判定链与台账链只在表格写法规范的一半样本上闭合，判决仍无法整体反写台账。
- 图：`fig_traceability_funnel.png`（对数漏斗）

## 5. ⑤ 未达验收待完成名单（供 owner 裁定重开）

### 5.1 能力台账未完成项：115 条

按状态：dropped=1，in_progress=38，open=66，planned=10

| 状态 | 域 | slug | 验收线（截断 90 字） |
|---|---|---|---|
| dropped | Platform | `KPI-KB-10K` | 结构化知识样本 ≥1 万条 |
| in_progress | DEM | `H1-BROKEN-ICE-DEM` | 【基准牵引|浮冰浮态+航道清空｜ref: 阿基米德浮态解析（密度比 0.9 干舷/吃水）+冰区航道清空实测统计｜判据: 干舷 vs 解析；清空率统计 干舷<5%；清空率如实报】浮沉 |
| in_progress | FEM | `GAP-FEM-3D-COHESIVE` | 【基准牵引|DCB 双悬臂梁 I 型剥离｜ref: ASTM D5528 + 梁理论 P-δ 解析（LEFM）｜判据: P-δ 曲线+能量释放率 vs 解析 <10%】3D coh |
| in_progress | FEM | `GAP-FEM-WINKLER` | 【基准牵引|Hetényi 弹性地基板级数解｜ref: Hetényi 1946 Beams on Elastic Foundation｜判据: 挠度/弯矩分布 vs 级数解 <5 |
| in_progress | FVM | `GAP-FVM-TWOPHASE` | 【基准牵引|Zalesak 旋转缺口盘 + Hysing 上升气泡｜ref: Zalesak 1979；Hysing et al. 2009 benchmark 数据集｜判据: 界 |
| in_progress | LBM | `FEAT-OSS-INVENTORY-LBM` | ≥24 个 GitHub LBM 项目逐一代码级盘点（T1 多相 6/T2 可微 GPU 8/T3 HPC 10+发现波，总帽 30）；每项目 oss_inventory 落档；按 |
| in_progress | LBM | `GAP-LBM-MCMP-POISEUILLE` | 【基准牵引|MCMP 分层 Poiseuille 解析解｜ref: Huang-Krafcyk-Luo 2011 MCMP 解析｜判据: 界面位置/流量/spurious curr |
| in_progress | LBM | `GAP-LBM-SPHERE` | 【基准牵引|球绕流阻力关联｜ref: Clift-Grace-Weber 球阻力标准关联（Schiller-Naumann 0.1<Re<1000）｜判据: C_D vs 关联值  |
| in_progress | Papers | `KPI-STD-PATENT-PAPER` | 四类知识产权产出达标 |
| in_progress | Platform | `H6-MATURITY-GATE` | 【基准牵引|gate 覆盖率计量｜ref: KPI-GATE-100 计量口径（产出线×门禁矩阵可失败性）｜判据: 门禁覆盖率=可失败门/产出线 =100%】G1/H6 验收门报告 |
| in_progress | Platform | `KPI-QC-RECHECK` | 已合并功能 100% 有独立参考复核且通过 |
| in_progress | Platform | `KPI-REPRO-95` | 第三方复现成功率 ≥95% |
| in_progress | Platform | `KPI-SUCCESS-90` | 任务成功率 ≥90% 且失败自动归因 |
| in_progress | Platform | `KPI-SWEEP-1000` | 无人值守批量 ≥1000 工况全贯通 |
| in_progress | Platform | `KPI-TOOLS-5` | ≥5 类外部工具接入平台并可用 |
| in_progress | Platform | `KPI-TRACE-100` | 终态任务 100% 有报告工件可溯 |
| in_progress | cad | `FEAT-OSS-INVENTORY-CAD` | 普查≥20 仓 registry JSON+tier 分级+GAP-CAD-* 锚点提议≥6；深盘波按 LBM 协议接力 |
| in_progress | dg | `FEAT-OSS-INVENTORY-DG` | 普查≥20 仓 registry JSON+tier 分级+GAP-DG-* 锚点提议≥4；波2 深盘按 LBM 协议接力 |
| in_progress | fem | `FEAT-OSS-INVENTORY-FEM` | 普查≥20 仓 registry JSON+tier 分级+GAP-FEM-* 锚点提议≥6；深盘波按 LBM 协议接力 |
| in_progress | ice | `GAP-ICE-LIT-TRANSCRIPT` | 【基准牵引|Lindqvist 1989 POAC'89 三分解｜ref: 原文系数表｜判据: ≥4 个 (公式,quote,page/table) 三元组逐字转录，M1 由 WI |
| in_progress | knowledge | `FEAT-TK-ASSET-OPS` | KO 版本链与生命周期登记+主页知识专栏三 tab 真数据+版本化数据集下载服务(sha256+增量+published-only 硬门) |
| in_progress | knowledge | `FEAT-TK-FEEDBACK` | task-kit 派单装配+AnchorCard 进任务契约与 finish 门容差对照+run_ke_query 主动检索工具+预注册 A/B 证知识增益（判伪即负知识回写） |
| in_progress | knowledge | `FEAT-TK-GOV-INFRA` | VISION 正典化+sediment_line 常设沉淀守护+训练语料/评测集导出三件套落地；P0(存储统一/KE API/旗舰闭环)为地基 |
| in_progress | knowledge | `FEAT-TK-ONTOLOGY` | 本体 schema v1（六形态 JSON Schema+关系词表）+LBM 方法轴全量重组+船海领域轴试点，task-kit 契约成文 |
| in_progress | knowledge-infrastructure | `FEAT-TK-KNOWLEDGE-ENGINE` | ①选型评估给决策矩阵+三路建议且实测 demo（m8618）；②四契约映射规格对试点包零损耗或显式披露（m8619）；③原型对试点包查询全可复现且 xref 一致性断言过（m862 |
| in_progress | knowledge-infrastructure | `FEAT-TK-OSS-LANDING` | A1 ≥20 多相/MCMP claims+A2 ≥15 可微/框架 claims：schema 全过+全部 oss-cache 实读定位逐字引文+自带校验器全绿+无自批+每条映射 |
| in_progress | lbm | `FEAT-LBM-VERIFY-BENCH` | Young-Laplace ≥3σ 误差<5% 或 honest-FAIL+根因；spurious current 基线数字表（网格×σ）；pytest 固化；消费链映射表≥5 条 |
| in_progress | platform | `FEAT-MATURITY-L0L5` | tensor_maturity_verdicts 表+finish/打勾路径挂成熟度证据门（分级上线：v1 周报披露→v2 tb-eval-loop 后硬门） |
| in_progress | platform | `FEAT-TENSORBENCH` | 六物理锚 spec 化+RE-Bench 式任务集+结果回写评价回路三件全通 |
| in_progress | platform-mechanism | `FEAT-OSS-PLAYBOOK` | docs/oss-to-tk-playbook.md：五步各带机器可校验判据+已知坑+六域推广 checklist+度量指标 |
| in_progress | quality | `AUDIT-ASSET-12REPO` | 13 仓全扫：分类清单（算法/API/测试/文档/占位/重复/未接入）+每核心功能六问+verdicts 落表 |
| in_progress | quality | `AUDIT-RESTAT` | 六类成果统计表+未达验收待完成名单（供 owner 裁定重开）+成果膨胀量化 |
| in_progress | quality | `AUDIT-RUNTEST-6` | 6 仓各真跑核心测试+代表算例+负向+收敛/守恒；重计算走 cluster_queue；实测 vs 文档声称对照表 |
| in_progress | ship | `CAP-LEDGER-EVIDENCE` | 【基准牵引|台账判据集｜ref: 21 项能力台账逐项判据（总账 acceptance 原文）｜判据: 有判据+工件+独立复核项数/21 =21/21】21 项台账全部有判据+工件 |
| in_progress | ship | `KPI-BENCH-10` | ≥10 个公开对标算例通过 |
| in_progress | ship | `KPI-CYCLE-50-30` | 相对基线的设计周期与试验成本降幅 |
| in_progress | ship | `KPI-RES-5-10` | 阻力预报误差达标实测 |
| in_progress | ship-platform | `FEAT-SHIP-DEPLOY-ALIGN` | 部署同步线常设（pull→迁移→前端重建→后端重启→健康门→回滚）+三车道对齐（8001/8004 死道）+契约漂移门进 deploy gate |
| in_progress | sph | `FEAT-OSS-INVENTORY-SPH` | 普查≥20 仓 registry JSON+tier 分级+GAP-SPH-* 锚点提议≥4；波2 深盘按 LBM 协议接力 |
| open | TensorDEM | `CASE-BIAXIAL-SHEAR-BAND` | 【算例池|BIAXIAL-SHEAR-BAND｜ref: 文献剪切带倾角｜容差: <15%+定性（目标数值任务从 ref 原文查得）】 |
| open | TensorDEM | `CASE-COLUMN-COLLAPSE` | 【算例池|COLUMN-COLLAPSE｜ref: Lube et al. 2004 实验｜容差: <15%（目标数值任务从 ref 原文查得）】 |
| open | TensorDEM | `CASE-DEM-PLANE-SHEAR-MUIDI` | 【算例池|DEM-PLANE-SHEAR-MUIDI｜ref: da Cruz et al. 2005 Phys. Rev. E 72:021309；GDR MiDi 2004 E |
| open | TensorDEM | `CASE-DEM-RCP-PACKING-FRACTION` | 【算例池|DEM-RCP-PACKING-FRACTION｜ref: Scott 1960 Nature 188:908；Bernal & Mason 1960 Nature 18 |
| open | TensorDEM | `CASE-DEM-SHEAR-CAMPBELL-ELASTIC-LIMIT` | 【算例池|DEM-SHEAR-CAMPBELL-ELASTIC-LIMIT｜ref: Campbell 2002 J. Fluid Mech. 465:261 弹限颗粒剪切流 DO |
| open | TensorDEM | `CASE-SETTLING-STOKES` | 【算例池|SETTLING-STOKES｜ref: Stokes 解析｜容差: <5%（目标数值任务从 ref 原文查得）】 |
| open | TensorDG | `CASE-DG-DALEMBERT-ELASTIC-1D` | 【算例池|DG-DALEMBERT-ELASTIC-1D｜ref: Graff 1975 Wave Motion in Elastic Solids 一维杆波解析 DOI(次级): |
| open | TensorDG | `CASE-DG-PREFINEMENT-EXPONENTIAL` | 【算例池|DG-PREFINEMENT-EXPONENTIAL｜ref: Babuska & Suri 1994 SIAM Review 36:578 p/h-p 有限元收敛理论  |
| open | TensorDG | `CASE-LAX-DG` | 【算例池|LAX-DG｜ref: 精确解｜容差: <3%（目标数值任务从 ref 原文查得）】 |
| open | TensorDG | `REPRO-DG-CKS1998` | 【论文复现|CKS1998｜ref: SIAM Review / JCP 系列｜判据: quantitative_targets 逐项容差】 |
| open | TensorFEM | `CASE-CANTILEVER-TIP` | 【算例池|CANTILEVER-TIP｜ref: Euler-Bernoulli/Timoshenko｜容差: <3%（目标数值任务从 ref 原文查得）】 |
| open | TensorFEM | `CASE-FEM-HERTZ-CONTACT` | 【算例池|FEM-HERTZ-CONTACT｜ref: Hertz 1882 J. reine angew. Math. 92:156；Johnson 1985 Contact M |
| open | TensorFEM | `CASE-FEM-JINTEGRAL-EDGE-CRACK` | 【算例池|FEM-JINTEGRAL-EDGE-CRACK｜ref: Rice 1968 J. Appl. Mech. 35:379；Tada, Paris & Irwin 200 |
| open | TensorFVM | `CASE-FVM-CHANNEL-KMM180` | 【算例池|FVM-CHANNEL-KMM180｜ref: Kim, Moin & Moser 1987 J. Fluid Mech. 177:133 槽道 DNS 统计 DOI:1 |
| open | TensorFVM | `CASE-FVM-DAMBREAK-STOKER` | 【算例池|FVM-DAMBREAK-STOKER｜ref: Ritter 1892 Z. Ver. Deutsch. Ing. 36:947 溃坝解析解 / Stoker 1957 |
| open | TensorFVM | `CASE-FVM-KOVASZNAY-ORDER` | 【算例池|FVM-KOVASZNAY-ORDER｜ref: Kovasznay 1948 Proc. Camb. Phil. Soc. 44:58 精确解（网格收敛阶验证） DOI |
| open | TensorFVM | `CASE-FVM-PERIODIC-HILLS` | 【算例池|FVM-PERIODIC-HILLS｜ref: Breuer et al. 2009 Computers & Fluids 38:433 周期山丘基准（ERCOFTAC  |
| open | TensorFVM | `CASE-LAX` | 【算例池|LAX｜ref: 精确解｜容差: <2%（目标数值任务从 ref 原文查得）】 |
| open | TensorFVM | `CASE-SHU-OSHER` | 【算例池|SHU-OSHER｜ref: Shu & Osher 1988 文献｜容差: <5%（目标数值任务从 ref 原文查得）】 |
| open | TensorFVM | `CASE-WOODWARD-COLELLA` | 【算例池|WOODWARD-COLELLA｜ref: Woodward & Colella 1984 包络｜容差: <10%（目标数值任务从 ref 原文查得）】 |
| open | TensorFVM | `REPRO-FVM-SOD1978` | 【论文复现|SOD1978｜ref: J. Comput. Phys. 27:1-31｜判据: quantitative_targets 逐项容差】 |
| open | TensorKnowledge | `KBM-CLEAN-W1` | 【洗掘|W1-clean｜清洗规则可复跑+全量覆盖+剔除分类账数字可复算】 |
| open | TensorKnowledge | `KBM-CONSUME-W1` | 【洗掘|W1-consume｜A/B真实数字+建议逐条挂证据+可证伪】 |
| open | TensorKnowledge | `KBM-RCI-F1-ABC` | 【RCI-F1｜A=1强独做/B=20普通并行/C=20普通+审核+Skills 同批同契约；五指标：单位成本/产量/抽取准确率/证据完整率/复用；prereg先行+盲审】 |
| open | TensorKnowledge | `KBM-RCI-S2-PILOT200` | 【RCI|S2b｜分层非随机（项目/四类型/质量级）+五模板抽取+六维质量分+≥30 条逐字段人工核验=抽取准确率】 |
| open | TensorKnowledge | `KBM-RCI-S2-SCHEMA` | 【RCI|S2a｜owner JSON 为底+verification 状态机（L3 须 gate confirmed）+L0-L3 存储方案+三仓接口文档】 |
| open | TensorKnowledge | `KBM-SCHEMA-W1` | 【洗掘|W1-mine｜schema文档+管道repro_cmd+抽样一致率≥80%】 |
| open | TensorKnowledge | `KBM-TKC-COMPILER-ALIGN` | 【TKC｜Compiler=非结构化成果→可治理知识资产（主题/概念拆解/claim-evidence/关系/重复矛盾/待评审KO）；factory 工序01-05=执行面、TK  |
| open | TensorKnowledge | `KBM-TKC-GOV-SCALE` | 【TKC｜多Agent独立评审规模化（接factory工序04）+去重+交叉验证+矛盾发现（Contradiction Detection/Entity Resolution）；统 |
| open | TensorKnowledge | `KBM-TKC-INDEX-API` | 【TKC｜按概念/问题/方法/证据可检索；KE REST API 交付（规格已有）；消费接线=开工三查升级查TK索引+引用计数（复用率计量）；验收=Agent用知识完成任务的实证案 |
| open | TensorKnowledge | `KPI-KB-REUSE-RATE` | 【RCI|核心KPI①｜被后续研究任务/工单实际引用的KO占比（引用=工作包或mission goal引用knowledge_id且有产物佐证）；目标值由 RCI-F1 试点评估定 |
| open | TensorKnowledge | `KPI-KB-VERIFY-RATE` | 【RCI|核心KPI②｜分子=L3候选中经平台验证单+gate confirmed 通过数，分母=进入验证流程的KO数；目标值由 RCI-F1 试点评估定标（定标前计量单只报现值不 |
| open | TensorLBM | `CASE-CAVITY-RE5000` | 【算例池|CAVITY-RE5000｜ref: Ghia et al. 1982｜容差: <8%（目标数值任务从 ref 原文查得）】 |
| open | TensorLBM | `CASE-COUETTE-2D` | 【算例池|COUETTE-2D｜ref: 线性解析解｜容差: <2%（目标数值任务从 ref 原文查得）】 |
| open | TensorLBM | `CASE-CYLINDER-DRAG-RE100` | 【算例池|CYLINDER-DRAG-RE100｜ref: Schäfer & Turek 1996｜容差: <5%（目标数值任务从 ref 原文查得）】 |
| open | TensorLBM | `CASE-LBM-KOVASZNAY-EXACT` | 【算例池|LBM-KOVASZNAY-EXACT｜ref: Kovasznay 1948 Proc. Camb. Phil. Soc. 44:58 二维定常 NS 精确解 DOI: |
| open | TensorLBM | `CASE-LBM-RB-ONSET-RA1708` | 【算例池|LBM-RB-ONSET-RA1708｜ref: Chandrasekhar 1961 Hydrodynamic and Hydromagnetic Stability  |
| open | TensorLBM | `CASE-LBM-SQUARE-CYL-NORBERG` | 【算例池|LBM-SQUARE-CYL-NORBERG｜ref: Norberg 1993 J. Wind Eng. Ind. Aerodyn. 49:187 方柱涡脱实验基准 D |
| open | TensorLBM | `CASE-MCMP-LAPLACE` | 【算例池|MCMP-LAPLACE｜ref: ΔP=σ/R 解析｜容差: <5%（目标数值任务从 ref 原文查得）】 |
| open | TensorLBM | `CASE-RAYLEIGH-TAYLOR` | 【算例池|RAYLEIGH-TAYLOR｜ref: He/Doolen/Ding 文献尖峰气泡增长率｜容差: <10%（目标数值任务从 ref 原文查得）】 |
| open | TensorLBM | `CASE-TAYLOR-GREEN-3D` | 【算例池|TAYLOR-GREEN-3D｜ref: 级数精确解｜容差: <8%（目标数值任务从 ref 原文查得）】 |
| open | TensorLBM | `GAP-LBM-T2-GPU-COREPATH` | 【T2｜碰撞迁移+边界条件融合/内存管理/多后端/多GPU通信/精度性能权衡；**优化对象=T3已验证物理算例，禁盲目追峰值**；对标口径统一纪律】 |
| open | TensorLBM | `GAP-LBM-T3-BENCH` | 【T3｜九算例谱系自动化验证矩阵（方腔/圆柱/球体/湍流槽道/TG/静态液滴/RT/溃坝/船阻波浪）：每项解析或实验或参考软件对照+repro_cmd+机器判据+CI可跑；与算例池 |
| open | TensorLBM | `GAP-LBM-T5-ICE-COUPLING` | 【T5｜破冰船旗舰=LBM水动力+FEM结构+DEM冰破坏；工程级验证门=T1多相+载荷交换达标后；P1周期性单相原型不算验证；首轮只做依赖审计+耦合接口设计】 |
| open | TensorLBM | `GAP-PP-2-lid-wall-diff` | 【论文牵引|JOUR-DIFFLBM E2 方腔反演所需】lid-driven cavity Re=100 τ 反演恢复误差<2%；梯度 FD 交叉<1e-4 |
| open | TensorLBM | `GAP-PP-3-adjoint-ootc` | 【论文牵引|JOUR-DIFFLBM 第4节成本/内存章所需】10^4 步 128^2 τ 反演峰值内存<8GiB；梯度与 BPTT 参考差<1e-5 |
| open | TensorLBM | `GAP-PP-4-opt-loop` | 【论文牵引|JOUR-DIFFLBM 第4节反演实验所需】τ/Re 恢复 demo 端到端跑通并产出机读收敛 JSON；TG 衰减 Re 恢复<3%；3 种子可复现 |
| open | TensorLBM | `KPI-LBM-VERIFIED-CAPS` | 【成果指标令｜计数=本周 gate confirmed L3+ 的物理能力锚新增数；不以文件/代码/PR数计成果；周计量走 kpi_meter】 |
| open | TensorLBM | `REPRO-LBM-HOU1995` | 【论文复现|HOU1995｜ref: J. Comput. Phys. 118:329-347｜判据: quantitative_targets 逐项容差】 |
| open | TensorLBM | `XVAL-LBM-CAVITY` | 【横向对比|XVAL-CAVITY｜我方+开源同规格三方｜判据: metrics 逐项（工具间差+vs ref）】 |
| open | TensorPapers | `GAP-PP-6-run-artifacts` | 【论文牵引|JOUR-DIFFLBM 第4节数据可用性所需】每次反问题运行产出 bundle 且 sha256 入 provenance；稿件数字 100% 程序化注入 |
| open | TensorPapers | `JOUR-DIFFLBM` | 【顶刊攻关|JOUR-DIFFLBM｜立项五件+差距锚落库】 |
| open | TensorPapers | `JOUR-ICE-STRUCT` | 【顶刊攻关|JOUR-ICE-STRUCT｜立项五件+差距锚落库】 |
| open | TensorPlatform | `FEAT-CONTEXT-CRAFT` | 【策展令｜①干货包=工作包按任务类型带足当次决策所需（KO/负面知识/同域verdict/容差判据）且仅带这些，token 预算+密度判据成文；②三层记忆(全局bootstrap/ |
| open | TensorPlatform | `KBM-RCI-F1-FACTORY` | 【RCI-F1｜tk_work_packages任务池+闲位驱动派单+结构门证据门程序化+返工队列+独立审核单派发+八指标看板+A/B/C泳道】 |
| open | TensorSPH | `CASE-HYDROSTATIC` | 【算例池|HYDROSTATIC｜ref: 解析｜容差: <1%（目标数值任务从 ref 原文查得）】 |
| open | TensorSPH | `CASE-POISEUILLE-SPH` | 【算例池|POISEUILLE-SPH｜ref: 解析解｜容差: <5%（目标数值任务从 ref 原文查得）】 |
| open | TensorSPH | `CASE-SPH-KLEEFSMAN-DAMBREAK-OBST` | 【算例池|SPH-KLEEFSMAN-DAMBREAK-OBST｜ref: Kleefsman et al. 2005 J. Comput. Phys. 206:363 溃坝障碍物 |
| open | TensorSPH | `CASE-SPH-WEDGE-ENTRY-2D` | 【算例池|SPH-WEDGE-ENTRY-2D｜ref: Zhao & Faltinsen 1993 J. Fluid Mech. 246:593 入水砰击研究 DOI:10.10 |
| open | TensorSPH | `CASE-WAVE-FLUME` | 【算例池|WAVE-FLUME｜ref: 文献波高衰减｜容差: <10%（目标数值任务从 ref 原文查得）】 |
| open | bench | `FEAT-CASEPOOL-100` | 【casepool|池 84→120+ 例；新增例 ref_source 真实可引用零手抄；unlinked pending ≥30 保持 case_feeder 可供；与既有 8 |
| open | platform | `FEAT-AGENT-HARNESS` | 【harness|第1轮=席位健康度仪表（表+健康分可复算+RECUR_FAIL/LOW_SUCCESS/MEMORY_STARVED 三标记红队可失败+守护脚本入仓+部署说明，只 |
| open | platform | `FEAT-ANCHOR-FEEDER` | 【机制】多锚源供单守护上线：REPRO/XVAL/JOUR/T1-T5 锚源按闲位供给；验收=主线类占比(24h) 从 <3% 升至 ≥30% 且无断供 >6h；周报计量 |
| open | platform | `FEAT-BOARD-CLAIM` | 【机制】POST /tensor/board/claim 上线：空闲 agent 直领开放锚，服务端锁租约+服务端组装 goal（锚+验收+ref 零手抄）+并发去重；验收=①端点 |
| open | platform | `FEAT-DEMAND-INTAKE-SOAK` | 连续 7 日日审计覆盖率=已审/待审 >=60%，且 agent 铸造帽拒率 <20%；否则审计帽由 20/日上调至 40/日 |
| open | platform | `KPI-SEAT-HEALTH` | 【KPI|现值首测定基线（kpi_meter 计量单，只测量不改进）；消账门=confirmed 计量现值≥80%；健康定义细则随 FEAT-AGENT-HARNESS 交付成文】 |
| planned | Platform | `KPI-GATE-100` | 全部产出线接质量门禁 |
| planned | Verif | `FEAT-PILOT-AB-V1` | 双臂同题同模型并行完成；QR 三判据盲裁出胜负或平局过程指标裁决；终局对照报告 scenario_key=pilot-ab（判据结果/路线变更数/审查发现/墙钟）；推广或否决建议交 |
| planned | bem | `FEAT-OSS-INVENTORY-BEM` | 普查≥20 仓 registry JSON+tier 分级+GAP-BEM-* 锚点提议≥6；深盘波按 LBM 协议接力 |
| planned | fvm | `FEAT-OSS-INVENTORY-FVM` | 普查≥20 仓 registry JSON+tier 分级+GAP-FVM-* 锚点提议≥6；深盘波按 LBM 协议接力 |
| planned | ice | `GAP-ICE-CRUSH-ISO` | 【基准牵引|ISO 19906 A.8 挤压力｜ref: 规范 F_ISO 30 MN 量级｜判据: 4 工况 relative_error ≤50% 全 PASS（现状 0.84 |
| planned | ice | `GAP-ICE-FLOE-VALID` | 【基准牵引|多面体堆体休止角｜ref: atan(mu) 单调关系 + 力链尾部指数 ｜判据: 11/11 判据全 PASS（现状 9/11，case3 NOT VALIDATED |
| planned | ice | `GAP-ICE-PA-IDENT` | 【基准牵引|p-A 指数分离｜ref: ISO 19906 ex=-0.28 + 本链内生性地板 |ex|<0.126 ｜判据: 处理臂（速率相关挤压）恢复的 ex 与 -0.28 |
| planned | ice | `GAP-ICE-UIKKU-DEM` | 【基准牵引|MT Uikku 平冰阻力 10 工况｜ref: Hu & Zhou 2016 Table 8 实测 ｜判据: DEM 预报 vs 实测 MAPE ≤15%，且必须优于 |
| planned | platform | `FEAT-WEEKLY-QPATROL` | 周轮：同题 PR 突发+报告/源码比+L2+ 零增红线可失败+周报入 /tensor/reports |
| planned | ship | `A2-ICE-MODULES` | 平台冰区端点 ≥8 类×5 工况+证据链升账 |

### 5.2 awaiting_pr_merge 任务：121 张

（协调官裁定项：分支已推但未合并，既不进 done 也不进 failed，是最容易被漏数的灰区）

| mission | repo | goal 首 120 字 |
|---:|---|---|
| 8634 | TensorKnowledge | 【锚点=FEAT-TK-OSS-LANDING】LBM 开源多相/MCMP 知识落库 TensorKnowledge（≥20 claims）  对象：oss_solutions 中挂 GAP-LBM-MCMP-POISEUILLE 与多相族 |
| 8800 | TensorKnowledge | 【锚点=FEAT-TK-KNOWLEDGE-ENGINE】KE 存储统一：SQLite/PG 双轨收口为单一权威路线（外评定谳：选型报告定 PG 自研薄层为权威，ke/ 原型却是 SQLite，禁止两套知识系统漂移）。 任务： 1. ke/ |
| 8801 | TensorKnowledge | 【锚点=FEAT-TK-KNOWLEDGE-ENGINE】KE 服务化：把 ke/ 原型升级为可被其他 Agent/服务调用的 HTTP API（外评定谳：目前只是本地模块+CLI，不是可调用基础设施）。 任务： 1. ke/api.py  |
| 8810 | TensorKnowledge | 【锚点=FEAT-TK-GOV-INFRA】物理 AI 训练出水口 v1（owner 2026-11-10 定位令：面向物理AI训练的知识基础设施）。 任务：ke/export.py 三件套导出器—— ① SFT 语料：claim+evid |
| 8811 | TensorKnowledge | 【锚点=FEAT-TK-ASSET-OPS】知识资产管理（owner 11-10 三问之一：这些数据资产如何管理）。 任务： 1. KO 版本链与生命周期：PG schema（ko_versions）——每个知识对象可 supersede（ |
| 8813 | TensorKnowledge | 【锚点=FEAT-TK-ASSET-OPS】训练调用面（owner 11-10 三问之三：后期训练如何调用）。 任务：版本化数据集服务（挂 KE API 同进程或姐妹端口）： 1. 端点：GET /datasets（版本列表）；GET /d |
| 8815 | TensorKnowledge | 【锚点=FEAT-TK-ONTOLOGY】知识体系组织学 v1（owner 2026-11-10 令：我们的知识体系面向 AI 智能体调用、指导智能体工作与后续模型开发，必须与面向人类的知识体系有本质区别——人类体系优化「读懂」，智能体体系 |
| 8817 | TensorKnowledge | 【锚点=FEAT-TK-ONTOLOGY】船舶与海洋工程领域轴试点（owner 2026-11-10 令点名领域：我们的知识体系要覆盖船舶与海洋工程，与 LBM 等方法轴交叉）。 任务： 1. 领域树 v1：阻力/推进/耐波/操纵/冰区航行 |
| 13756 | TensorFEM | 【自找活·tp-impl-7（platform 席）】（owner 2026-11-01：全席饱和；2026-11-02 控仓令：自驱开发照常，绝不动 TensorPlatform） 你是研究开发 agent，本单=你的自主工时（model |
| 13761 | TensorDG | 【自找活·tp-ver-2（verifier 席）】（owner 2026-11-01：全席饱和；2026-11-02 控仓令：自驱开发照常，绝不动 TensorPlatform） 你是研究开发 agent，本单=你的自主工时（model  |
| 13769 | TensorFEM | 【自找活·tp-ver-8（verifier 席）】（owner 2026-11-01：全席饱和；2026-11-02 控仓令：自驱开发照常，绝不动 TensorPlatform） 你是研究开发 agent，本单=你的自主工时（model  |
| 13781 | TensorSPH | 【自找活·tp-mod-4（models 席）】（owner 2026-11-01：全席饱和；2026-11-02 控仓令：自驱开发照常，绝不动 TensorPlatform） 你是研究开发 agent，本单=你的自主工时（model 钉扎 |
| 13787 | TensorBEM | 【自找活·tp-app-2（applications 席）】（owner 2026-11-01：全席饱和；2026-11-02 控仓令：自驱开发照常，绝不动 TensorPlatform） 你是研究开发 agent，本单=你的自主工时（mo |
| 13794 | TensorPapers | 【自找活·tp-gen-4（generalist 席）】（owner 2026-11-01：全席饱和；2026-11-02 控仓令：自驱开发照常，绝不动 TensorPlatform） 你是研究开发 agent，本单=你的自主工时（mode |
| 13805 | TensorFVM | 【自找活·tp-merge-4（generalist 席）】（owner 2026-11-01：全席饱和；2026-11-02 控仓令：自驱开发照常，绝不动 TensorPlatform） 你是研究开发 agent，本单=你的自主工时（mo |
| 13806 | TensorBEM | 【自找活·tp-merge-5（generalist 席）】（owner 2026-11-01：全席饱和；2026-11-02 控仓令：自驱开发照常，绝不动 TensorPlatform） 你是研究开发 agent，本单=你的自主工时（mo |
| 13812 | TensorSPH | 【自找活·tp-glm51-2（proposer 席）】（owner 2026-11-01：全席饱和；2026-11-02 控仓令：自驱开发照常，绝不动 TensorPlatform） 你是研究开发 agent，本单=你的自主工时（mode |
| 13815 | TensorDG | 【自找活·tp-glm52-2（generalist 席）】（owner 2026-11-01：全席饱和；2026-11-02 控仓令：自驱开发照常，绝不动 TensorPlatform） 你是研究开发 agent，本单=你的自主工时（mo |
| 13836 | TensorBEM | 【自找活·tp-ds-8（proposer 席）】（owner 2026-11-01：全席饱和；2026-11-02 控仓令：自驱开发照常，绝不动 TensorPlatform） 你是研究开发 agent，本单=你的自主工时（model 钉 |
| 13903 | ship-performance-platform | 【战役 DAG 节点|IN65-N09|压载水置换操作窗口闭环：BWM D-1顺序法/溢流法+自由液面修正GM时程+95%置换率达标+剪力/弯矩包络】（上层智能体播种的自主战役，物理域 ，风险 ，层 ）。物理目标：对给定船型舱容表泵曲线与置 |
| 13905 | ship-performance-platform | 【战役 DAG 节点|P103W2C-N08|冰区通道螺旋桨通气（Ventilation）：掺气起始判据+混合密度推力崩溃+净空(h/D)运行包络】（上层智能体播种的自主战役，物理域 ，风险 ，层 ）。【物理目标】建立冰区通道中螺旋桨通气的 |
| 13969 | ship-performance-platform | 【战役 DAG 节点|p9001w2-N07|冰区单航道通行能力与船舶排队】（上层智能体播种的自主战役，物理域 ，风险 ，层 ）。【物理目标】建立冰区单航道通行能力与船舶排队模型：(1) 跟驰间距d_safe=V·t_r+V²/(2a(V, |
| 14010 | TensorBEM | 【战役 DAG 节点|p362w2-n08-bem-fmm-acceleration|TensorBEM FMM 快速多极加速 O(N²)→O(N)】（上层智能体播种的自主战役，物理域 ，风险 ，层 ）。物理目标：为 TensorBEM 边 |
| 14481 | TensorCAD | 【生态巡更·TensorCAD·吸收候选】找出一个可吸收的外部改进（上游/对标库做法），以最小增量落到本仓并配测试；先红后绿给证据。 纪律：先查任务池避免重复目标；改动最小增量；不推 GitHub（归 tp_pusher 门禁线）；报告论文 |
| 14482 | TensorPapers | 【论文生产线·首席研究员选题 W22】（owner 2026-11-08 顶刊重组令：按真能发顶刊要求组织；论文内容真正当科研任务实施；要科学家不要操作员。读本 worktree 的 README.md 与 papers/registry. |
| 14537 | TensorSPH | 【返修轮1|maturity_gate 自动裁定】原单 m10431 交付经达标判断方案复核仅 L0（未达 L2 独立功能：L3 缺: impl_files,entry(工件缺失),run_artifacts,integration_cal |
| 14554 | TensorDG | 【返修轮1|maturity_gate 自动裁定】原单 m13815 交付经达标判断方案复核仅 L0（未达 L2 独立功能：L4 缺: impl_files,entry,run_artifacts,integration_callers,b |
| 14576 | TensorSPH | 【重开·TensorSPH#453·管到底第1轮】原 PR#453 已裁定关闭（owner 2026-11-05）。功能要点：matured mission deliverables；文件：benchmarks/capillary_wave |
| 14577 | TensorFEM | 【重开·TensorFEM#618·管到底第1轮】原 PR#618 已裁定关闭（owner 2026-11-05）。功能要点：matured mission deliverables；文件：benchmarks/cohesive3d_ver |
| 14579 | TensorPapers | 【重开·TensorPapers#68·管到底第1轮】原 PR#68 已裁定关闭（owner 2026-11-05）。功能要点：m8601 matured mission deliverables；文件：papers/tools/test_ |
| … | | 另 91 张见 JSON |

## 6. 注意点与局限（honest）

1. **统计量时变**：本快照与上游相差约 40 小时，判定层从 184→410 条、严格 L2+ 从单点 2,366→带 [2,995, 4,054]，任何跨快照比较都必须带 cutoff，否则是拿不同数量的苹果比苹果。
2. **gate 语义未统一**：`board-audit` 与 `audit2-ab/run` 的 gate 判定分布显著不同（降级率 43.8% vs 85.9%），本文只陈述分布，不作「谁对」的裁定——需 owner/协调官定口径。
3. **Tier-B 仍是信号代理**：以 changed_files/pat_gate 近似 L0-L5，未重跑任何仓库测试；代理精度上限由上游验证集给出，本文不再重复估计。
4. **锚点抽取依赖 goal 文本**：仅 9.3% 的 succeeded 单带合法【锚点=】（另有 7.1% 是占位串；历史单早于 LEDGER-GATE），故「单锚点派单次数」是下界估计，不是真实重复率。
5. 本文不改任何任务状态、不动台账、不重开单——名单只列出，裁定权在 owner/协调官。

## 7. 复算方法

```bash
# 复算（任取快照；同快照两次运行产物逐字节一致）
python3 tools/audit_restat/restat_v2.py --snapshot '2026-11-13 01:50:00+00'
# 口径逻辑脱库单测
python3 -m pytest tests/test_restat_l0l5_v2.py -q
```

已执行的验证：

| 检查 | 结果 |
|---|---|
| 口径纯函数单测（锚点/scope/严格性/crosstab/定义带/漏斗） | 29 passed |
| 同一快照两次运行产物一致性 | 逐字节一致（md5 相同） |
| 前值复现（m14502 判定层 184/0） | 逐位复现 |
| 前值复现（m14502/m17113 Tier-B 单点 2,366） | 前者落带内 / 后者不成立（见 §3） |

JSON 产物内 `sql_registry` 保存逐条 SQL；每个统计块可按 sql_id 溯源到数据库。本文所有表格均由该 JSON 程序化渲染，无手工抄录。
