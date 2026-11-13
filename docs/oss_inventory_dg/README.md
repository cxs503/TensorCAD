# FEAT-OSS-INVENTORY-DG 首波：开源 DG 技术盘点（registry + tier + 7 条缺口锚点）

> 程序化生成骨架 + 人工撰写的六要素章节。数字全部来自 `registry.json`（`wave.py build` 生成，
> 2026-11-13T01:51:06+00:00），不在正文手抄。
> **锚点**：`FEAT-OSS-INVENTORY-DG`　**目标仓**：`TensorDG` @ `281ba5c6cf6777ced440b8df756ad4402287d696`
> **产出**：`registry.json`（机器可解析）· `registry.md`（视图）· `gap-proposals.md`（锚点提案）
> · `evidence/target_grep.json`（现状证据）· `figures/`（2 图）

---

## 一、案例介绍（背景 / 范围 / 参数）

间断 Galerkin（DG）把单元内高阶多项式与单元间数值通量组合起来，同时拿到高阶精度、
局部守恒与 h/p 自适应的可能，是可压缩流、浅水、弹性波、等离子体的通用高阶离散。
对本生态而言 DG 的落点是 `TensorDG`：它已有 1D/2D/3D 对流、Burgers、Euler、
浅水、线弹性、LDG/SIPG 椭圆、人工粘性与 TVB/WENO 限界等算子。

本波的**盘点对象不是物理算例，而是外部技术面**：哪些开源 DG 实现在哪些能力上已经成熟、
其中哪一份可以在许可证允许的前提下被吸收。范围与门槛：

| 维度 | 覆盖 |
|---|---|
| 方法族 | DGSEM · 节点 DG（GLL/GL）· 通量重构 FR/CPR · 谱/ hp 有限元 · 单纯形 DG |
| 能力面 | 14 个功能点（基/求积、数值通量、限界、保界、熵稳定、h/p 自适应、单纯形、曲线单元、时间推进、隐式求解、并行、后验误差、多物理、色散分析） |
| 目标仓证据 | `TensorDG` @ `281ba5c6cf6777ced440b8df756ad4402287d696`，scope=`tensordg/`，逐条 grep（模式/命中数/file:line 样例全部入 `registry.json`） |
| 取数方式 | GitHub REST：每仓 1 次 repo 元数据 + 1 次递归 tree（默认分支），SPDX 许可证由 API 实测，不靠记忆 |

普查规模：池 **23** 条（上游 **22** 条 + 1 条负对照），
其中 T1 **13** / T2 **8** / T3 **1**；
取到 live 元数据 **20** 条，取数失败 **2** 条（如实记录，见 §5）。

---

## 二、怎么做的（方法 / 流程 / 工具链）

1. **协议代码化**：把 LBM/FVM/CAD 三波各自手写的盘点流程收敛成一个零领域语义的引擎
   `tools/oss_inventory_wave/wave.py`（`survey` / `build` / `plot` / `verify`），
   域内容全部外置到 `domains/dg/{wave,pool,gaps}.json`（仓池、功能点签名、锚点、叙事模板）。
   引擎里不含任何 DG 词条——这是「域数据经配置注入」纪律的直接落地。
2. **实测而非回忆**：`survey` 对每个仓取 `stars / license.spdx_id / pushed_at / archived`
   与**递归文件树**；"仓 X 实现功能点 F"由 X **自己的路径名**支持（功能点签名命中计数，
   剔除 tests/docs/vendor），每条保留最多 3 个路径样例作为可引用指针。
3. **负对照**：池里放一条故意不存在的仓（`negctl`），其失败必须出现在 registry 里——
   证明普查会记录失败而不是悄悄丢弃条目。
4. **许可证判定**（本波新增，见 §3.3）：把实测 SPDX 依保守口径映射成
   `port_ok / adapt_at_file_level / read_only_reimplement / read_only_until_resolved`，
   直接决定「可以搬码」还是「只能照论文重写」。
5. **目标仓现状证据**：每条 `GAP-DG-*` 的现状都用 grep 探针表述（`expect: absent|present`
   + 正则 + 命中数 + file:line 样例），在 `TensorDG` 的固定 commit 上跑；
   `verify` 重跑探针并检查 markdown 与 JSON 的往返一致性（禁止手改视图）。
6. **自纠错**：首轮探针被引擎判否两次，都是**我写错了**——
   `hanging` 被 `sign-changing` 误命中、`two-point` 被收敛率的「两点估计」误命中；
   更关键的是首轮把「二维保界」判为 `absent`，重 grep 发现 1D 已有 **Zhang-Shu 缩放限制器**
   （`euler_dg_1d.py:98 _positivity_scale`），于是该锚点从 P0/absent 改判 P1/partial，
   现状描述改为「1D 有 theta 缩放、2D 起只有整格塌到均值的事后过滤」。这条是机器门抓出来的，
   不是靠自觉。

```
python tools/oss_inventory_wave/wave.py survey --domain dg        # 实测取数（缓存可重跑）
python tools/oss_inventory_wave/wave.py build  --domain dg --target /path/TensorDG
python tools/oss_inventory_wave/wave.py plot   --domain dg
python tools/oss_inventory_wave/wave.py verify --domain dg        # 探针重跑 + 视图往返
pytest -q tests/test_oss_inventory_wave_dg.py
```

---

## 三、计算结果（数字表）

### 3.1 Tier 分级（tier 语义见 `wave.json`）

| Tier | 数量 | 含义 |
|---|---:|---|
| `T1` | 13 | 方法标杆实现：读源码即可重实现，许可证允许（或仅需重写） |
| `T2` | 8 | 使能机制库（AMR 森林 / matrix-free 核 / Krylov-IMEX / MPI 抽象） |
| `T3` | 1 | 仅背景：邻接范式（生产级 FVM）或不可复用许可证 |

### 3.2 功能点覆盖（实测路径命中）

矩阵全表见 `registry.md` §5 与 `figures/coverage_heatmap.png`。要点：14 个功能点在池内都有
至少一个仓以自身路径名支持；TensorDG 侧则集中在基/通量/限界/曲线单元/多物理，
在 **单纯形、h/p 自适应、熵稳定、并行、隐式** 五处为零命中。

### 3.3 许可证可用性（决定「能否搬码」）

| 判定 | 数量 | 处置 |
|---|---:|---|
| `port_ok` | 5 | 可直接搬码/改编 |
| `adapt_at_file_level` | 2 | 文件级 copyleft，可改编但需隔离 |
| `read_only_reimplement` | 3 | 只读：照论文重写，禁止拷码 |
| `read_only_until_resolved` | 10 | SPDX 未定（仓库多许可证），先读 `LICENSES/` 再定 |
| `no_metadata` | 2 | 不在 GitHub / 取数失败 |

### 3.4 缺口锚点（本波 7 条，P0 2 条，全部证据探针满足 = 7/7）

| # | 锚点 | 优先级 | 现状 | 参考实现 |
|---|---|---|---|---|
| 1 | `GAP-DG-001` 非结构单纯形单元 | P0 | absent | nektar, ngsolve, mfem, dealii, dolfinx, dune_fem |
| 2 | `GAP-DG-002` h/p 自适应 + 悬挂节点 | P0 | absent | p4est, amrex, philip, trixi, parthenon, ngsolve |
| 3 | `GAP-DG-003` 熵稳定分裂型 | P1 | absent | trixi, flexi, exadg |
| 4 | `GAP-DG-004` 二维保界（Zhang-Shu / 子单元） | P1 | partial（1D 有 theta 缩放） | trixi, flexi, gkeyll, mfem |
| 5 | `GAP-DG-005` MPI / 加速器后端 | P1 | absent | nek5000, nekr, pyfr, libceed, mfem, amrex |
| 6 | `GAP-DG-006` 隐式/IMEX + p-多重网格 | P1 | partial（仅 SIPG 块 Jacobi CG） | exadg, petsc, sundials, dealii, libceed |
| 7 | `GAP-DG-007` 目标导向后验误差（DWR） | P2 | partial（仅伴随灵敏度） | philip, firedrake, dealii |

每条锚点的「现状证据 / 参考实现 / 机器可裁验收判据 / 外部锚（文献）」全文见 `gap-proposals.md`。

---

## 四、结果图

1. `figures/coverage_heatmap.png` — 仓 × 功能点覆盖热图（单元格=该仓源码路径命中数，
   log 色标），直接支持 §3.2；
2. `figures/tier_survey.png` — 按 tier 的池规模与实测成功数（柱内标注 `实测/计划`），
   直接支持 §3.1 与 §5。

**无场量云图的原因**：本波是盘点波，不跑求解器、不产生物理场；所有数字都是
API 实测元数据、路径命中计数与目标仓 grep 探针，因此不给场量云图（下一波 3.4 的
GAP 落地才有物理量可比）。图由 `wave.py plot` 从 `registry.json` 直接渲染，无手绘。

---

## 五、注意点（踩过的坑 / 局限 / 适用范围）

1. **未认证 API 配额**：本波可用 60 次/小时（实测已用 46→14），因此池控制在
   23 条、且每个锚点的参考实现只**实测路径命中**、**不读上游源码**。
   `evidence_policy.curated` 明确标注：tier/focus/upstream_note 属 `curated_knowledge`，
   第 2 波（深盘）必须补读源码，否则「上游怎么做」只能算经验判断。
2. **两个上游不在 GitHub**：`nektar/nektar`（GitLab）与 `dune-fem`（DUNE GitLab）按
   `not_on_github` 记录，不给元数据、不参与 stars 统计——不假装它们不存在，也不编。
3. **第一次判定的错误**（诚实记录）：初版把「二维保界」判 `absent`，探针被
   `Zhang-Shu` 命中打回；改判 `partial`，并把判据从「有没有」改成「2D 触发塌陷比例
   <=1% 且误差退化 <=10%」。这说明**锚点的现状必须让机器重跑**，人工读代码容易只读到一半。
4. **路径命中 ≠ 能力存在**：功能点覆盖是「仓里有没有同名文件」的证据，不是「实现是否正确」
   的证据。它只用于「找候选参考实现」，判级仍以 tier/curated 说明为准。
5. **SPDX 不代表全仓**：10 个仓 API 报 `NOASSERTION`（多许可证/许可证文件非标准），
   一律按 `read_only_until_resolved` 处理，禁止在未读 `LICENSES/` 前搬码。
6. **目标仓在盘点期间前进**：`TensorDG` 从 `348e767` 前进到 `281ba5c6cf6777ced440b8df756ad4402287d696`——锚点的
   commit 钉扎在后者，`verify` 会对比 registry 与实时 commit，漂移即报错。

---

## 六、误差判定

本波不产出物理量，**"关键量 vs 基准误差 <3%" 这道门在本单不适用**（无场量、无解析解可比），
不拿它充数。可裁的门是盘点自身的验收线，逐项核对：

| 验收项（`FEAT-OSS-INVENTORY-DG`） | 门 | 实测 | 判定 |
|---|---|---|---|
| 普查仓数 | ≥20 | 22 条上游 + 1 负对照 | 通过 |
| tier 分级 | 三级且每仓唯一 | T1 13 / T2 8 / T3 1，`verify` 校验取值合法 | 通过 |
| `GAP-DG-*` 锚点 | ≥4 | 7 条（P0 2） | 通过 |
| 锚点现状证据 | 全部探针满足 | 7/7（重跑 `verify` 一致） | 通过 |
| 参考实现可引用 | 每条锚点 ≥2 个 | 每条 ≥3（`verify` 强制 ≥2） | 通过 |
| 视图不可手改 | markdown 往返一致 | `registry.md` / `gap-proposals.md` 与 JSON 重渲染逐字节相等 | 通过 |
| 失败如实记录 | 负对照必须失败入库 | `negctl` 记录为 `HTTP 404` | 通过 |
| 元数据实测率 | — | 20/22 上游；2 条 `not_on_github`、1 条负对照 | 记实 |

**局限（未完成项，如实报）**：上游源码未读（第 2 波深盘）；`registry` 未给每仓 commit 钉扎
（配额不足，仅 `pushed_at` 实测）；未做许可证文本级复核（只取 API SPDX）；锚点验收判据
虽机器可裁但尚未有实现验收（本波只提案，不实现）。