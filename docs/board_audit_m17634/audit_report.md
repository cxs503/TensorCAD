# Board audit m17634 — m17595 `XVAL-CAD-M17483-NATIVE-SIX-FORMAT-MATRIX` (TensorCAD)

交差核验单（非重做）。审计者：`board-audit`。审计单：m17634。被审单：m17595（`status=succeeded`,
`origin=agent`, `task_key=anchor-xval-cad-m17483-native-six-format-matrix`, finished 2026-11-13T07:05:58Z）。

被审交付的**唯一权威指针**：分支 `agent/m17595`，commit `4f5b492`（`xval(m17595): cross-validate native
six-format mesh matrix against meshio/numpy-stl/trimesh/h5py`, 115 files / +276281）。本审计不重做，只复核。

---

## 1. 复核对象与判据（先判据后动手）

| 项 | 内容 |
|---|---|
| 被审能力 | TensorCAD 原生六格式网格矩阵（STL/OBJ/VTK/SU2/CGNS/Gmsh）的交叉对标（XVAL）证据 |
| 机器裁判 | `SELECT case_name,metric,tolerance FROM tensor_benchmark_cases WHERE anchor_slug='XVAL-CAD-M17483-NATIVE-SIX-FORMAT-MATRIX'` → **0 rows**（见 §4 F2） |
| 交付自设判据 | `results.json.acceptance`：`implementation_threshold=1e-3`（体积相对误差）、`owner_threshold=0.15`；拓扑=顶点双射+面多重集相等；定向=相等定向面多重集、flipped=0 |
| 级别判据 | `maturity_gate.SCHEME`（L2=impl_files+entry+run_artifacts；L3 另需 integration_callers；L4 另需 benchmark+repro_cmd；L5 另需 external_anchor） |

## 2. 独立复核方法（方法与独立性）

为了避免"用被审代码验证被审代码"，本次复核**不 import** `tensorcad.mesh_xval`、不运行其断言来取数：

1. 从 `agent/m17595` 克隆出干净副本 `/tmp/m17595_verify`（`git clone --branch agent/m17595`，未改动原工作树）。
2. 复跑其自身门禁：`pytest tests/test_mesh_xval.py -q` → **18 passed in 2.58s**；
   `python3 tools/xval_mesh_six_format.py --check` → **`OK: results.json reproduces`**（exit 0，且 `git status`
   全干净 → 图与 results.json 生成确定可复现）。
3. 用**本审计自写的读者**重算（`docs/board_audit_m17634/audit_independent.py`，第三方库仅 meshio 5.3.5 /
   trimesh 5.0.0 / numpy-stl / h5py 3.15.1 + 自写 SU2 读者 + 自写 CGNS 读者两种方言）：
   读 TensorCAD 写出的 18 个文件 + 18 个 upstream（m17483）文件，独立算体积/三角形数，与解析体积
   （cube=8.0 精确）或六格式中位体积比较。
4. 独立核 `manifest.json` 的 sha256（110 个工件全量）。
5. 独立用 h5py 走 SIDS-to-HDF5 映射规则数 CGNS 违规（不采信其 `sids_violations` 列表）。
6. 查 `tensor_benchmark_cases` 注册面、查 `agent/m17595` 是否已并入生态仓 main。

产物：`audit_independent.py`、`audit_evidence.json`（全部机器数）、`snapshot/`（被审工件与源码字节副本
+ `SHA256SUMS`）、`figures/audit_volume_deviation.png`。

## 3. 三件套核验结果（run_artifacts / benchmark ref vs got / repro_cmd）

| 三件套 | 结论 | 证据 |
|---|---|---|
| run_artifacts 真存在 | **真** | `manifest.json` 声明 110 文件，复核 110/110 sha256 一致、0 缺失、0 不符；doc `docs/xval_mesh_six_format.md` sha256 亦一致；图 3 张（84.6/63.3/48.0 kB） |
| benchmark ref vs got | **真，且在容差内** | 63 个 reader 检查（matrix A 33 + matrix B 30）：ok 48、unreadable 12、nonconformant 3；带数值的 51 个全部 ≤1e-3（自设门）且 ≤0.15（owner 线）；最差 **1.3111681e-8**（matrix A，icosphere2/STL/meshio）；重算复现其 headline "1.311e-8"（相对偏差 <1e-3） |
| repro_cmd 可执行 | **真** | `pytest tests/test_mesh_xval.py -q` → 18 passed；`python3 tools/xval_mesh_six_format.py --check` → 复算无漂移（审计者在干净克隆内实跑） |

审计者**独立重算**（不经过被审代码）与交付数字对照：

| 几何 | 独立重读体积（六格式） | 参考 | 审计重算最大相对误差 | 交付声称 |
|---|---|---|---|---|
| cube | 8.0（STL/OBJ/VTK/Gmsh/CGNS/SU2 全同，解析精确） | 8.0 | **0.000e+00** | 0.000e+00 |
| icosphere2 | 4.047044626915（STL float32）… 4.047044679978（其余五格式） | 4.047044679978 | **1.311e-08** | 1.311e-08 |
| torus64x32 | 1.762301534515…1.762301535676 | 1.762301534515 | **6.912e-10** | 6.912e-10 |

上游 m17483 文件（matrix B）审计侧重读：`obj/msh/su2/CGNS(name-addressed)` 三者体积完全一致
（cube 8.0、icosphere 4.04704468、torus 1.7623015345）→ **几何正确**；`vtk`（POLYDATA）meshio 读不了，
与交付标注的 `unreadable` 一致。

## 4. 主要发现（按严重度排序）

**F1（minor，措辞精度）"独立实现重读六格式"把 4 个第三方读者与 2 个自写读者并称。**
矩阵 A 的 33 个检查里，真第三方几何 oracle 只覆盖 STL(3 个读者)/OBJ(2)/VTK(1)/Gmsh(1)；
CGNS 只有 `h5py (SIDS mapping checker)` 结构检查 + 其自写 `read_cgns_sids` 读几何，SU2 只有其自写
`read_su2_strict`（meshio 5.3.5 无 SU2 reader；其 CGNS reader 硬编码 `Base/Zone1/GridElements` 且只支持
四面体，见 `meshio/cgns/_cgns.py:20-38`，被审方亦有测试 `test_meshio_cgns_reader_is_a_stub` 钉住）。
减轻情节：表头 reader 列写得准确、正文 §5/§6 有差异清单，只有 headline/摘要把"六格式"整体称作
OSS 交叉读数。审计者用自写 CGNS/SU2 读者独立重算出同一体积（8.0 / 4.04704468 / 1.7623015345），
**数字本身不依赖被审代码即可复现**，故非阻断。建议：摘要改为"4/6 格式第三方 oracle + 2 格式（CGNS/SU2）
自写读者并由审计侧独立重导"。

**F2（moderate，裁判注册面缺失）该锚点在 `tensor_benchmark_cases` 无任何行（0 rows）。**
因此本单的"唯一裁判"空置：1e-3/0.15 判据只活在交付自己的 `results.json` 与报告里，机器无法独立断言
通过与否。另有旁证：兄弟锚点 `GAP-CAD-B2-MESH-FORMAT-MATRIX` 下确实登记了
`CASE-CAD-B2-ROUNDTRIP-18`（`created 2026-11-13`，tolerance 文本里就写着 "worst measured 1.311e-08, STL
float32"，metric 由 `tensorcad.mesh_formats.compare_meshes` 内部自检产生）——同一 18 格矩阵、同一最差数，
但挂在**另一个锚点**上且未与本 XVAL 锚点互链。后果：XVAL 系列的"基准裁判"实际上是被审交付的自述。
建议：XVAL feeder 落单时同步登记 case（metric 用第三方读者重读数，tolerance 用 1e-3），并补 run 时回填
`passed_at`。

**F3（moderate，集成/主线性）交付未并入生态仓 main，能力对仓库读者不可见。**
`4f5b492` 只在 `agent/m17595` 上（`git branch --contains` 仅此分支）；生态仓 `/root/tensor-ecosystem/TensorCAD`
main（`928f3e3`）无 `tensorcad/mesh_xval.py`，`grep -rn mesh_xval --include=*.py` 在 main 0 命中；无 CI job、
无 docs 索引引用。仓库内唯一引用是 `tools/xval_mesh_six_format.py` 与 `tests/test_mesh_xval.py`（均在分支上）。
故 `integration_callers` 为空（见 §6 级别裁定的直接原因），且工作树被清理则该证据链只剩本审计的字节副本。

**F4（minor，规范精度）"validated reference dialect" 略高于实证。**
审计者独立按 SIDS-to-HDF5 规则走查：upstream `cube.cgns` 10 个节点 **10/10 缺 name、10/10 缺 label、
10/10 缺 order、10/10 把 label 塞在 `type` 槽**，根上亦无版本信息 —— 被审方的"映射违规、几何无损"
结论**独立成立**；但其自写方言只做到 name/label/order/type 齐全 + 根上 `format`/`version` **属性**，
**没有 SIDS 要求的 `CGNSLibraryVersion` 节点**，而其检查器只查了根属性存在性。建议：报告里把"validated
reference dialect"限定为"满足本仓 CGNS_RULES 三条规则"，并补 `CGNSLibraryVersion` 节点或写明豁免依据。

**F5（positive，honest-FAIL 处理正确）** 12 个 `unreadable` + 3 个 `nonconformant` 全部显式标注、不计入
pass；报告 §6.2 的 15% 误差表有 **两行显式 NO**（STL 在 1000 m 坐标偏置下坐标偏差 3.85e-5 m → 需局部坐标系；
上游 CGNS 12 违规/文件 → 结构问题而非几何）。审计者复核这两条均能独立重现（见 §5 未跑清单中的偏置项）。

**F6（minor，计数口径）** `results.json` 的 `matrix_b.summary.reader_checks=30 / green=21` 与
`matrix.md` 一致；但 `reader_checks_green` 把 3 个 CGNS `nonconformant` 计为绿（几何可读），报告若被
外部读者当作"全通过"会误读。建议 summary 里把 `nonconformant` 单列。

## 5. 未跑/未复核的检查（明确标 unknown，不算通过）

- **未跑**：`cgnslib`（本机无 `libcgns`/`CGNS` python 模块）→ CGNS 几何缺"真第三方消费者"一侧的验证；
  审计侧只能用自写 SIDS 读者重导，不能替代官方库。
- **未跑**：1000 m 坐标偏置敏感性重算（依赖其 `triangulate`/offset 链路；该项被审方已自报为 NO，
  审计侧未复现该条，故其具体数值 3.85e-5 m 标 **not_run/unknown**）。
- **未跑**：收敛阶梯（icosphere sub1…4、torus 16x8…128x64）的逐点重算；只核了它进入 `results.json`
  与 `matrix.md` 的一致性与其单调性叙述。
- **未跑**：HPC/CI 侧回归（该仓此能力无 CI 接线，见 F3）。
- **未复核**：`upstream_m17483/files/*` 是否与 m17483 原始工件逐字节相同（只用了交付内副本 + sha256 自洽）。

## 6. 级别裁定

**裁定 L2（独立功能实现），不升 L3/L4，也不降到 L1。**

- L1 通过：`tests/test_mesh_xval.py` 18 项（含 `test_compare_detects_wrong_orientation_and_missing_vertex`
  这类反例测试、`test_committed_results_match_the_report` 把提交工件与报告数字钉在一起）。
- L2 通过：`tensorcad/mesh_xval.py`（823 行，占位模式命中 1 处且为异常分支 `pass`）有真实现；
  入口 `tools/xval_mesh_six_format.py`（705 行，`--check` 漂移门）实跑可用；运行工件 110 个 sha256 全验。
- L3 **未达**：`integration_callers` 空——主求解器/平台侧无调用方，交付未并入 main（F3）。这不是"再做一点"
  的措辞问题，而是机器 grep 可复现的缺失。
- L4 **未达（按积分制被 L3 卡住）**：基准面本身其实是 L4 档（完整算例 18 格 × 2 矩阵、收敛/偏置/变体研究、
  误差分析、`--check` 可复现、第三方 oracle），唯一缺的是注册面（F2）与集成（F3）。
  若后续把工具接进仓内流水线（例如 CI 里跑 `--check`）并在 `tensor_benchmark_cases` 落 case，可直接冲 L4。

## 7. 处置（disposition）

**approve with notes**：交付为真、数字经独立重算成立、可复现、honest-FAIL 处理正确；级别按机器判据
落 **L2**（`integration_callers` 与裁判注册面缺失，已写成可复查的差异清单）。返修建议按 F1/F2/F3/F4/F6
逐条，均为"补证据/补接线"而非"改数字"。不主张 reopen：无造假、无虚标数值。

---

### 复现指令（审计侧）

```bash
git clone --branch agent/m17595 /root/tp-agent/work/m17595/.git /tmp/m17595_verify
cd /tmp/m17595_verify && python3 -m pytest tests/test_mesh_xval.py -q          # 18 passed
python3 tools/xval_mesh_six_format.py --check                                  # OK: results.json reproduces
cd /root/tp-agent/work/m17634
SUB=docs/board_audit_m17634
python3 $SUB/audit_independent.py     # → audit_evidence.json（独立重算，不 import 被审代码）
python3 $SUB/audit_figure.py          # → figures/audit_volume_deviation.png
```

### 误差判定（对 owner 15% 红线逐项回答"是否 <3%"）

| 关键量 | 参考 | 独立复核 got | 误差 | <3%? | <15%? |
|---|---|---|---|---|---|
| cube 体积 × 6 格式 | 8.0（解析精确） | 8.0 | 0.000e+00 | 是 | 是 |
| icosphere 体积（最差格式 STL） | 4.047044679978 | 4.047044626915 | 1.311e-08（1.3e-6 %） | 是 | 是 |
| torus 体积（最差格式 STL） | 1.762301534515 | 1.762301535676 | 6.912e-10 | 是 | 是 |
| 上游 m17483 CGNS 几何（name-addressed 重读） | 上游 obj/msh 体积 | 完全一致 | ≤1e-9 | 是 | 是 |
| 上游 m17483 CGNS 结构合规 | 0 违规 | 10/10 节点缺 name/label/order | 结构违规 | **否**（非几何量，见 F4） | 不适用 |
| STL 在 1000 m 偏置下坐标偏差 | 0 m | 3.85e-5 m（被审方自报，审计未复现） | not_run | **未知** | 不适用 |

报告图：`figures/audit_volume_deviation.png`（A：审计侧重算的六格式体积相对误差热图；B：每格式的
几何 oracle 归属——第三方 vs 自写，正是 F1 的图证）。