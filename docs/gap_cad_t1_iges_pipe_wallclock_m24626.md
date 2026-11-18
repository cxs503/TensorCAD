# GAP-CAD-T1-IGES-PIPE — m24626 墙钟实测：KVLCC2 过门，ONRT 精确定位（criterion 2 欠账闭环）

**Anchor** `GAP-CAD-T1-IGES-PIPE` · **Mission** m24626 · **Repo** TensorCAD · **Branch** agent/m24626
**日期** 2026-11-17 · **上游基线** 0e944ca7（origin/main）

## 1. 本单干了什么

m23992 §2.2 与 m24043 §3 把同一笔欠账留给了下次自找活：

> 「142 字段位解析修复与 NURBS 导数 memoize 已合并，但 KVLCC2/ONRT 是否已落入 120 s 门内
> 未在集群完成测量（集群节点无资产字节权限）。下次自找活：在有资产权限的节点跑
> `test_real_iges_gate.py` + `test_iges_pipeline_gate.py` 计时。」

本单在有资产权限的环境（控制面本机 + 外协方 loong64 节点 n330/n315）把这笔墙钟欠账
**实测闭环**：

1. **KVLCC2 转绿**：`import_cad(kvlcc2.igs)` 墙钟 **17.5 s**（loong64 n330）/ 17.8 s（本机
   x86），19/19 面全部 tessellate、0 拒绝、最大反演残差 0.0，G1/G2/G5 全过。此前状态是
   「>58 s 未完成、无结果」（inversion_speedup 文档 §5）——现在是有结果且在 120 s 门内。
2. **ONRT 精确定位**（此前状态是「>120 s 无结果」）：墙钟 **414.1 s**（超 120 s 门），
   132 面中 **12 面被拒**，G1/G2 不过、G5 过。12 个拒绝逐面点名（见 §3），失败模式收敛为
   两类：(a) 反演残差 1.9e-04…3.2e-04 略超 transfer budget 1.0e-04（2–3 倍）；
   (b) 反演后的环自交，耳切法报「not a simple polygon」。
3. **回归钉死**：新增 `tests/test_real_iges_wallclock_gate.py`——三个真实 IGES 逐文件断言
   120 s 墙钟门 + G1/G2/G5；资产缺失自动 skip（与 `test_real_iges_gate.py` 同约定）。
   KCS/KVLCC2 两例本机实测 22 s 通过；ONRT 例在本机/集群当前均不过门，故**不**写入断言
   （honest：不粉饰），ONRT 留作新总账条目（§5）。

## 2. 实测数字（程序化，零手抄）

测量工具 `work/m24626_wallclock/measure_wallclock.py`（逐文件计时 + 门控摘录 → JSON）。
集群作业 **cj5181**（n330，loong64，stage_dir 直传 19 MB 工件包）产出
`work/m24626_wallclock/cj5181/cj5181/repo/wallclock.json`；本机 x86 复测 KVLCC2 见
`work/m24626_wallclock/wallclock.json`。

| 文件 | 面数 | 拒绝 | 最大反演残差 | 墙钟（loong64） | 120 s 门 | G1 | G2 | G5 |
|---|---:|---:|---:|---:|---|---|---|---|
| KCS_hull_SVA.igs | 55 | 0 | 2.6e-09 | **1.6 s** | ✅ | ✅ | ✅ | ✅ |
| kvlcc2.igs | 19 | 0 | 0.0 | **17.5 s** | ✅ | ✅ | ✅ | ✅ |
| onr5613-…mergedR.igs | 132 | **12** | 5.1e-03 | **414.1 s** | ❌ | ❌ | ❌ | ✅ |

本机 x86 交叉验证：KVLCC2 17.8 s（与 loong64 17.5 s 一致，±2%）；KCS 1.6 s。

### 2.1 KVLCC2 过门细节（criterion 2 的欠账主体）

- 单位门：文件声明 MM、`scale_applied=0.001`（坐标实为 m，与 m24043 单位门结论一致）；
- 19/19 面全部走 142 字段位解析后的精确参数环路线（`reader_refused` 里 19 条记录全部是
  「3-D 曲线仅作 cross-check、参数空间像被采用」的**说明性**记录，不是拒绝）；
- 103 915 三角（修复前）→ 103 717（修复后）；
- G3/G4/G6 仍不过（3 自由边、37 非流形边、27 连通块）——这是 trimmed 路线的**已知**
  面间缝问题（`docs/iges_pipe_limitations.md`），水密船体由 fleet_baseline 半宽场路线
  交付（m24043：体积 312 568 m³ vs 公开 312 622 m³ = −0.017%，VERIFIED）。

## 3. ONRT 失败画像（honest FAIL，逐面点名）

12/132 面被拒，两类失败模式：

| 模式 | 面 | 证据（report.import.refused） |
|---|---|---|
| 反演残差略超 transfer budget | 5, 23, 32, … | `inversion residual 1.9e-04…3.2e-04 > transfer budget 1.0e-04 (exact budget 1.0e-06)` |
| 反演环自交 | 46, 60, … | `ear clipping found no ear: the trim loop is not a simple polygon` |

完整 12 条拒绝记录与逐面统计：`work/m24626_wallclock/cj5190/…/onrt_full_report.json`
（集群作业 cj5190）。

## 4. 回归与不回归

| 检查 | 结果 |
|---|---|
| `pytest tests/test_real_iges_wallclock_gate.py -k "not ONRT"`（本机，资产齐全） | **2 passed**（22 s） |
| `pytest tests/test_real_iges_gate.py tests/test_iges_gold_regression.py`（本机） | **13 passed, 1 skipped**（基线不变） |
| 金标回归 DTMB5415（criterion 3） | 绿（`test_iges_gold_regression.py` 9 测试，含在本行） |
| 依赖锁定（criterion 4） | `requirements-iges-pipe.lock` 未动 |

本单**未改任何 `tensorcad/` 源码**——纯验证欠账闭环 + 回归钉死，零算法语义变更。

## 5. 给下次自找活的建议线索

1. **首选（新总账条目）**：`GAP-CAD-T1-IGES-PIPE-ONRT`，所属域 CAD 进口。
   验收 = ONRT trim_mesh 路线 132/132 面过 G1/G2 且墙钟 <120 s；
   两条子线：(a) 残差 2–3× 超预算的面的种子策略（seeded walk 落点离解太远）；
   (b) 自交环的环简化/重采样后再耳切。
2. 次选：trimmed 路线面间缝缝合（KCS 795 自由边 / KVLCC2 3 自由边 + 37 非流形边），
   让 G3/G4/G6 对真实文件可达——这是「全过门」的完整语义，工作量大于 ONRT 单点。
3. 若 ONRT 修绿：criterion 2 整体关闭，锚点可转 done。
