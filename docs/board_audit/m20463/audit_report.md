# Board Audit Report: m20389 REPRO-GEOM-FLEET-REFS

**审计任务**: m20463  
**被审任务**: m20389  
**锚点**: REPRO-GEOM-FLEET-REFS（六标模文献交叉参考值库：每船 ≥2 独立来源）  
**仓库**: TensorCAD  
**被审commit**: 1f5763d (branch agent/m20389)  
**判级**: **L2**  
**门禁建议**: **reopened**

---

## 审计范围

1. 取证：核验三件套（run_artifacts / benchmark ref vs got / repro_cmd）
2. 判级：L0-L5
3. 落表：INSERT INTO tensor_maturity_verdicts

## 三件套核验

### run_artifacts ✅ 存在
- `benchmarks/verified/fleet_refs/crossref.json` (1588行, schema完整)
- `docs/fleet_refs_crossref.md` (115行, 六段报告)
- `docs/figs/fleet_refs/` (3张PNG: closure_check/coverage_matrix/displacement_ref_vs_got)

### benchmark ref vs got ⚠️ 部分通过
| 船 | ref ∇ (m³) | got ∇ (m³) | 误差 | ref判定 | benchmark |
|---|---|---|---|---|---|
| KCS | 52,030 | 52,188 | **+0.30%** | CROSS_VERIFIED | ✅ VERIFIED |
| KVLCC2 | 312,622 | — | — | SINGLE_SOURCE | WITHHELD |
| JBC | 178,370 | 145,769 | -18.28% | SINGLE_SOURCE | WITHHELD |
| ONRT | 8,300 | 93 | — | CONFLICT | WITHHELD |
| DTC | — | — | — | MISSING | N/A |
| Series60 | — | — | — | MISSING | N/A |

**仅KCS(1/6)通过benchmark**。

### repro_cmd ⚠️ 部分可执行
- 文档命令完整
- 测试：**12 passed, 2 failed** (0.46s)
- 失败项为实现bug：
  1. `test_fallback_column_is_flagged_and_does_not_vote`：fallback列未排除投票
  2. `test_conflict_is_reported_not_averaged`：family_values键用source key而非family

## 报告六段 ✅ 齐全
1. 案例介绍 ✅
2. 方法 ✅
3. 计算结果 ✅
4. 结果图 ✅ (3张matplotlib)
5. 注意点 ✅
6. 误差判定 ✅

## 判级理由

**L2**：独立功能实现（真代码1142行：HTML解析→谱系分类→交叉验证→闭锁自检→报告渲染）

不达L3：无求解器/平台集成调用者（仅tests和tools调用）

不达L4：
- 2/14测试失败（实现bug未修复）
- 仅1/6船满足≥2独立来源要求
- 3/6船benchmark WITHHELD

## 正面发现
- 诚实度高：主动暴露KVLCC2 Bwl=580排版错误、ONRT质量/体积混淆
- SHA256钉扎出版商字节，每格可溯源
- 谱系家族判定合理（SIMMAN 2014/2020正确合并为同一源）
- Cb闭锁自检为独立数学验证

## 负面发现
- 2/14测试失败未修复即交差
- 仅KCS(1/6)满足锚点要求
- 门禁应reopened

## 产物
- `docs/board_audit/m20463/insert_verdict.sql` — 落表SQL
- `docs/board_audit/m20463/verdict_evidence.json` — 完整证据
- `docs/board_audit/m20463/audit_summary.png` — 审计图
- 报告已发布: report_id=18622
