# 落表回执 · m17405（board-audit on m17357）

| 项 | 值 |
|---|---|
| 表 | `tensor_maturity_verdicts` |
| id | **543**（INSERT 0 1, 2026-11-13） |
| repo / scope | `TensorCAD` / `m17357 GAP-CAD-D3-GEOM-CONTRACT` |
| mission_id | 17405（本审计单） |
| level | **2**（L2 独立功能实现；L3/L4 明确拒绝，理由见 `evidence.level_refused`） |
| verdict_by | `board-audit` |
| gate_state | `pending` → maturity_gate 守护（≤30 min 周期）裁定 |
| gate 预演（真裁判函数 `/root/tp-agent/maturity_gate.py:judge`） | `judge(2, evidence, 'TensorCAD', 17405)` → **gate_level 2, confirmed**，无 fail reason |
| 六要素报告 | `tensor_reports` id=**16381**（published, figures=contract_audit.png） |
| 证据包 | `evidence.json`（20 键：8 机器判据键 + pointers + checks_not_run + commit_state + sha256） |
| 字节快照 | `snapshot/`（6 文件 + MANIFEST.json，对象未提交交付物的 sha256 留痕） |

复现链：`audit_probe.py`（探针）→ `make_evidence.py`（证据）→ `insert_verdict.py`（落表，含 gate 预演）
→ `make_audit_doc.py` / `make_figure.py`（报告与图）→ `post_report.py`（发布）。
