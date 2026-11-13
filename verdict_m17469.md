# 交差审计报告：m17469 GAP-CAD-D1-MESH-QUALITY-LIB

**审计任务**：m17578（board-audit 交差核验）
**审计对象**：m17469（GAP-CAD-D1-MESH-QUALITY-LIB·网格质量指标库）
**审计日期**：2026-11-13
**审计人**：board-audit

## 执行摘要

m17469 交付了 TensorCAD 的 `tensorcad/meshquality` 包（1566 行代码），实现 tri/quad/tet/hex 四类单元的质量指标（scaled_jacobian、min_sicn、Gmsh gamma/eta/gmsh_icn、角度/面积/体积/边比），含统计量、机器门、JSON 交换格式与 CLI。26 个单元测试全部通过，基准测试与解析值/独立实现/trimesh 三方交叉验证通过。

**判级：L3**（系统集成）

## ① 取证核验

### 三件套

| 项目 | 状态 | 证据 |
|---|---|---|
| run_artifacts | ✅ | `artifacts/mesh_quality/{results.json, hull_hex.json, hull_surface_tri.json, manifest.json, plot.png}` |
| benchmark ref vs got | ✅ | `run_benchmark.py` → criteria.passed=true；analytic 1e-12, reference 1e-7~1e-9, trimesh 1e-9 |
| repro_cmd | ✅ | `python3 examples/mesh_quality/run_benchmark.py` (exit 0); `pytest tests/test_meshquality.py -v` (26 passed) |

### 报告六段：全部齐全

## ② 代码审查

- **实现文件**：6 个文件 1566 行，非占位，定义清晰
- **测试**：26/26 通过，覆盖解析值/独立实现对齐/CLI/门/边界
- **集成**：validation.py 导入 tri_metrics；CLI 可作 CI 门
- **基准**：9 个网格，所有指标差异远小于容差

## ③ 判级

| 级别 | 判定 | 理由 |
|---|---|---|
| L0 | ✅ | 文档齐全 |
| L1 | ✅ | 26 单元测试 |
| L2 | ✅ | 1566 行真代码 |
| L3 | ✅ | validation.py 集成 + CLI CI 门 |
| L4 | ⚠️ | 完整算例+误差分析+可复现 ✅，但 tensor_benchmark_cases 表不存在 |
| L5 | ❌ | |

**最终判级：L3**

## ④ verdict 落表

- verdict id=44
- `INSERT INTO tensor_maturity_verdicts (repo, scope, mission_id, level, evidence, verdict_by) VALUES ('TensorCAD', 'm17469 GAP-CAD-D1-MESH-QUALITY-LIB', 17578, 3, '<JSON>', 'board-audit')`

## ⑤ 关键发现

1. 代码质量高，Gmsh/Verdict 出处可追溯
2. 交叉验证严谨，未安装工具诚实标记 not_run
3. tensor_benchmark_cases 表缺失（平台缺口）
4. meshquality 未导出到 tensorcad/__init__.py

## ⑥ 误差判定

| 对比项 | max diff | 容差 | 通过 |
|---|---|---|---|
| analytic | 0 | 1e-12 | ✅ |
| reference_tri | 1.4e-14 | 1e-9 | ✅ |
| reference_tet | 1.9e-12 | 1e-9 | ✅ |
| reference_hex | 1.3e-08 | 1e-7 | ✅ |
| trimesh | 7.1e-15 | 1e-9 | ✅ |

**所有误差均 << 容差。**

## ⑦ maturity_gate

maturity_gate API 端点不存在（404），无法自动裁定。verdict 已落表，待平台侧实现 gate 逻辑后自动处理。
