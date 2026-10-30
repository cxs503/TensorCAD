# m585: 变形族可视化出图 + 报告升级 — 完成报告

## 任务概述

对 m542 (DTMB 5415) 和 m546 (KCS) 艏部参数化变形族进行可视化出图与报告升级，将原 500 字符 stub 报告 498/499 升级为论文级报告。

## 交付物清单

### 1. 可视化图（6 张，全部本单真跑）

| 图号 | 文件名 | 内容 | 大小 |
|------|--------|------|------|
| Fig 1 | `artifacts/fig01_side_overlay.png` | DTMB 5415 & KCS 侧视叠加图（基线 vs 6 代表变体，标变形区） | 120 KB |
| Fig 2 | `artifacts/fig02_top_overlay.png` | DTMB 5415 & KCS 俯视叠加图（水线面半宽曲线） | 156 KB |
| Fig 3 | `artifacts/fig03_non_deform_deviation.png` | 非变形区偏差验证（30+32 变体全零偏差） | 88 KB |
| Fig 4 | `artifacts/fig04_param_distribution.png` | 排水量分布 + 参数空间散点图 | 253 KB |
| Fig 5 | `artifacts/fig05_gate_heatmap.png` | 三门禁热力图（62 变体 × 3 门禁 = 全绿） | 76 KB |
| Fig 6 | `artifacts/fig06_3d_comparison.png` | 3D 线框对比（基线 vs v015） | 978 KB |

### 2. 变体库状态

| 船型 | 变体数 | 三门禁全通过 | 通过率 | 门禁结果文件 |
|------|--------|-------------|--------|-------------|
| DTMB 5415 (m542) | 30 | 30 | 100% | `deform5415/family_summary.json` + 每变体 `gates.json` |
| KCS (m546) | 32 | 32 | 100% | `kcs_variants/variant_summary.json` |

### 3. 升级报告

- **报告 ID**: 通过 tensor_reports API 提交（mission_id: 585）
- **Markdown 长度**: 3624 字符（≥2500 要求 ✅）
- **内嵌图**: 6 张 base64 编码 PNG（figures 非空 ✅）
- **本地副本**: `docs/report_m585.json`

### 4. 复现脚本

- `generate_figures.py`: 完整可视化脚本，可从 m542/m546 产物复现全部 6 张图

## 判据验证

| 判据 | 要求 | 实际 | 状态 |
|------|------|------|------|
| 图数量 | ≥4 | 6 | ✅ |
| 图来源 | 本单真跑 | generate_figures.py 从 m542/m546 STL 真跑 | ✅ |
| 变体门禁结果 | 在 | family_summary.json + variant_summary.json | ✅ |
| 报告长度 | ≥2500 | 3624 | ✅ |
| 报告 figures | 非空 | 6 张 base64 PNG | ✅ |
| ≥30 变体三断言全绿 | 是 | 5415: 30/30, KCS: 32/32 | ✅ |

## 关键数字（零手抄，程序化注入）

- DTMB 5415 基线体积: 0.481435 m³
- DTMB 5415 基线 Cb: 0.480043
- KCS 基线排水量: 49440.17 t
- KCS 基线 Cb: 0.61812
- 非变形区最大偏差: <1e-6 m（全部 62 变体）
- 门禁通过率: 100%（186/186 检查项）

## 局限

1. STL 网格分辨率约 6000-8000 三角面，CFD 应用需加密
2. 变形区边界 C¹ 连续但 C² 可能有微小跳跃
3. 基线水线面面积计算为零（z=0 离散问题）
