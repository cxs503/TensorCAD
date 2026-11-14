# GAP-CAD-T1-IGES-PIPE: IGES/STEP → 水密网格进口管道

## 状态总结

**锚点**: `GAP-CAD-T1-IGES-PIPE`
**验收标准**: 管道 CLI+库函数；3 个真实 IGES 全过门；金标回归测试绿；依赖锁定

## 当前能力

### 路线 A: fleet_baseline（半宽场 → 构造水密）

| 船型 | 状态 | 体积误差 | 水密 | 备注 |
|------|------|----------|------|------|
| KCS | ✅ VERIFIED | +0.30% | ✅ | 55 个 type-128 面片，半宽场重建 |
| KVLCC2 | ❌ WITHHELD | N/A | N/A | ValueError: fewer than two usable waterline slices |
| JBC | ❌ WITHHELD | -18.28% | ✅ | IGES layout guard fired at -18.3% |
| ONRT | ❌ WITHHELD | -98.91% | ✅ | model-scale units issue |

### 路线 B: trim_mesh（面片 tessellation + repair）

| 船型 | 面片数 | 拒绝数 | 水密 | 自由边 | 备注 |
|------|--------|--------|------|--------|------|
| KCS | 55 | 0 | ❌ | 798 | 面片间缝隙无法通过焊接修复 |

### 金标回归

- DTMB5415_Hull_exact: ✅ 9/9 测试通过

### 依赖锁定

- 核心包 `dependencies = []` ✅
- IGES 阅读器 (iges.py) 仅使用标准库 ✅

## 缺口分析

### KVLCC2 水线切片 bug

`fleet_baseline.py` 的 `slice_profile` 函数对 KVLCC2 返回少于两个可用切片。
原因可能是 KVLCC2 的 IGES 文件面片布局与 KCS 不同（38 vs 55 面片），
导致半宽场重建时某些水线层没有数据。

### JBC 体积误差 -18.28%

JBC 的 IGES 文件触发了 layout guard，说明面片控制点布局解析有问题。
JBC 文件包含 type-126（多项式 B-spline）面片，而 KCS 只有 type-128（有理 B-spline）。

### trim_mesh 路线水密性

trim_mesh 路线对面片独立 tessellation，相邻面片的边界采样不匹配，
导致 T-vertex 和裂缝。焊接容差无法解决结构性缝隙。

## 测试

新增测试文件: `tests/test_iges_pipeline_gate.py`

```bash
pytest tests/test_iges_pipeline_gate.py -v
# 7 passed, 1 skipped (KVLCC2 expected failure)
```

## 下一步建议

1. **KVLCC2 水线切片修复**: 检查 `slice_profile` 对 KVLCC2 面片的采样逻辑
2. **JBC type-126 支持**: 完善 `iges.py` 对 type-126 面片的控制点布局解析
3. **trim_mesh 面片缝合**: 实现基于边界匹配的面片缝合算法，而非简单焊接

## 工件

- `tests/test_iges_pipeline_gate.py` — 端到端门控测试
- `docs/iges_pipeline_status.md` — 本文档
