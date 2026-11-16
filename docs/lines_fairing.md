# 型值表三向光顺（曲率能量最小化）

See [`docs/lines_fairing/lines_fairing.md`](lines_fairing/lines_fairing.md) for the method, the acceptance evidence (generated from `results.json`), the figures and the honest limitations.

Quick start:

```python
from tensorcad.offset_forms import offset_table
from tensorcad.lines_fairing import fair_offset_table

table = offset_table('wigley_iii', n_stations=21, n_waterlines=17)
result = fair_offset_table(table, tolerance=0.01 * 0.5 * table.beam)
print(result.flips_before, '->', result.flips_after, result.within_tolerance)

# a power-law section (zeta**(1/m)): declare the near-keel band a design feature, otherwise
# its singular curvature spends the whole deviation budget and the rest barely moves
s60 = offset_table('series60_cb060', n_stations=21, n_waterlines=17)
faired = fair_offset_table(
    s60, tolerance=0.008 * 0.5 * s60.beam, features='auto'
)
print(faired.feature_waterlines)   # (1,) -- detected, not assumed
```
