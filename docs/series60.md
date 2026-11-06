# Series 60 Hull Family

The Series 60 (Todd 1951) is the standard systematic displacement-hull benchmark.
Five parent forms cover Cb = 0.60–0.80 with closed-form hydrostatic coefficients.

## Geometry

Half-breadth: `y(xi, zeta) = (B/2) * (1 - |2*xi-1|^n) * zeta^m`

Exponents derived from target Cb via published Cp values:
- Cb=0.60 → Cp=0.625, n=1.667, m=24.0, Cm=0.960
- Cb=0.65 → Cp=0.670, n=2.030, m=32.5, Cm=0.970
- Cb=0.70 → Cp=0.715, n=2.509, m=46.7, Cm=0.979
- Cb=0.75 → Cp=0.760, n=3.167, m=75.0, Cm=0.987
- Cb=0.80 → Cp=0.805, n=4.128, m=160.0, Cm=0.994

## Closed Forms

| Quantity | Formula |
|---|---|
| Volume | `L·B·T·n·m/((n+1)(m+1))` |
| Waterplane area | `L·B·n/(n+1)` |
| Midship area | `B·T·m/(m+1)` |
| VCB above keel | `T·(m+1)/(2·(m+2))` |
| Cb | `n·m/((n+1)(m+1))` |
| Cp | `n/(n+1)` |
| Cm | `m/(m+1)` |
| Cwp | `n/(n+1)` |

## Usage

```python
from tensorcad.series60 import Series60Parameters, generate_series60, analytic_series60_hydrostatics

p = Series60Parameters(length=10, beam=1.5, draft=0.6, cb=0.60)
model = generate_series60(p)
ref = analytic_series60_hydrostatics(p)
```

## Note on Discretisation

Series 60 sections are nearly rectangular (U-shaped) with large m values.
The default 41×13 grid gives ~95% volume error because the vertical taper
concentrates breadth near the waterline. Use `estimate_series60_discretisation_error`
to find the resolution needed for your accuracy target. At 161×41 the error
drops to ~5%; at 641×101 it reaches <1%.

This is a property of the hull form, not a defect: U-section hulls inherently
require finer waterline resolution than V-section forms like Wigley.
