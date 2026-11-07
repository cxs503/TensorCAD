# Series 60 Hull Family

The Series 60 (Todd 1951) is the standard systematic displacement-hull benchmark.
Five parent forms cover Cb = 0.60–0.80 with closed-form hydrostatic coefficients.

## Geometry

Half-breadth: `y(xi, zeta) = (B/2) * (1 - |2*xi-1|^n) * zeta^(1/m)`

`zeta` is height above the keel over the draft. The exponent is `1/m`, not `m`:
`int_0^1 zeta^(1/m) dzeta = m/(m+1) = Cm`, which is what the closed forms below
(and the published `Cp` values that fix `n` and `m`) require. With the exponent
inverted the lofted sections collapse to a sliver over the bottom 90 % of the
draft and the mesh integrates ~21 times less volume than its own analytic
reference -- a failure that survives refinement (see "Measured accuracy" at the
end of this page, and `tests/test_series60.py::TestVerticalTaper`).

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
| VCB above keel | `T·(m+1)/(2·m+1)` |
| Cb | `n·m/((n+1)(m+1))` |
| Cp | `n/(n+1)` |
| Cm | `m/(m+1)` |
| Cwp | `n/(n+1)` |

VCB is the first moment of `zeta^(1/m)` over its own area:
`int zeta*y dzeta / int y dzeta = (m+1)/(2m+1)`. The earlier `(m+1)/(2(m+2))`
was 6.5 % low at m=24 (0.481 T against 0.510 T) and passed every "0 < VCB < T"
sanity check, which is why it needed an independent recomputation to find.

## Usage

```python
from tensorcad.series60 import Series60Parameters, generate_series60, analytic_series60_hydrostatics

p = Series60Parameters(length=10, beam=1.5, draft=0.6, cb=0.60)
model = generate_series60(p)
ref = analytic_series60_hydrostatics(p)
```

## Fidelity to the published parents

`Cwp` equals `Cp` identically in this family and `LCB` sits at `L/2`; both are identities of
the single-exponent taper, and against the tabulated parents that is −10.7 % on `Cwp` at
`Cb = 0.60`. Pass `Series60Parameters(cwp=..., lcb=...)` to set them independently, or read
`docs/series60_fidelity.md` for the measurement and the verdict.

## Measured accuracy

Mesh displacement against the closed form, `L=10 m, B=1.5 m, T=0.6 m`, with the
`zeta^(1/m)` section (numbers from `examples/geometry_audit/run_audit.py` and the
grid study in `docs/m5283_geometry_audit.md`):

| Grid (nx × nz) | Cb=0.60 volume error | Cb=0.70 volume error |
|---|---|---|
| 41 × 13 | -3.68 % | -3.95 % |
| 81 × 21 | -2.14 % | -2.32 % |
| 161 × 41 | -1.04 % | -1.14 % |

The residual is the loft's resolution of the steep section near the keel: the
profile `zeta^(1/m)` with m = 24–47 rises to most of its height within the first
few percent of the draft, so the trapezoidal-in-z sample of the loft under-fills
the bottom. It converges at first order in `1/nz`; `estimate_series60_discretisation_error`
reports the convergence rather than asserting a tolerance, and
`Series60Parameters(nx=161, nz=41)` is the resolution used for the 2 % claim in
the audit.

Note that the parent form has a **zero-radius bilge**: the section turns from the
flat of bottom to the flank in one mesh cell, so the Gaussian curvature field is
singular along the keel line. The geometry audit's fairness gate reports that (as
it should) for the fuller parent forms; a production hull needs a finite bilge
radius.