# Transom-stern hull family

> Status: implemented in `tensorcad/transom.py` (0.14.0), tested in
> `tests/test_transom.py` (96 cases), exercised end to end by
> `examples/transom_stern/run_example.py` (figures + CSV under `outputs/transom_stern/`).

## 1. What was missing and why it matters

The parametric hull in `tensorcad/hull.py` tapers to a **point at both ends**:
`b(x, ζ) = (B/2)·sin(πx/L)^0.7·√ζ`, so the aft body closes on the centreline. That is
a legitimate closed surface and a good self-verifying benchmark, but it is not a stern.
Every displacement vessel above `Fn ≈ 0.25` carries an **immersed transom** — a flat
plate where the shell is cut off — and that plate changes the wetted surface, the
pressure field and the wave-making that a downstream BEM/resistance code sees first.
The backlog entry was *"P0: 非对称船体 + 艉板"*; this module closes the **transom**
half of it (asymmetry between port and starboard remains open, see §6).

## 2. Geometry

Parameters (`TransomHullParameters`):

| symbol | field | meaning |
|---|---|---|
| `L` | `length` | **overall length, transom to stem (LOA)** |
| `B` | `beam` | maximum moulded breadth at the design waterline |
| `T` | `draft` | keel to waterline |
| `τ` | `transom_fraction` | station of the cut, as a fraction of the *parent* body, `0 ≤ τ < 0.5` |
| `p` | `breadth_exponent` | exponent of the sine taper, `0.25 ≤ p ≤ 3.0` |
| `nx`, `nz` | — | longitudinal stations / vertical levels |

The parent body is `L_parent = LOA/(1-τ)` long, and the mesh spans `x ∈ [0, LOA]`
with the **transom at `x = 0`** (the frame of `hull.py`, "x from the aft end").

Construction: one closed section ring per station (starboard waterline → keel → port
waterline, the keel point shared), lofted into quads, a waterline cap strip, and — for
`τ > 0` — a **planar fan** over the aft ring closing the stern. The section is convex
(a parabolic keel profile under a flat waterline), so the fan needs no new vertices and
no triangulation library. `τ = 0` skips the plate and reproduces the pointed hull
(same faces; vertices equal to 1e-12 — see the test).

Topology is invariant across `τ ∈ {0, 0.01, 0.05, 0.2, 0.35, 0.49}`: watertight, edge
manifold, `V - E + F = 2`, outward oriented, no zero-area faces.

## 3. Closed-form references, and how they are validated

Truncation turns the complete integral `I = ∫₀¹ sin(πu)^p du` into an incomplete one,
`J(τ) = ∫_τ¹ sin(πu)^p du`, which has no elementary closed form for general `p`. The
reference is kept independent of the mesh by splitting the integrand:

```
sin(πs)^p = (πs)^p · (sinc(πs))^p ,   (sinc(πs))^p = Σ_k c_k s^{2k}
```

The first factor is the only singularity and integrates exactly; the second is even and
therefore a power series in `s²` with factorial convergence. Term by term

```
∫₀^a s^m sin(πs)^p ds = π^p Σ_k c_k · a^(m+p+2k+1)/(m+p+2k+1)
```

(`sinc_power_coefficients`, 24 terms). The result is a reference accurate to
**round-off**, and it is checked three ways in `tests/test_transom.py`:

1. at `τ = 0` it must equal `taper_integral(p) = Γ((1+p)/2)/(√π Γ(1+p/2))` — an
   independent gamma-function derivation — to `1e-14` for `p ∈ [0.25, 3]`;
2. the first moment must satisfy `M(0)/J(0) = 1/2` exactly (the pointed body's LCB is
   midship), to `1e-14`;
3. for `τ > 0` it must agree with `tanh_sinh_quadrature` (the double-exponential rule
   already used for the hull's endpoint singularities) to `1e-13`.

Measured worst disagreement over the published ladder
(`outputs/transom_stern/reference_ladder.csv`): **2.2e-16** against the gamma form and **1.1e-15** against
tanh-sinh.

From `J(τ)` and `M(τ) = ∫_τ¹ u sin(πu)^p du`:

| quantity | exact value |
|---|---|
| displaced volume | `(2/3)·B·T·L_parent·J` |
| waterplane area | `B·L_parent·J` |
| LCB = LCF (fwd of transom) | `L_parent·(M/J - τ)` |
| VCB | `-0.4·T` (independent of `τ`, `p`) |
| `Cb` over the LOA | `(2/3)·J/(1-τ)` |
| `Cp` | `J/(1-τ)` |
| `Cwp` | `J/(1-τ)` |
| midship section area | `(2/3)·B·T` |
| immersed transom area | `(2/3)·B·T·sin(πτ)^p` |

## 4. Mesh integration against the references

`outputs/transom_stern/convergence.csv`, reproduced (`V` in m³, errors relative to the
analytic value):

| τ | 41×13 | 81×25 | 161×41 | 321×81 |
|---:|---:|---:|---:|---:|
| 0.00 | −0.880 % | −0.308 % | −0.136 % | −0.049 % |
| 0.05 | −0.813 % | −0.286 % | −0.129 % | −0.045 % |
| 0.20 | −0.785 % | −0.278 % | −0.126 % | −0.045 % |
| 0.35 | −0.770 % | −0.274 % | −0.125 % | −0.044 % |

Waterplane area and LCB converge faster (the truncation removes the singular aft
station): at 41×13 the worst waterplane error is −0.18 % and the worst LCB error is
−0.065 %; at 321×81 they are −0.005 % and −0.002 %. Convergence rate is ≈ 1.4–1.5 order in spacing,
not 2, because the `√ζ` keel profile has an infinite slope — the same signature the
pointed hull shows, and the reason the analytic reference (not self-refinement) is the
acceptance criterion. **All key quantities are inside 3 % at every grid tested, and
inside 0.2 % from 161×41 up.**

The transom plate itself is checked separately: the planar fan under-reads the
parabolic section (piecewise-linear chords inside a convex curve), −0.71 % at `nz = 13`,
−0.12 % at `nz = 41`, −0.043 % at `nz = 81` — monotone, and reported as a warning-level
check rather than silently ignored.

## 5. Design-space use

The immersed transom area ratio is `A_t/A_m = sin(πτ)^p`, and
`transom_fraction_for_area_ratio` inverts it in closed form, so a design loop can ask
for "a 10 % immersed transom" instead of guessing `τ`:

| `p` | `A_t/A_m` | `τ` | `Cb` | `Cp = Cwp` | LCB/LOA |
|---:|---:|---:|---:|---:|---:|
| 0.7 | 0.05 | 0.0044 | 0.4715 | 0.7073 | 0.4979 |
| 0.7 | 0.15 | 0.0212 | 0.4784 | 0.7177 | 0.4905 |
| 1.5 | 0.10 | 0.0691 | 0.3965 | 0.5947 | 0.4653 |
| 3.0 | 0.15 | 0.1783 | 0.3386 | 0.5079 | 0.3989 |

Two things to read off this table. First, with the default fine run (`p = 0.7`) a
*realistic* immersed transom (a few percent of the midship section, which is what
practice aims at so that the plate is only just wetted at the design condition) needs
`τ ≈ 0.005–0.02`, i.e. a cut of only 5–20 cm on a 10 m hull; the sine run is that
fine. Fuller bodies (`p` larger) reach the same area ratio with a visibly larger
transom — that is the trade the two parameters expose. Second, `Cb` grows with `τ` at
fixed LOA (the removed slice is thin, the LOA shrink is not), flattening past
`τ ≈ 0.2`; `p` moves `Cb` over the much wider range 0.31–0.48, so it is the primary
block-coefficient knob, not `τ`.

## 6. Limitations (honest list)

- **`Cp = Cwp` identically.** The family has a constant midship section and a
  waterplane equal to the section-area distribution scaled by `2/3`, so prismatic and
  waterplane coefficients coincide. Real hulls decouple them (parallel middle body,
  fore/aft asymmetry). Do not use this family to study `Cp`/`Cwp` trade-offs.
- **No port/starboard asymmetry**, no flare/tumblehome, no bulbous bow, no skeg or
  appendages; the other half of backlog item "非对称船体 + 艉板" stays open.
- **Transom is vertical and planar at one station.** Raked transoms, buttock flow-off
  and the actual immersed-transom drag behaviour are not modelled — this is geometry
  plus hydrostatics, no resistance.
- **`τ ≥ 0.5` is rejected by construction**: the cut would remove the maximum section
  and invalidate `Cm = 2/3`, on which every reference above depends.
- Convergence is first-order-dominated by the `√ζ` keel edge; a body with a flat keel
  or chine would converge faster. The default 41×13 grid is a 0.8 % displacement
  estimate, acceptable for trends, not for acceptance against a 3 % target.
- No stability, structure, vibration or rule checks. Wetted surface is available from
  the mesh (includes the plate) but has no closed form here — `hydrostatics_analytic`
  covers only the pointed hull.

## 7. Reproduce

```bash
python3 -m pytest tests/test_transom.py -q
PYTHONPATH=. python3 examples/transom_stern/run_example.py --output outputs/transom_stern
```

Artefacts: `reference_ladder.csv` (series vs Γ vs tanh-sinh), `convergence.csv`
(4 truncations × 4 grids), `design_table.csv`, and `convergence.png`,
`geometry.png`, `design_space.png`.