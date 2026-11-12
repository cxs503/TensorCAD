# Port/starboard asymmetric transom stern

`tensorcad.transom` truncates the pointed hull of `tensorcad.hull` with a flat
transom plate, but both sides of that body are mirror images. `docs/transom_stern.md` §6
says so and names the half of the backlog item *"asymmetric hull + transom"* that stayed
open: the **transverse** one. `tensorcad.hull`'s `fore_exponent`/`aft_exponent` asymmetry is
longitudinal, and `tensorcad.transom_asymmetric` is the transverse case — a hull whose
starboard side is wider than its port side, still closed by a transom plate.

## 1. The family

| symbol | field | meaning |
|---|---|---|
| `L` | `length` | overall length, transom to stem (LOA) |
| `B` | `beam` | **overall** waterline breadth, starboard half + port half |
| `T` | `draft` | keel to waterline |
| `tau` | `transom_fraction` | station of the cut as a fraction of the parent body, `0 < tau < 0.5` |
| `p` | `breadth_exponent` | exponent of the sine taper, `0.25 <= p <= 3.0` |
| `s` | `skew` | breadth skew, `|s| <= 0.75` |
| `nx`, `nz` | — | longitudinal stations / vertical levels |

The section is the transom family's separable profile, with one half breadth per side:

```
y_s(x,z) = +(B/2)(1+s) f(x) sqrt(1 + z/T)      (starboard)
y_p(x,z) = -(B/2)(1-s) f(x) sqrt(1 + z/T)      (port)
f(x)     = sin(pi (tau + x/L_parent))**p,      L_parent = LOA/(1-tau)
```

`y_s + |y_p| = B f(x)` at every station, so **the asymmetry slides the body, it does not
widen it**: the bounding box breadth is `B` for every `s` (checked against the declared
dimensions by `check_dimensions` and by a dedicated `breadth_span_preserved` check).
`s = 0` reproduces `tensorcad.transom` bit for bit — same vertices, same faces, same
analytic reference (`tests/test_transom_asymmetric.py::test_zero_skew_reproduces_the_symmetric_family_exactly`).

## 2. Closed forms

With `J1 = int_tau^1 f du`, `M = int_tau^1 u f du`, and `Jk = int_tau^1 f**k du`
(`Jk = truncated_taper_integral(tau, k*p)` — the same incomplete-integral series the transom
module validates against the gamma-function closed form at `tau = 0`):

| quantity | closed form |
|---|---|
| displacement volume | `(2/3) B T L_parent J1` |
| waterplane area | `B L_parent J1` |
| LCB / LCF (fwd of transom) | `L_parent (M/J1 - tau)` |
| VCB | `-2/5 T` |
| `Cb` / `Cwp` | `(2/3) J1/(1-tau)` / `J1/(1-tau)` |
| `y_B` (transverse centre of buoyancy) | `(3/8) s B J2/J1` |
| `y_F` (transverse centre of flotation) | `(1/2) s B J2/J1` |
| `I_T` about the centreline | `(B_s**3 + B_p**3)/24 L_parent J3`, `B_s,p = B(1 ± s)` |
| `I_T` about the floating centre | `I_T - A_w y_F**2` |

The derivation of `y_B` is one line of section moments: the horizontal slice at level `z`
has length `(y_s + y_p)` and centroid `(y_s - y_p)/2`, so its `y` first moment is
`(y_s**2 - y_p**2)/2`; integrating `(1 + z/T)` over the draft gives `T/2`, leaving
`int y dA = T (B_s**2 - B_p**2) f(x)**2 / 16` per station — the `f**2` is why `J2` appears.

Every quantity that is an integral of the *sum* of the half breadths is **skew-free to
round-off**: volume, waterplane area, LCB, LCF, VCB, `Cb`, `Cp`, `Cwp`, maximum section
area and immersed transom area (metrics: |Δ| = 0 exactly, not "small", in
`tests/test_transom_asymmetric.py::test_skew_cannot_change_the_volume_or_its_longitudinal_centres`).
`y_F/y_B = 4/3` identically, for any `tau`, `p` and grid — a shape invariant, not a fit.

## 3. Results

The reference body is `L = 10 m`, `B = 2 m`, `T = 0.8 m`, `tau = 0.02`, `p = 0.7`.
Tables are generated, not typed: `docs/transom_asymmetric/tables.md` (written by
`examples/transom_asymmetric/run_example.py`), numbers in
`artifacts/transom_asymmetric/results.json`, sha256 of every artifact in `manifest.json`.

### 3.1 Mesh against the closed forms (`skew = 0.4`)

| nx | nz | V [m³] | V error | Aw error | LCB error | y_B [m] | y_B error |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 41 | 13 | 7.584546 | -0.824 % | -0.120 % | -0.020 % | 0.244181 | +0.042 % |
| 81 | 25 | 7.625484 | -0.289 % | -0.035 % | -0.008 % | 0.243731 | +0.019 % |
| 161 | 41 | 7.637668 | -0.129 % | -0.010 % | -0.003 % | 0.243548 | +0.010 % |
| 321 | 81 | 7.644081 | -0.046 % | -0.003 % | -0.001 % | 0.243427 | +0.004 % |

Observed orders (log-log, measured): volume 1.51 / 1.16 / 1.51, waterplane 1.78 / 1.77 / 1.76,
`y_B` 1.12 / 0.93 / 1.33. **The rate is first-order, not second**, and that is not a defect
of the asymmetry: the section's `sqrt` keel profile has an infinite slope at the keel, so the
piecewise-linear chord under-reads the section by `O(h)`. `docs/transom_stern.md` reports the
same behaviour for the symmetric family. The asymmetry adds no new error term: the `y` moments
converge at the same first order.

### 3.2 The skew sweep (closed forms; Simpson cross-check in the JSON)

| skew | B_s/2 [m] | B_p/2 [m] | V [m³] | Aw [m²] | y_B [m] | y_F [m] | y_F/y_B | I_T cl [m⁴] | I_T cg [m⁴] |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| -0.60 | 0.4000 | 1.6000 | 7.647568 | 14.339189 | -0.365020 | -0.486693 | 1.333333333 | 6.941675 | 3.545150 |
| -0.30 | 0.7000 | 1.3000 | 7.647568 | 14.339189 | -0.182510 | -0.243346 | 1.333333333 | 4.238426 | 3.389295 |
| +0.00 | 1.0000 | 1.0000 | 7.647568 | 14.339189 | +0.000000 | +0.000000 | &mdash; | 3.337344 | 3.337344 |
| +0.30 | 1.3000 | 0.7000 | 7.647568 | 14.339189 | +0.182510 | +0.243346 | 1.333333333 | 4.238426 | 3.389295 |
| +0.60 | 1.6000 | 0.4000 | 7.647568 | 14.339189 | +0.365020 | +0.486693 | 1.333333333 | 6.941675 | 3.545150 |

Read the two inertias together: at `skew = ±0.6` the centreline inertia is 2.08× the
symmetric value, but the *centroidal* inertia is only 1.06× — most of the extra `I_T` is the
`A_w y_F**2` shift of an offset waterplane, not a stiffer waterplane. Using `I_T` about the
centreline as `BM_T` for an asymmetric hull would overstate the metacentric height by 96 % at
this skew. That is the trap this pair of numbers exists to expose.

## 4. Figures

* `artifacts/transom_asymmetric/convergence_and_skew.png` — log-log convergence of volume
  and `y_B` against the closed forms; `y_B`, `y_F` linear in the skew and the Simpson route
  landing on the same line.
* `artifacts/transom_asymmetric/geometry.png` — the waterline plan of a `skew = +0.5` body
  (starboard half wider, same span) and the transom/midship sections.

## 5. Error judgement against the references

| quantity | reference | grid | error | < 3 %? |
|---|---|---|---:|---|
| displacement volume | closed form (series `J1`) | 41×13 | -0.824 % | **yes** |
| displacement volume | closed form | 321×81 | -0.046 % | **yes** |
| waterplane area | closed form | 41×13 | -0.120 % | **yes** |
| LCB | closed form (`M/J1`) | 41×13 | -0.020 % | **yes** |
| `y_B`, `skew = 0.4` | closed form (`3/8 s B J2/J1`) | 41×13 | +0.042 % | **yes** |
| `y_B`, `skew = 0.4` | closed form | 321×81 | +0.004 % | **yes** |
| `I_T` (centreline) | closed form (`J3`) | 161×41 | < 0.5 % (test) | **yes** |
| volume, closed form | independent composite Simpson (4000 panels) | — | 8.4e-8 worst | **yes** |
| `y_B`, closed form | independent Simpson route | — | 2e-8 worst (test) | **yes** |
| `y_F/y_B` | exact `4/3` | any | 0 (to 1e-14) | **yes** |
| volume / `Aw` / LCB vs `skew` | invariance (must not move) | — | 0 exactly | **yes** |

No quantity above is a comparison against a *measured* ship, because none exists for this
family: the references are closed forms and an independent quadrature of the same analytic
surface. That is a verification of the implementation, not of the geometry family's realism.

## 6. Points to watch

1. **Upright only.** Nothing here is a stability or resistance claim. An asymmetric waterplane
   has a heel angle at which the body *prefers* to float; this module does not compute it, and
   `GZ` from `tensorcad.gz` is not valid for a body whose upright attitude is not its
   equilibrium.
2. **One taper, two breadths.** Port and starboard share `p`, the deadrise and the profile.
   A hull with a sponson or an asymmetric skeg is not in this family.
3. **`tau > 0` is required.** With no transom plate the degenerate aft ring would only add
   zero-area faces to check, and the pointed family already exists in `tensorcad.hull`.
4. **`I_T` about the centreline is not `BM_T`.** Subtract `A_w y_F**2` (see §3.2).
5. **The offsets table is one number per station.** `build_asymmetric_mesh` returns the
   *overall* breadth `y_s - y_p`; the split is `AsymmetricTransomParameters.half_breadths()`.
   A lines plan of a real asymmetric hull needs both halves, so a consumer that wants them
   must ask for them — the offsets tuple is deliberately the same shape as the symmetric
   family's.
6. Converged at 41×13 for trends (−0.8 % displacement), needs 161×41 or finer for anything
   that is compared to a 3 % target.

## 7. Reproducing

```bash
python3 examples/transom_asymmetric/run_example.py         # artifacts + generated tables
python3 -m pytest tests/test_transom_asymmetric.py -q      # 50 cases
```

The script exits non-zero if any of its eight published claims stops holding
(`claims` block in `artifacts/transom_asymmetric/results.json`).