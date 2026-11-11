# Resistance: three friction lines, a Prohaska form factor, and the trip to full scale

TensorCAD answers questions about a body -- how much water it displaces, where the
waterplane is, what the righting arm is at 20 degrees. This page documents the module
that answers the question the body exists to be asked: **at this speed, how hard does it
pull?** -- and its mirror image, addressed to a towing tank rather than a hull: *I
measured this on a 10 m model; what will the ship pay?*

Everything on this page is measured by `python examples/resistance/run_example.py`,
which writes `artifacts/resistance/results.json`, `plot.png` and `report.md` and exits
non-zero if any claim stops holding. `tests/test_resistance.py` pins the same claims (60
cases). Numbers below are read from that artifact, not copied by hand.

## 1. What the module is, and what it deliberately is not

`tensorcad/resistance.py` implements the **1966 ITTC method, first stage** (Prohaska):
a friction line, a form factor `(1 + k)` fitted from the model's own resistance curve,
and the **1978 ITTC** two-dimensional scaling of a model test to the ship. Units are SI.

The residual (wave + separation) resistance is **not predicted here**. That number is a
property of the free surface; this package's position is that it comes from a model
test, a panel code (`TensorBEM`) or a RANS solver, not from an algebraic guess. The API
therefore has *no default* for it, and no default for the form factor either: the
defaults (`1 + k = 1`, `CR = 0`) read as **"not measured"**, and the returned record
carries both so a consumer can see whether a resistance was predicted or merely bounded.

| Public name | What it does |
|---|---|
| `ittc57_friction_coefficient(Re)` | `0.075 / (log10 Re - 2)^2`, the 1957 ITTC line |
| `schoenherr_friction_coefficient(Re)` | `0.242/sqrt(CF) = log10(Re CF)`, solved by bisection -- the *independent route* |
| `prandtl_schlichting_friction_coefficient(Re)` | `0.455 / (log10 Re)^2.58`, a third estimate |
| `friction_coefficient(Re, line=...)` | the same, by name (`FRICTION_LINES`) |
| `prohaska_fit(Fn, CT, *, length=... \| friction_coefficients=...)` | `1 + k` and the `Fn^4` slope, with residuals and standard errors |
| `bare_hull_resistance(...)` / `mesh_bare_hull_resistance(mesh, ...)` | `R_F`, `(1+k)R_F`, `R_R`, `R_T`, `P_E` for one body at one speed |
| `ship_from_model_test(...)` | Froude scaling + `CT_s = (1+k)CF_s + CR + CA` |

Explicitly absent, each with the number it would move: wave resistance (no Michell
integral, no Rankine source method -- `TensorBEM`), a Holtrop-Mennen style regression
over ship particulars (a shipyard correlation with its own validity box, not a geometry
operation), a propulsive-coefficient or propeller model, a roughness allowance model
(`CA` is a caller input), and seakeeping added resistance.

## 2. The friction lines, and the honest distance between them

| Re | `CF` ITTC-57 | `CF` Schoenherr | `CF` Prandtl-Schlichting | ITTC-57 vs Schoenherr | ITTC-57 vs P-S |
|---:|---:|---:|---:|---:|---:|
| 1e4 | 0.018750 | 0.013074 | 0.012726 | +43.41 % | +47.33 % |
| 1e5 | 0.008333 | 0.007179 | 0.007156 | +16.07 % | +16.45 % |
| 1e6 | 0.004687 | 0.004409 | 0.004471 | +6.31 % | +4.85 % |
| 1e7 | 0.003000 | 0.002934 | 0.003004 | +2.24 % | −0.12 % |
| 1e8 | 0.002083 | 0.002072 | 0.002128 | +0.55 % | −2.11 % |
| 1e9 | 0.001531 | 0.001531 | 0.001571 | −0.02 % | −2.55 % |

Three facts, all asserted in the tests rather than left as prose:

* **In the ship range the agreement is what the 1957 line was fitted for**: worst
  |Δ| against Schoenherr for `Re >= 1e7` is **2.24 %**, and against Prandtl-Schlichting
  **2.55 %**. Both sit under the 3 % accuracy bar this repository applies to a key
  quantity, so treating the 1957 line as a proxy for Schoenherr is defensible where it
  is actually used (models and ships, `Re >= 1e7`).
* **Below `1e6` the lines genuinely part** -- 16.1 % at `1e5`, 43.4 % at `1e4`. This is
  a property of the published lines, not an implementation defect, and the test asserts
  the divergence *exists* (0.15 < Δ < 0.17 at `1e5`): a future change that made them
  agree there would be a different model, not a bug fix.
* The Schoenherr solve is checked against **its own implicit equation**, not against a
  table: the worst residual over the grid is **9.5e-14**, which is a statement about the
  bisection, independent of whether the line itself is right.

## 3. The wetted area comes off the mesh -- and the deck is not wetted

Wigley reference hull (`L = 10 m`, `B = 1 m`, `T = 0.625 m`), three grids. Shell area is
`mesh.surface_area`; the deck patch is the triangles lying in the `z = 0` waterline
plane; the reference is the analytic flank integral
`tensorcad.wigley.analytic_wetted_surface_area`, which quadratures the analytic surface
and never touches the mesh.

| nx x nz | shell [m²] | deck [m²] | deck share | wetted [m²] | analytic [m²] | error |
|---:|---:|---:|---:|---:|---:|---:|
| 41 x 13 | 21.538575 | 6.662500 | 30.93 % | 14.876075 | 14.879063 | −0.0201 % |
| 81 x 23 | 21.543616 | 6.665625 | 30.94 % | 14.877991 | 14.879063 | −0.0060 % |
| 161 x 43 | 21.545201 | 6.666406 | 30.94 % | 14.878795 | 14.879063 | −0.0016 % |

The wetted area converges to the analytic flank at second order (−0.020 % → −0.0016 %),
which is the geometry half of the resistance statement: the `S` that multiplies `0.5
rho V^2 CF` is a measured area, not a coefficient someone remembered.

**The trap this table exists to document**: the closed hull's *shell* is 45 % larger
than its wetted surface, because the families here close the body at the design
waterline and that cap is a real deck patch. `mesh_bare_hull_resistance` uses the shell
as given and says so in its docstring; a hull mesh capped at the waterline needs no
correction, and a freeboarded/decked mesh must have the patch subtracted (or
`wetted_surface` passed explicitly). Using 21.54 m² where 14.88 m² belongs overstates
friction by 44.8 %.

## 4. Prohaska: the form factor recovered from the curve

Prohaska's transform, `CT/CF = (1 + k) + y Fn^4/CF`, is a straight line in
`x = Fn^4/CF`. Synthetic "model test" totals were generated from an **injected**
`1 + k = 1.15` and `y = 3.0e-3` at `L = 10 m`, `Fn = 0.10 ... 0.24` (eight points):

| quantity | injected | recovered | relative error |
|---|---:|---:|---:|
| `1 + k` | 1.150000 | 1.150000000000 | 0.0 |
| `y` (wave coefficient) | 3.0000e-3 | 2.9999999999999e-3 | 3.1e-16 |
| max relative residual of the line | -- | 1.9e-16 | -- |

An exact recovery is the correct result *and* the only meaningful test of the method's
arithmetic: the transform and the least-squares fit are inverses of the generator, so
anything above round-off is a coding error. The tests also pin the negative: a 20 %
different injection (`1.05`) is recovered as `1.05`, i.e. the fit is not "close enough"
in a way that would swallow a real error; and scattered data (1 % alternating noise)
raises `max_relative_residual` above `1e-3` while reporting the standard errors of
intercept and slope instead of pretending the line is exact.

![plot](plot.png)

Left: the three flat-plate lines over `Re = 1e4 ... 1e9` (log-log). Right: the Prohaska
plot for the synthetic model test; the fitted intercept is the form factor.

## 5. Model to ship: a 10 m Wigley carried to 100 m

Model `L = 10 m` at `Fn = 0.30` (`V = 2.971 m/s`, `Re = 2.498e7`), ship `L = 100 m`,
geometrically similar wetted surface, `1 + k = 1.15`, `CA = 0`, sea water both sides.

| quantity | model | ship |
|---|---:|---:|
| `Re` | 2.4982e+07 | 7.9000e+08 |
| `CF` (ITTC-57) | 2.5743e-03 | 1.5764e-03 |
| `(1 + k) CF` | 2.9604e-03 | 1.8128e-03 |
| `CR` (from the model total) | 2.43e-5 | 2.43e-5 |
| `CT` | 2.9847e-03 | 1.8371e-03 |
| wetted surface | 14.879 m² | 1487.9 m² |
| `R_T` | 201.06 N | 123 753 N |
| `P_E = R_T V` | -- | 1.163 MW |

Three consequences, each asserted:

* **Speed scales as `sqrt(lambda)`**: `V_s / V_m = 3.1623`, and the two Froude numbers
  agree to 0.0 absolute -- Froude scaling preserves the wave pattern and nothing else.
* **`CT` falls by a factor 1.62** between model and ship (`scale_effect = 1.6247`).
  That width is the whole reason a ship is not `lambda^3` times its model: the naive
  ratio is `1000`, the Froude-scaled resistance ratio is **615.5**.
* **`CR` is recovered from the model total** by `CR = CT_m - (1 + k) CF_m`, and putting
  it back reproduces `CT_m` to 1e-16 -- the round trip is asserted in the tests, so the
  inverse of the scaling cannot silently drift. The same ship resistance computed by a
  second, direct route (build the ship body with `residual_coefficient = CR` from
  scratch) agrees to **0.0 relative**.

Effective power and total resistance over the model's speed range (same `(1+k)`, `CR`
scaled as `y Fn^4`):

| V [m/s] | Fn | R_T [N] | P_E [W] |
|---:|---:|---:|---:|
| 1.485 | 0.150 | 55.94 | 83.1 |
| 2.228 | 0.225 | 117.84 | 262.6 |
| 2.971 | 0.300 | 201.06 | 597.3 |
| 3.714 | 0.375 | 306.94 | 1139.9 |
| 4.456 | 0.450 | 439.44 | 1958.3 |

Monotone in speed, as `R_T ~ V^2 CF` with `CF` decreasing slowly in `Re`; the test
asserts the monotonicity rather than the individual values (the values are the
artifact's, not the contract's).

## 6. Error judgement against the reference values

| key quantity | reference | measured | error | < 3 %? |
|---|---|---:|---:|---|
| `CF` ITTC-57 at `Re >= 1e7` | Schoenherr line (independent, implicit) | worst 2.24 % | +2.24 % | **yes** |
| `CF` ITTC-57 at `Re >= 1e7` | Prandtl-Schlichting line | worst 2.55 % | +2.55 % | **yes** |
| `CF` ITTC-57 at `Re = 1e5` | Schoenherr line | 16.07 % | -- | **no** (stated, not hidden; the lines are not interchangeable at model-scale-low `Re`) |
| wetted area, 161 x 43 | analytic flank integral (quadrature of the analytic surface) | −0.0016 % | 0.0016 % | **yes** |
| wetted area, 41 x 13 | same | −0.0201 % | 0.020 % | **yes** |
| `1 + k` from Prohaska | injected `1.15` | 1.150000000000 | 0.0 | **yes** |
| Schoenherr root | its own implicit equation | residual 9.5e-14 | -- | **yes** |
| ship `R_T`, model-to-ship route | direct construction with the same `CR` | 0.0 relative | 0.0 | **yes** |

The third row is the one that matters for honesty: the module's friction line does
**not** agree with the other two below `Re ~ 1e6`, and no amount of tabulating changes
that. What is claimed is narrower and checkable -- *in the ship range, where the line is
used, the three agree to under 3 %*.

## 7. Points to watch

1. **Shell is not wetted area.** 30.9 % of a closed Wigley's shell is the waterline deck
   patch; using it inflates friction 44.8 %. `mesh_bare_hull_resistance` takes the mesh
   as given on purpose, and the docstring says which meshes need no correction.
2. **`CR` can come out negative** on a well-formed hull at low Froude number (the model
   measures less than the flat plate the method subtracts). It is reported, not clipped
   -- but a negative value on a *ship* prediction usually means the form factor or the
   water temperature is wrong, so it is worth reading before it is trusted.
3. **`CA` is not invented.** It is a caller input and is reported separately from `CR`,
   because a wrong allowance moves the answer by more than a wrong form factor does.
4. **Prohaska needs at least three points**; with two the line fits exactly and reports a
   zero residual, which is not evidence. The function refuses two.
5. **The transform assumes `R_W ~ Fn^4`**, so the fitted `1 + k` is only as good as that
   assumption holds -- conventionally `Fn` below about 0.4, and only on a hull whose
   wave resistance has not grown a hump.
6. The fluid constants (`SEA_WATER_*`, `FRESH_WATER_*`) are conventional 15 °C inputs and
   are not verified against an ITTC document in this repository. Every claim above is a
   ratio of like quantities and therefore independent of them.

## 8. Reproducing

```bash
python examples/resistance/run_example.py          # writes artifacts/resistance/
python -m pytest tests/test_resistance.py -q       # 60 cases
```

`artifacts/resistance/plot.png` is the figure above; `results.json` carries every number
in this page, including the six-row friction table, the three wetted-area grids and the
full prediction record of both bodies.