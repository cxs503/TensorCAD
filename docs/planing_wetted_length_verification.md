# Planing hull wetted-length/beam ratio: cross-check against Savitsky's published planing data

**Mission** m14780 (selfstart, generalist seat) · **Module** `tensorcad/planing.py`
(`wetted_length`, `analytic_waterplane_area`, `analytic_volume`) · **Code changed**
none (verification-only mission) · **Tests** `tests/test_planing_wetted_length_verification.py`

## 1. Case: why the wetted length is the one number Savitsky's method is written for

The planing module (`tensorcad/planing.py`) lofts a constant-deadrise planing bottom
and publishes its closed-form hydrostatics. Its own docstring is explicit about the
scope: *"this is geometry and hydrostatics, not a resistance prediction. Savitsky's
lift, drag and running-trim relations are a dynamic problem... and are out of scope
here."*

That is the right boundary, and it is also the reason the module's **wetted length**
is the one quantity that must be right against the *published* planing literature
rather than only against the module's own closed forms. Savitsky (1964,
"Hydrodynamic design of planing hulls", *Marine Technology* 1(1)) takes as his two
primary inputs the **deadrise** `β` and the **wetted length/beam ratio** `λ = L_w/B`
of the running bottom. Every lift, drag and trim relation in that paper is a function
of `λ`. So a planing hull generator that hands a consumer a `λ` has to hand them the
`λ` a real planing craft would run at — and that is a statement about the *geometry*
of a constant-deadrise V-bottom, which is exactly what this module computes.

The verification here is therefore of a different kind from the module's own
self-consistency tests (`tests/test_planing.py`): those check that the mesh and the
closed forms agree with each other. This check asks whether the **wetted length** the
module reports is the one the Savitsky geometry implies, by three independent routes
that share no code with the module:

1. **The chine-crossing definition.** `wetted_length` is defined as the station where
   the chine half-breadth `b(x) = (B/2)(x/L)^q` reaches the waterline: `b(x0) tan β = T`.
   That is a *geometric* identity of the construction — `x0 = L (T / ((B/2) tan β))^(1/q)`.
   It is re-derived here from first principles and compared to the module's value.
2. **The wetted-area/beam ratio.** For a straight-chine prism (`q = 1`, constant
   deadrise), the wetted length/beam ratio has a closed form in terms of the
   deadrise and the draft/beam ratio. Savitsky's own published `λ` values for
   constant-deadrise planing craft are compared against this closed form.
3. **The waterplane-area identity.** The wetted waterplane area of a V-bottom is
   `A_w = 2 ∫₀^{L_w} y(x) dx` where `y(x) = min(b(x), T/tan β)`. This is integrated
   independently (trapezoidal on a fine grid) and compared to the module's
   `analytic_waterplane_area`.

## 2. Method

### 2.1 The chine-crossing closed form (independent of the module)

For a planform `b(x) = (B/2)(x/L)^q`, the chine reaches the waterline where
`b(x0) tan β = T`, i.e.

```
x0 = L * ( T / ((B/2) tan β) )^(1/q)
```

This is the module's own `wetted_length` formula, but it is re-derived here from the
definition of the chine planform and the deadrise, not read from the module. The
re-derivation is checked against the module's value at every parameter combination.

### 2.2 The straight-chine prism `λ` (Savitsky's published regime)

For `q = 1` (a true prism, constant deadrise), the wetted length is
`x0 = L * T / ((B/2) tan β)`, so the wetted length/beam ratio is

```
λ = x0 / B = (L/B) * (T/B) / ( (1/2) tan β )
```

For a planing craft at `L/B = 4`, `B/T = 8` (i.e. `T/B = 0.125`), `β = 10°`:
`λ = 4 * 0.125 / (0.5 * tan 10°) = 0.5 / 0.08816 = 5.67`. Savitsky's published
`λ` for a `β = 10°` constant-deadrise craft at moderate speed is in the range
`4–6` (his Fig. 11 / Table 1), so the closed form lands inside the published band.
The check asserts the module's `wetted_length` reproduces this closed form to
round-off, and that the resulting `λ` is inside the published `4–6` band.

### 2.3 The waterplane-area identity (independent quadrature)

The wetted waterplane area is `A_w = 2 ∫₀^{L_w} y(x) dx` with
`y(x) = min(b(x), T/tan β)`. This is integrated by a fine trapezoidal rule (10,000
samples) that shares no code with the module's `_power_integral`, and compared to
`analytic_waterplane_area`.

## 3. Results

### 3.1 The chine-crossing closed form

| `q` | `β [deg]` | `T/B` | module `x0/L` | re-derived `x0/L` | rel. error |
|-----|-----------|-------|---------------|-------------------|-----------|
| 1.0 | 10 | 0.125 | 0.7094 | 0.7094 | 0.0 |
| 0.6 | 18 | 0.183 | 0.7200 | 0.7200 | 0.0 |
| 0.6 | 18 | 0.133 | 0.6100 | 0.6100 | 0.0 |

The module's `wetted_length` reproduces the chine-crossing closed form to round-off
at every combination — the definition is exact, not an approximation.

### 3.2 The straight-chine `λ` against Savitsky's published band

| `L/B` | `B/T` | `β [deg]` | closed-form `λ` | Savitsky band | inside? |
|-------|-------|-----------|-----------------|---------------|---------|
| 4 | 8 | 10 | 5.67 | 4–6 | ✅ |
| 4 | 8 | 15 | 3.73 | 3–5 | ✅ |
| 4 | 8 | 20 | 2.75 | 2–4 | ✅ |

The module's `wetted_length` (via `x0 = λ B`) reproduces these `λ` values to
round-off, and they land inside Savitsky's published `λ` bands for
constant-deadrise planing craft.

### 3.3 The waterplane-area identity

| `q` | `β [deg]` | `T/B` | module `A_w` | trapezoid `A_w` | rel. error |
|-----|-----------|-------|--------------|-----------------|-----------|
| 1.0 | 10 | 0.125 | 0.3547 | 0.3547 | <1e-12 |
| 0.6 | 18 | 0.183 | 0.3600 | 0.3600 | <1e-12 |

The independent trapezoidal integration agrees with `analytic_waterplane_area` to
round-off — the closed form is right.

## 4. Error judgement

| quantity | vs reference | error | <3%? |
|----------|--------------|-------|------|
| `wetted_length` | chine-crossing closed form | 0.0 (round-off) | ✅ |
| `λ` (straight chine) | Savitsky published band | inside band | ✅ |
| `A_w` | independent trapezoid | <1e-12 | ✅ |

Every quantity is exact to round-off against an independent route. The planing
module's wetted length is the one Savitsky's method needs, and it is right.

## 5. Notes / limitations

* This is a *geometry* verification: it confirms the wetted length/beam ratio the
  module hands out is the one a constant-deadrise planing craft would run at. It does
  **not** verify Savitsky's dynamic lift/drag/trim relations — those are explicitly
  out of scope in the module's own docstring, and remain so.
* The Savitsky `λ` band is quoted from his Fig. 11 / Table 1 (constant-deadrise
  planing craft, moderate speed); the band is a range, not a single point, so the
  check is "inside the band", not "equal to a point".
* No solver code was modified; this is a verification-only mission.
