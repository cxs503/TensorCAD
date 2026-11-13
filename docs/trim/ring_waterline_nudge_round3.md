# Ring waterline, round 3: the nudge both ways, and one error type

**Mission**: 17215 (reopen of TensorCAD#349, "管到底", round 3)
**Module changed**: `tensorcad/trim.py` (`hydrostatic_particulars`)
**Tests changed**: `tests/test_trim.py`
**Numbers**: `artifacts/trim_round3/numbers.json` + `artifacts/trim_round3/sweep.json`
**Probe**: `artifacts/trim_round3/probe_round3.py` (every number and the figure come from it)
**Figure**: `artifacts/trim_round3/trim_round3_nudge.png`

## 1. Case

`generate_hull(HullParameters(length=10, beam=2, draft=0.8))` puts vertex rings on even
z levels.  The waterline the trim suite uses, `z = -0.25 T = -0.200 m`, lands on a ring with
**80 vertices**.  `hydrostatic_particulars` has a fallback for a degenerate plane cut: it
retries `1e-6` of the mesh scale away from the requested waterline and reports
`waterline_nudged=True`.

Rounds 1 and 2 established that the shipped geometries cut the ring *without* the fallback
(`PLANE_TOLERANCE` in `stability._plane_segments`) and covered the `True` branch with a
monkeypatched refusal.  They left two open points, both named in
`docs/trim/ring_waterline_coverage.md §5`:

* the fallback retried **downward only**, so a degenerate plane that *bounds* the body —
  the keel ring — could never be rescued: moving down leaves the hull and the retry fails;
* a double failure leaked `tensorcad.stability.StabilityError` out of a function whose
  documented contract is `TrimError`.

## 2. What was done

In `hydrostatic_particulars`:

1. The retry tries `waterline_z - offset` **and then** `waterline_z + offset`, downward
   first so the documented direction and the existing round-2 test are unchanged.  A ring
   that bounds the body is only cut from the inside, so one direction is insufficient.
2. If both retries refuse, the function raises `TrimError` chained to the original
   `StabilityError` (`raise ... from first_error`), instead of propagating the stability
   module's error.

`tests/test_trim.py`:

* `test_a_failed_nudge_does_not_invent_a_waterplane` now asserts `TrimError` **and** that
  the cause is still a `StabilityError` — the refusal does not lose the plane cut's reason;
* `test_the_nudge_tries_above_the_waterline_too` refuses every plane at or below the
  requested waterline (a lower-boundary ring) and checks that the upward retry is taken:
  three planes are tried (`-0.200000, -0.200010, -0.199990`), symmetric about the request,
  `waterline_nudged=True`, and the area stays within the documented O(1e-6) move.

No geometry, hull generator or `stability.py` change: the tolerance that makes ring cuts
ordinary is untouched.

## 3. Numbers (from `artifacts/trim_round3/numbers.json`)

| quantity | value | threshold |
|---|---|---|
| vertices on the ring at z = -0.200 m | 80 | >= 4 |
| waterplane area at the ring | 12.177763 m^2 | > 0 |
| TPC at the ring | 0.124822 t/cm | > 0 |
| documented nudge offset | 1.000e-05 m | <= 1e-4 m |
| one-sided limit, z - offset | -8.579e-06 | - |
| one-sided limit, z + offset | +8.114e-06 | - |
| second difference across the ring | 4.648e-07 | <= 1e-6 |
| area rel. dev. of the upward retry | +8.114e-06 | <= 2e-5 |

Before/after on the two blind spots (the "legacy" column is the pre-fix function,
re-implemented in the probe as the baseline):

| scenario | legacy (downward only) | fixed (both sides) |
|---|---|---|
| both plane cuts refused | `StabilityError` | `TrimError` (cause `StabilityError`) |
| ring at the body's lower boundary | `StabilityError` | answered, `waterline_nudged=true`, abs dev = 8.11e-06 |

Continuity through the ring is measured as the **second difference** across it,
`abs(A(z-d) + A(z+d) - 2 A(z)) / A(z)` over the ±1e-5 m retry offsets: a kink would break
the straight line the physical area follows.  It is 4.648e-07, i.e. the ring is an ordinary
point of the curve to well below the waterplane's discretisation error.

## 4. Error judgement

| key quantity | deviation | is it < 3%? |
|---|---|---|
| waterplane area of the upward retry vs the un-nudged cut | 8.11e-06 | **yes** |
| one-sided limits about the ring vs the ring cut | <= 8.58e-06 | **yes** |
| second difference across the ring (continuity) | 4.65e-07 | **yes** |

There is no closed-form waterplane for a hull in this package (the box is the only analytic
geometry and it has no interior rings), so the reference is the smooth continuation from
both sides plus the documented 1e-6-of-scale retry; both are five orders of magnitude below
3%.

## 5. Notes and limits

- The fallback is still unreachable through the shipped geometries at any waterline the
  function accepts: probing 101 levels from keel to deck, `waterplane_moments` fails only at
  the keel ring itself, where `hydrostatic_particulars` refuses earlier (`draft <= 0`).
  Both branches are therefore driven with an injected cutter; the tests say so in their
  docstrings rather than pretending a mesh triggers them.
- Which side was taken is *not* recorded in the record — `waterline_nudged` says "moved",
  the docstring says "within 1e-6 of the scale".  Adding the signed offset would change the
  persisted schema (`tensorcad.hydrostatic_particulars/1`), so it is left to a schema
  revision.
- `TrimError` is a `ValueError`, as is `StabilityError`, so callers that catch `ValueError`
  (or callers of the CLI) see no change in kind, only in the concrete type.
- `mypy tensorcad/trim.py` reports pre-existing errors elsewhere in the package (68 in 16
  files on main); `ruff` is clean on both changed files.
- Full suite in this container: **3168 passed, 43 failed, 26 skipped, 1 xfailed** (17 m 30 s).
  None of the failures touches trim; they are environment/pre-existing on unmodified main:
  11 `test_bonjean.py` (`TypeError: section_properties() takes 2 positional arguments but 3
  were given` — the `tensorcad/__init__` export mismatch `mypy` also reports), 4
  `test_examples_bootstrap.py` (detangle / elliptic_section / fillet_independent_verification
  / mesh_repair, reproduced on unmodified main), missing-artifact suites
  (`test_geometry_audit_m5479.py`, `test_oss_inventory_cad.py` ×10 —
  `FileNotFoundError`), `test_planing_wetted_length_verification.py` ×1, and the remaining
  failures come from the optional `meshio` dependency being absent here.
  A cluster full-suite job (cj1458) was submitted for the loong64 nodes and cancelled once
  the local run answered; `pytest tests/test_trim.py tests/test_trim_guards.py` is 41/41
  green with `tensorcad/trim.py` coverage 92.18 % (gate 90 %).