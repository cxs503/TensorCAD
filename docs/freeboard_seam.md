# The freeboard seam: a waterline that contains mesh edges

*Closed verification debt from `docs/wigley_gz_verification_gap.md` (opened by m6594,
closed here). Two defects, one root cause each, both found by asking a Wigley hull for its
metacentric height with freeboard modelled.*

## 1. What was wrong

A hull is generated only up to its design waterline; a GZ curve needs freeboard, so
`tensorcad.gz.with_freeboard` removes the waterline cap, extrudes a vertical band upwards
and closes it with a deck. The band is vertical, so **the waterplane cannot change** --
that is the whole point. Measured:

| | waterplane area at z = 0 |
|---|---|
| Wigley 81x21, capped at the waterline | 6.665625 m^2 |
| same hull + 0.2 m freeboard band | **-13.33125 m^2** (`StabilityError: bounds no net area`) |

`-13.33125` is exactly `-2 x 6.665625`. **Root cause**: `stability._plane_segments` cut the
surface with the plane and took one segment per face. A plane that *contains* a mesh edge
-- the seam ring where the band meets the hull -- is met by **both** of that edge's faces,
so every boundary segment was emitted twice. `_loops` chained the doubled contour into two
coincident loops, and the hole-parity test in `waterplane_moments` (meant for the waterline
through a twin-hull tunnel) saw each loop's first point "inside" the identical other loop,
marked both odd-depth and **subtracted** them: `-A - A = -2A`.

The same loading exposed a second defect. `initial_stability` integrated the *whole closed
mesh* via `hydrostatics.compute_hydrostatics`, which equals the displacement only for a hull
capped at its waterline. With the band on:

| quantity (Wigley 81x21, KG = 0.3 m) | capped hull | banded hull, before | banded hull, after | closed form |
|---|---|---|---|---|
| displacement [m^3] | 2.775608 | **4.108733** (+48 %) | 2.775608 | 4LBT/9 = 2.777778 |
| KB [m] | 0.390674 | 0.499150 | 0.390674 | 5T/8 = 0.390625 |
| GM_transverse [m] | 0.227874 | **0.291834** (+28 %) | 0.227874 | 0.227768 |

Silently: no error, no warning -- precisely the failure mode the module's own docstrings say
they exist to prevent ("a zero waterplane silently turns into an infinite GM two calls
later").

## 2. What was changed (minimal, 3 sites)

1. `tensorcad/stability.py::_plane_segments` -- keep one copy of a segment whose endpoint
   identity pair repeats (a face pair sharing an edge that lies in the plane). One geometric
   curve, one segment; the contour stays a simple ring. Documented in the docstring.
2. `tensorcad/stability.py::initial_stability` -- if the mesh rises above the waterline,
   take the immersed volume and `KB` from `righting_arm.heeled_displacement` at zero heel
   (the same integral the GZ curve uses) instead of the whole-body closed-surface integral.
   A hull capped at its waterline keeps the original fast path **bit for bit**, so every
   previously verified number (box, hemisphere, hull, barge GM to 2e-16) is untouched.
3. `examples/intact_stability/run_example.py` -- the stale comment documenting the
   workaround ("the band duplicates the waterline ring, which the plane cut will not
   accept") now records that both routes agree.

Verification lives in `tests/test_stability.py` (synthetic knuckle ring: two boxes welded at
`z = 0`, where every reference is an exact rational -- area `LB`, `I_T = LB^3/12`,
immersed volume `LBd`, `KB = d/2`, `GM = d/2 + B^2/(12d) - KG`) and in
`tests/test_wigley_gz_verification.py::TestFreeboardSeamWaterplane`.

## 3. Result: GM and the GZ curve from one mesh

Finest grid 161x41, KG = 0.3 m. `GM` is measured on the **banded** mesh -- the same mesh the
curve is computed on -- so the two halves of the check are independent integrals (immersed
volume + waterplane moments vs. heeled volume integral):

| heel [deg] | GZ [m] | GZ/sin(phi) [m] | vs GM = 0.227794 m |
|---|---|---|---|
| 2.0 | 0.007950578 | 0.227813533 | +0.0084% |
| 5.0 | 0.019873787 | 0.228026131 | +0.1017% |
| 10.0 | 0.039731256 | 0.228803184 | +0.4429% |

Convergence of the banded-mesh route to the closed forms:

| grid | capped A [m^2] | banded A [m^2] | rel diff | loops | nabla [m^3] | KB [m] | GM [m] | GM vs closed form |
|---|---|---|---|---|---|---|---|---|
| 41x13 | 6.662500 | 6.662500 | 0.0e+00 | 1 | 2.771222 | 0.390761 | 0.228028 | +0.114% |
| 81x21 | 6.665625 | 6.665625 | 0.0e+00 | 1 | 2.775608 | 0.390674 | 0.227874 | +0.047% |
| 161x41 | 6.666406 | 6.666406 | 0.0e+00 | 1 | 2.777235 | 0.390637 | 0.227794 | +0.012% |

The seam invariance (banded vs capped) is exact to machine precision at every grid; GM
error against the analytic `KB + I_T/nabla - KG` falls 0.114 % -> 0.047 % -> 0.012 %
under refinement. **All key quantities are within 3 % of the reference**; the largest
deviation at the finest grid is 0.012 % (GM), and the 10-degree slope deviates from GM by
0.443 %,
which is the expected wall-sided `(BM/2) tan^2 phi` term, not an error.

Reproduce: `PYTHONPATH=. python examples/freeboard_seam/run_example.py --output
examples/freeboard_seam/outputs` (writes `results.json`, `report.md`, `plot.png`; exits
non-zero if a claim breaks).

## 4. Notes and limits

- The dedup keys on endpoint identity, so it also fixes the general "plane contains a mesh
  edge" case (a chine or a bilge knuckle exactly at the waterline), not just freeboard.
- A mesh whose seam vertices are *not welded* (coincident but distinct vertex indices) is
  still outside the contract: the two contours cannot be chained at all. `with_freeboard`
  reuses the hull's indices, so the intended path is safe.
- The `initial_stability` route switch is decided by "does any vertex rise above the
  waterline", a speed heuristic; the correctness comes from using the immersed integral on
  that branch. A capped hull therefore keeps the pre-existing numbers exactly, and a banded
  hull stops being silently wrong.
- `KB` on the banded branch is exact to the immersed-integral tolerance (1e-9 relative
  below), not to round-off: it comes from the clipped body rather than a closed-form sum.
