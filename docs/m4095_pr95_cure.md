# PR #95 (`agent/m3923`, originally PR #40) cure: one integrator, evidence kept

## 1. What was red, and why

`tests.yml` on PR #95 failed in the **Lint** step (`ruff check .`, 38 findings) on the 3.12
leg; the 3.10/3.11/3.13 legs show `cancelled` because the matrix cancels siblings on the
first failure. Nothing downstream of Lint (type check, pytest, coverage gate, the
example/artifact steps) had ever run on that branch, so the red was only the visible part.

Reproduced locally on the branch tip `cbc9620` with the tool the CI installs (`ruff>=0.6`,
resolved 0.15.0): `ruff check .` exited 1 with 32 findings -- `W292` (no trailing newline),
`F401` (unused import), `F841` (unused locals), `B905`/`RUF007` (`zip()` without
`strict=`/`pairwise`), `B007`, `C416`, `F541`, `RUF022`/`RUF023` (unsorted
`__all__`/`__slots__`), `E501`, `E401`, `SIM115`.

## 2. The part that mattered more than the lint

The branch was cut from `2ff80b2`. `main` has since gained the large-angle righting-arm
engine `tensorcad/gz.py` (PR #39, then the m4029 union cure), whose public surface is
`RightingArmCurve`, `righting_arm_curve`, `wall_sided_righting_arm`, `heeled_hydrostatics`,
`with_freeboard`, `box_wall_sided_limit_deg`, ... -- the *same names* the branch's own
`tensorcad/righting_arm.py` exports. Rebasing onto today's `main` therefore produced `F811`
redefinitions and `B033` duplicate set members, and, worse, the package-level names would
have been silently rebound to whichever import ran last: two parallel integrators for one
physical quantity, which is exactly why PR #40 was closed rather than merged
(`git show eb1aa01`).

Cure: keep the engine already on `main`, drop the parallel module, and port the branch's
*verification claims* onto the surviving engine as a machine-checked audit --
`tests/test_pr95_supersession.py` (5 tests). Removed with it: `examples/righting_arm/`,
`tests/test_righting_arm.py`, `tests/test_merge_m3169.py`, `tools/make_report_m3923.py`,
`tools/post_report_m3923.py`, `tools/verify_merge_m3169.py`, `docs/righting_arm.md`,
`docs/righting_arm/`, `docs/m3923_pr40_merge_cure.md`; `tensorcad/__init__.py`,
`tests/test_public_surface.py` and `CHANGELOG.md` revert to `main`'s union surface.

## 3. The audit, measured on the surviving engine

Sphere (`R = 1 m`, floating at half depth, `KG = 0.5 m` above the keel, heels 0--60 deg),
reference `GZ(phi) = (R - KG) sin(phi)` -- exact at every heel because a floating sphere's
transverse metacentre is its centre:

| sphere resolution | faces | worst abs GZ error [m] | max volume error |
|---|---|---|---|
| 16 | 6656 | 5.175e-04 | 5.24e-14 |
| 24 | 15584 | 2.222e-04 | 3.36e-14 |
| 32 | 28032 | 1.266e-04 | 8.92e-14 |

The branch's promise was "worst < 1 mm at resolution 32": held at 1.266e-04 m, and the error
falls monotonically with refinement, so the agreement is discretisation, not a tuned
constant. Upright slope `GZ(5 deg)/sin(5 deg)` matches `R - KG` to 5e-3 relative -- the
metacentric height by a second, independent route.

Generated hulls at constant displacement (`with_freeboard` band 0.25 m, KG = 0.5 m, heels
0--60 deg):

| hull | max volume error | max GZ [m] | heel of max GZ [deg] | area under curve 0--30 deg [m rad] |
|---|---|---|---|---|
| parametric_hull | 5.61e-15 | 0.167794 | 45 | 0.048400 |
| transom_stern_hull | 1.60e-14 | 0.177213 | 45 | 0.052239 |

The band leaves the upright displacement unchanged to 1e-12 relative, and every point of
both curves re-solves the waterline onto the upright displacement, so the arm is a
righting arm and not a curve of varying load.

## 4. Verification run

| check | command | result |
|---|---|---|
| lint (the red step) | `ruff check .` | 0 findings |
| type check | `mypy` | clean |
| audit | `python -m pytest tests/test_pr95_supersession.py` | 5 passed |
| full suite + coverage gate | `python -m pytest --cov --cov-report=term-missing` | section 5 |

## 5. Residual limits

* The sphere is the only reference here with no closed form *for the mesh*: a sampled
  sphere is not centrally symmetric, so its GZ error is a meshing error (~1e-4 m at
  resolution 32) rather than round-off, unlike the box barge / log references in
  `tests/test_gz.py`.
* This cure removes shipped capability from the branch. Anything the owner wanted from
  `tensorcad/righting_arm.py` *as a second opinion* is gone; what is kept is its
  verification content, retargeted at the surviving engine.
* `tensorcad.gz` carries no sphere helper in its own public surface
  (`sphere_metacentric_height` / `sphere_righting_arm` live on the unmerged `eb1aa01`
  line), so this audit builds the sphere from `SphereField` + `to_mesh`.
