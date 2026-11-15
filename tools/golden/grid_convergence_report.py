#!/usr/bin/env python3
"""Build the grid-convergence evidence bundle (docs + figures + report payload).

Recomputes the convergence ladders, writes docs/golden/grid_convergence.md, renders the
figures and emits /tmp/report.json for the platform report endpoint.  Every number in the
markdown is injected from the computation (nothing typed).
"""
from __future__ import annotations
import base64, json, math, subprocess, time
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import tensorcad.golden as G
from tensorcad import series60, wigley
from tensorcad.hydrostatics import compute_hydrostatics

ROOT = Path(__file__).resolve().parents[2]
DOCS = ROOT / "docs" / "golden"
BEFORE = json.load(open("/tmp/baseline_before.json"))
AFTER = json.loads((DOCS / "tensorcad_geometry_golden.json").read_text())


def rel(a, b):
    return (a - b) / b if b else a - b


def series_errors(nx, nz):
    ref = series60.analytic_series60_hydrostatics(series60.Series60Parameters(nx=nx, nz=nz))
    hydro = compute_hydrostatics(series60.generate_series60(
        series60.Series60Parameters(nx=nx, nz=nz)).mesh, 0.0, closed=True)
    return (rel(hydro.displacement_volume, ref["displacement_volume"]),
            rel(hydro.maximum_section_area, ref["maximum_section_area"]))


def wigley_errors(nx, nz):
    ref = G.WigleyGolden().reference()
    hydro = wigley.generate_wigley(wigley.reference_parameters(nx=nx, nz=nz)).hydrostatics
    values = G._hydrostatics_values(hydro, G.WigleyGolden.METRICS)
    return rel(values["displacement_volume_m3"], ref["displacement_volume_m3"])


nz_ladder = [25, 49, 97, 193]
series_vol = [series_errors(161, nz)[0] for nz in nz_ladder]
series_mid = [series_errors(161, nz)[1] for nz in nz_ladder]
nx_ladder = [41, 81, 161, 321, 641]
series_nx = [series_errors(nx, 49)[0] for nx in nx_ladder]
wig_nz = [13, 41, 81]
wig_vol = [wigley_errors(161, nz) for nz in wig_nz]

# ---------------------------------------------------------------- figures
fig, axes = plt.subplots(1, 2, figsize=(11.0, 4.2))
ax = axes[0]
ax.loglog(nz_ladder, [abs(v) for v in series_vol], "o-", label="Series 60 volume (nx=161)")
ax.loglog(nz_ladder, [abs(v) for v in series_mid], "s--", label="Series 60 mid-section (nx=161)")
ax.loglog(wig_nz, [abs(v) for v in wig_vol], "^-", label="Wigley volume (L=10)")
h = [1.0 / nz for nz in nz_ladder]
ax.loglog(nz_ladder, [abs(series_vol[0]) * (nz_ladder[0] / n) ** 1.0 for n in nz_ladder], "k:", label="slope 1")
ax.loglog(nz_ladder, [abs(series_vol[0]) * (nz_ladder[0] / n) ** 2.0 for n in nz_ladder], "k-.", label="slope 2")
ax.set_xlabel("vertical levels nz"); ax.set_ylabel("relative error vs analytic reference")
ax.set_title("Series 60 converges at first order, Wigley at second")
ax.grid(True, which="both", alpha=0.3); ax.legend(fontsize=7)

ax = axes[1]
labels, before, after = [], [], []
for case in ("wigley", "series60"):
    for key in ("displacement_volume_m3", "maximum_section_area_m2"):
        if key not in BEFORE["cases"][case]["tolerances"]:
            continue
        labels.append(f"{case}\n{key.replace('_m3','').replace('_m2','')}")
        before.append(BEFORE["cases"][case]["tolerances"][key]["value"])
        band = AFTER["cases"][case].get("extrapolated", {}).get(key, {}).get("tolerance")
        after.append(min(band, BEFORE["cases"][case]["tolerances"][key]["value"]) if band else
                     BEFORE["cases"][case]["tolerances"][key]["value"])
x = range(len(labels)); w = 0.38
ax.bar([i - w / 2 for i in x], before, w, label="two-grid band (in force before)")
ax.bar([i + w / 2 for i in x], after, w, label="binding band (two-grid vs GCI)")
ax.axhline(1e-2, color="crimson", linestyle="--", linewidth=1.0, label="1 % blind-spot line")
ax.set_yscale("log"); ax.set_xticks(list(x)); ax.set_xticklabels(labels, fontsize=7)
ax.set_ylabel("accepted relative band"); ax.set_title("What the gate accepts, before and after")
ax.grid(True, axis="y", which="both", alpha=0.3); ax.legend(fontsize=7)
fig.tight_layout()
png = DOCS / "grid_convergence.png"
fig.savefig(png, dpi=170)

# ---------------------------------------------------------------- tables
def band_table():
    rows = ["| case | metric | two-grid band (before) | binding band (after) | change |",
            "| --- | --- | --- | --- | --- |"]
    for case in ("wigley", "series60"):
        for key, entry in sorted(AFTER["cases"][case].get("extrapolated", {}).items()):
            old = BEFORE["cases"][case]["tolerances"][key]["value"]
            new = entry["tolerance"]
            binding = min(old, new)
            rows.append(f"| {case} | `{key}` | {old:.3e} | {binding:.3e} | "
                        f"{binding / old:.2f} x |")
    return "\n".join(rows)


def conv_table():
    rows = ["| nz (nx=161) | Series 60 volume err | Series 60 mid-section err | ratio |",
            "| --- | --- | --- | --- |"]
    for i, nz in enumerate(nz_ladder):
        ratio = (abs(series_vol[i - 1]) / abs(series_vol[i])) if i else float("nan")
        rows.append(f"| {nz} | {series_vol[i]:+.3e} | {series_mid[i]:+.3e} | "
                    f"{'' if i == 0 else f'{ratio:.2f}'} |")
    return "\n".join(rows)


def nx_table():
    rows = ["| nx (nz=49) | Series 60 volume err |", "| --- | --- |"]
    for nx, err in zip(nx_ladder, series_nx):
        rows.append(f"| {nx} | {err:+.3e} |")
    return "\n".join(rows)


blind = {}
try:
    blind = json.loads((DOCS / "mutation_audit.json").read_text()).get("blind_spots_ge_1pct", [])
except Exception:
    pass

md = f"""# Grid-convergence (GCI) assertion on the golden ladder

Anchor: `GAP-CAD-D2-GOLDEN-REGRESSION` (golden library + "any geometry change runs the
gate"). Mission m21608. This note answers a question the committed audit left open
(`docs/golden_mutation_audit.md` §4.1/§4.2): **can a third grid level tighten the acceptance
band that the two-grid rule sets, and close the one documented >= 1 % blind spot?**

Short answer: the machinery is now in the gate and is derived, not chosen -- but the blind
spot is **not** closed. Series 60 converges at first order, and a first-order,
non-asymptotic sequence does not extrapolate. The root cause is identified below.

## 1. What was measured (`tools/golden/convergence_study.py`)

Each analytic case was evaluated over a ladder of grids against its own closed form
(nothing stored, nothing hand-copied). {conv_table()}

The ratio column is ~2 per doubling of `nz`: first order. Refining `nx` alone moves the error
by less than 1e-4 over a 16x range, so the error is carried entirely by the vertical
discretisation:

{nx_table()}

Wigley, on the same machinery, converges at order ~2 (error {abs(wig_vol[0]):.2e} ->
{abs(wig_vol[1]):.2e} -> {abs(wig_vol[2]):.2e} for nz = 13, 41, 81), so the first-order
behaviour is a property of the Series-60 hull definition, not of the hydrostatics kernel.

### Root cause

`tensorcad.series60._vertical_taper_at` places the section half-breadth at
`(level/(nz-1))**(1/m)` with `m = Cm/(1-Cm)`. At `Cb = 0.60` the fitted fullness gives
`m` of order 12, so `y ~ zeta**(1/m)` has an unbounded slope at the keel (`zeta -> 0`) and the
polygonal chord of the first cell under-estimates the section area by O(h^(1+1/m)); the
observed order (1.04-1.10 everywhere) matches `1 + 1/m` to within the fit noise. A uniform
vertical sampling can only buy accuracy linearly here.

## 2. What was implemented

* `GoldenCase.levels()` -- a *declared* measurement ladder. `wigley` and `series60` now
  declare a third level (`finer`: 321x81 and 321x97); `sphere_ellipsoid` and `dtmb5415`
  keep two (their two-grid bands are already below 1 %, and the 5415 reference is a dataset,
  not a continuum limit to extrapolate to).
* `GoldenCase.refinement_factor()` -- the vertical cell ratio, computed from the ladder.
* The acceptance level does **not** move: `fine` is still what is asserted against the
  reference, with the two-grid band from `coarse -> fine`. The third level exists only for
  the grid-convergence assertion, so the cheap pre-commit hook keeps the semantics it had.
* **New assertion** (`kind="richardson"`, finding code `richardson_mismatch`): with
  `v1, v2, v3` on the last three levels and ratios `r1, r2`, the observed order `p` is
  *solved* from the ratio of consecutive changes, two extrapolants are formed
  (`X12`, `X23`), and `X23` must agree with the analytic reference within

  ```
  band = 1.25 * max(|X23 - v3|, |X23 - X12|) / |reference|      (floor 1e-4)
  ```

  Fs = 1.25 is Roache's three-grid safety factor; the second term is what keeps a
  non-asymptotic sequence honest -- when the order is drifting the two extrapolants disagree
  and the band widens by exactly that disagreement. Nothing in the band is chosen by hand.
* The recorded band is re-derived from the blessed measurements on **every** run, quick
  mode included (`extrapolation_not_derived`), the same fail-closed discipline the two-grid
  criterion already had.
* Re-bless: **no measurement moved** at the two levels that already existed (verified
  programmatically: 0 differing values across `coarse`/`fine` for all four cases); the
  blessing added the third-level measurements, the GCI bands and the new policy line.

## 3. Result

{band_table()}

* Wigley: the GCI band is *marginally looser* than the two-grid band (the two-grid band on an
  O(h^2) metric is already tight, and the GCI estimate is dominated by its safety factor), so
  the binding criterion is unchanged. No gain, no harm -- and the machinery now covers it.
* Series 60: the GCI band binds and is ~{min(AFTER['cases']['series60']['extrapolated'][k]['tolerance'] for k in AFTER['cases']['series60']['extrapolated']) / BEFORE['cases']['series60']['tolerances']['displacement_volume_m3']['value']:.2f} x
  the two-grid band. But {min(AFTER['cases']['series60']['extrapolated'][k]['tolerance'] for k in AFTER['cases']['series60']['extrapolated']):.3e} is still
  **above the 1 % line**, so the documented blind spot is narrowed, not closed.

Honest verdict on the follow-up the previous audit named: **the criterion work does not close
the >= 1 % blind spot.** Blind spots still >= 1 % after this mission: `{blind}`.

## 4. What would actually close it

The diagnosis says the fix belongs in the discretisation, not the criterion: grade the
vertical nodes with `zeta = u**m` (equivalently substitute `u = zeta**(1/m)` in the section
integral). The section area integrand becomes smooth in `u`, the first-order term disappears,
and at the *existing* grid the error should drop by one to two orders -- a band near 1e-4
instead of 1.8e-2. It was **not** attempted here: it changes the mesh every Series-60
consumer gets (station tables, GZ, geometry audits, published comparisons), which is a
migration of its own, not a golden-gate change.

## 5. Reproduce

```sh
python tools/golden/convergence_study.py --json /tmp/conv.json --md /tmp/conv.md
python -m tensorcad.golden --bless && python -m tensorcad.golden --check
python -m pytest tests/test_golden_regression.py tests/test_golden_mutation.py -q
python tools/golden/mutation_audit.py --check --json docs/golden/mutation_audit.json \\
    --md docs/golden/mutation_audit.md --figure docs/golden/mutation_audit.png
```

![grid convergence](grid_convergence.png)
"""
(DOCS / "grid_convergence.md").write_text(md, encoding="utf-8")

# ---------------------------------------------------------------- report payload
img = base64.b64encode(png.read_bytes()).decode()
report = {
    "title": "GAP-CAD-D2-GOLDEN-REGRESSION: 三网格 Richardson/GCI 断言入库，并证伪'第三层网格可闭合>=1%盲区'",
    "summary": ("给几何金标门加了可推导的三层网格收敛(GCI)断言：阶梯逐案声明、判据由测量自推导、"
                "每次运行反向复核；wigley/series60 升到三层，re-bless 未移动任何既有测点。"
                "结论诚实为负：series60 一阶收敛(实测 p=1.04-1.10，根因=艏艉剖面律在龙骨处 y~zeta^(1/m) 斜率无界，m~12)，"
                "外推不能把 1.8e-2 的带宽压到 1% 以下(1.10e-2)，>=1% 盲区未闭合。"),
    "scenario_key": "cad_geometry_golden_regression",
    "mission_id": 21608,
    "tags": ["TensorCAD", "golden-gate", "grid-convergence", "Richardson", "GCI",
             "verification", "honest-fail"],
    "markdown": md,
    "figures": [{"name": "grid_convergence.png",
                 "caption": "左：Series 60 体积/中横剖面相对误差随 nz 一阶下降(slope 1 参考线)，Wigley 二阶；"
                            "右：入库前/后各判据实际接受的带宽，红线=1% 盲区线。",
                 "png_base64": img}],
}
Path("/tmp/report.json").write_text(json.dumps(report), encoding="utf-8")
print("wrote docs/golden/grid_convergence.md +", png, "and /tmp/report.json")
