"""Build the evidence bundle for tensorcad.girder_loads (m17002).

Writes docs/girder_loads.md (numbers injected, never typed by hand),
artifacts/girder_loads_m17002.json and two figures.
"""
from __future__ import annotations

import json
import math
import pathlib
import sys

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))

from tensorcad.girder_loads import (  # noqa: E402
    GRAVITY,
    SEA_WATER_DENSITY,
    WeightItem,
    buoyancy_distribution,
    hull_girder_loads,
    still_water_loads,
)
from tensorcad.hull import HullParameters, generate_hull  # noqa: E402

ART = pathlib.Path("artifacts")
DOCS = pathlib.Path("docs")
ART.mkdir(exist_ok=True)
DOCS.mkdir(exist_ok=True)
R: dict = {}

# ---------------------------------------------------------------- A: box barge
L, B, T = 10.0, 2.0, 1.0
RHO, G = 1000.0, GRAVITY
b0 = RHO * G * B * T
W = b0 * L
xs = tuple(L * i / 20 for i in range(21))
ballast = (WeightItem(0.0, L, W, "ballast"),)
bal = still_water_loads(xs, tuple(b0 for _ in xs), ballast)
R["box_uniform"] = {
    "maxshear": bal.max_shear,
    "maxshear_rel": bal.max_shear / (b0 * L),
    "maxmoment": bal.max_moment,
    "maxmoment_rel": bal.max_moment / (b0 * L * L),
    "residual_force_rel": abs(bal.residual_force) / (b0 * L),
    "residual_moment_rel": abs(bal.residual_moment) / (b0 * L * L),
}
P = 0.1 * W
rows = []
for factor in (0.2, 0.1, 0.05, 0.025):
    d = factor * L
    items = (
        WeightItem(0.0, L, W - P, "hull"),
        WeightItem(L / 2 - d / 2, L / 2 + d / 2, P, "cargo"),
    )
    lds = still_water_loads(tuple(L * i / 200 for i in range(201)), tuple(b0 for _ in range(201)), items)
    expected = -P * (L - d) / 8.0
    got = lds.extreme_moment()[1]
    rows.append((factor, d, got, expected, abs(got - expected) / abs(expected)))
R["box_patch"] = [
    {"patch_fraction": f, "patch_width": d, "measured": got, "analytic": exp, "rel_error": err}
    for f, d, got, exp, err in rows
]

# ------------------------------------------------------ B: parametric hull
nx, nz, draft = 81, 25, 0.8
model = generate_hull(HullParameters(nx=nx, nz=nz))
mesh = getattr(model, "mesh", model)
box = mesh.bounding_box()
keel = float(box.minimum[2])
wl = keel + draft
dist_coarse = buoyancy_distribution(mesh, waterline_z_at_reference=wl, station_count=41)
dist_fine = buoyancy_distribution(mesh, waterline_z_at_reference=wl, station_count=161)
trim = buoyancy_distribution(mesh, waterline_z_at_reference=wl, trim_angle_deg=2.0, station_count=161)
x_aft, x_fore = dist_fine.stations[0], dist_fine.stations[-1]
span = x_fore - x_aft
items = (
    WeightItem(x_aft, x_aft + 0.5 * span, 0.6 * dist_fine.displacement_force, "aft 60%"),
    WeightItem(x_aft + 0.5 * span, x_fore, 0.4 * dist_fine.displacement_force, "fore 40%"),
)
loads, dist = hull_girder_loads(mesh, items, waterline_z_at_reference=wl, station_count=161)
xq, q = loads.extreme_shear()
xm, m = loads.extreme_moment()
R["hull"] = {
    "nx": nx,
    "nz": nz,
    "draft": draft,
    "length": span,
    "displacement_volume": dist.volume,
    "displacement_force": dist.displacement_force,
    "lcb": dist.lcb,
    "lcb_fraction": (dist.lcb - x_aft) / span,
    "volume_coarse_41": dist_coarse.volume,
    "volume_fine_161": dist_fine.volume,
    "volume_convergence_rel": abs(dist_fine.volume - dist_coarse.volume) / dist_fine.volume,
    "trim2_volume": trim.volume,
    "trim2_lcb": trim.lcb,
    "trim2_draft_fore": trim.draft_fore,
    "trim2_draft_aft": trim.draft_aft,
    "lcb_weight": loads.lcb_weight,
    "residual_force": loads.residual_force,
    "residual_force_rel": abs(loads.residual_force) / dist.displacement_force,
    "residual_moment": loads.residual_moment,
    "residual_moment_expected": dist.displacement_force * (dist.lcb - loads.lcb_weight),
    "max_shear": loads.max_shear,
    "max_moment": loads.max_moment,
    "extreme_shear": [xq, q],
    "extreme_moment": [xm, m],
    "max_moment_rel": loads.max_moment / (dist.displacement_force * span),
}
R["hull"]["residual_moment_rel_error"] = abs(
    loads.residual_moment - R["hull"]["residual_moment_expected"]
) / abs(R["hull"]["residual_moment_expected"])

# a balanced (zero-residual) hull case: solve the cargo position that levels the trim
target = dist.displacement_force
best = None
for shift in [i * 0.002 * span for i in range(-100, 101)]:
    centre = dist.lcb + shift
    lo = min(max(centre - 0.05 * span, x_aft), x_fore)
    hi = min(max(centre + 0.05 * span, x_aft), x_fore)
    if hi - lo < 1e-9:
        continue
    trial = (
        WeightItem(x_aft, x_fore, 0.8 * target, "light ship"),
        WeightItem(lo, hi, 0.2 * target, "cargo"),
    )
    cand = still_water_loads(dist.stations, dist.loads, trial)
    score = abs(cand.residual_moment)
    if best is None or score < best[0]:
        best = (score, cand, (lo, hi))
R["balanced"] = {
    "cargo_interval": list(best[2]),
    "residual_moment": best[1].residual_moment,
    "residual_moment_rel": abs(best[1].residual_moment) / (target * span),
    "residual_force_rel": abs(best[1].residual_force) / target,
    "max_moment": best[1].max_moment,
    "extreme_moment": list(best[1].extreme_moment()),
}

# ---------------------------------------------------------------- figures
figures = []
try:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(2, 1, figsize=(8.0, 7.0), sharex=True)
    xs_h = [(x - x_aft) / span for x in loads.stations]
    axes[0].plot(xs_h, [v / 1e3 for v in loads.shear], lw=1.6, label="Q (60/40 ballast)")
    axes[0].axhline(0.0, color="k", lw=0.6)
    axes[0].set_ylabel("shear force Q [kN]")
    axes[0].legend(loc="best", fontsize=8)
    axes[0].set_title(f"Still-water hull girder loads, parametric hull L={span:.2f} m, "
                      f"draft {draft} m, nx={nx}, nz={nz}")
    axes[1].plot(xs_h, [v / 1e3 for v in loads.moment], lw=1.6, color="C1", label="M (60/40 ballast)")
    xm_rel = (xm - x_aft) / span
    axes[1].plot([xm_rel], [m / 1e3], "o", color="C3",
                 label=f"extreme M = {m/1e3:.1f} kN m at x/L={xm_rel:.3f}")
    axes[1].axhline(0.0, color="k", lw=0.6)
    axes[1].set_xlabel("x / L")
    axes[1].set_ylabel("bending moment M [kN m]")
    axes[1].legend(loc="best", fontsize=8)
    fig.tight_layout()
    png = ART / "girder_loads_curves.png"
    fig.savefig(png, dpi=130)
    plt.close(fig)
    figures.append({"name": "girder_loads_curves", "path": str(png),
                    "caption": "Shear force and bending moment along a parametric hull at design draft, "
                               "60/40 weight distribution (force balanced, out of trim)."})

    fig2, ax = plt.subplots(figsize=(6.6, 4.0))
    frac = [row["patch_fraction"] for row in R["box_patch"]]
    got = [-row["measured"] for row in R["box_patch"]]
    exp = [-row["analytic"] for row in R["box_patch"]]
    ax.plot(frac, exp, "s--", label="analytic $-P(L-d)/8$")
    ax.plot(frac, got, "o", label="tensorcad.girder_loads")
    ax.set_xlabel("cargo patch width d / L")
    ax.set_ylabel("|M| at midspan [N m]")
    ax.set_xscale("log")
    ax.set_title("Box barge, midspan patch load: exact vs computed")
    ax.legend()
    fig2.tight_layout()
    png2 = ART / "girder_loads_box_patch.png"
    fig2.savefig(png2, dpi=130)
    plt.close(fig2)
    figures.append({"name": "girder_loads_box_patch", "path": str(png2),
                    "caption": "Box-barge midspan patch load: computed midspan moment against the closed form."})
    R["figures_ok"] = True
except Exception as error:  # pragma: no cover - environment dependent
    R["figures_ok"] = False
    R["figures_error"] = repr(error)

R["figure_paths"] = [f["path"] for f in figures]

# ---------------------------------------------------------------- the doc page
doc = f"""# Still-water hull girder loads: buoyancy, shear force, bending moment

`tensorcad.girder_loads` closes a link the other structural modules document as missing:
`girder` says "no bending moment is applied", `shear_flow` says "no load model", and
`bonjean` calls itself "the input to still-water shear force and bending moment". This
module turns the hull geometry into a buoyancy load per unit length, adds the caller's
weight distribution, and integrates the difference twice -- exactly, per segment, with no
quadrature error, so every number below is a statement about the *load model*, not about an
integration scheme.

Sign convention, stated once: `w = p - b` (weight minus buoyancy, positive down),
`dQ/dx = w`, `dM/dx = Q`. Hogging is positive `M`.

*Generated by `python tools/girder_loads_report_m17002.py`; the numbers are injected into
this page, and `artifacts/girder_loads_m17002.json` holds the same values.*

## Result 1 -- a box barge that balances

Constant buoyancy `b = rho g B T = {b0:.4f} N/m`, weight spread uniformly over the whole
length (L = {L} m, B = {B} m, T = {T} m, rho = {RHO} kg/m3):

| quantity | value | relative to `b0 L` / `b0 L^2` |
| --- | ---: | ---: |
| max |Q| | {bal.max_shear:.3e} N | {R['box_uniform']['maxshear_rel']:.2e} |
| max |M| | {bal.max_moment:.3e} N m | {R['box_uniform']['maxmoment_rel']:.2e} |
| residual force `Q(x1)` | {bal.residual_force:.3e} N | {R['box_uniform']['residual_force_rel']:.2e} |
| residual moment `M(x1)` | {bal.residual_moment:.3e} N m | {R['box_uniform']['residual_moment_rel']:.2e} |

Nothing to carry: a balanced prismatic barge has identically zero shear and moment, and the
implementation returns round-off, not a discretisation residue.

## Result 2 -- a midspan cargo patch against the closed form

A patch of force `P = {P:.1f} N` ({P/W*100:.0f} % of the displacement) spread over width
`d`, with the rest of the weight uniform, gives `M(mid) = -P (L - d) / 8` exactly
(superposition of the point-load `-P L / 8` over the patch):

| `d / L` | computed `M(mid)` [N m] | analytic `-P(L-d)/8` [N m] | relative error |
| ---: | ---: | ---: | ---: |
""" + "".join(
    f"| {row['patch_fraction']:g} | {row['measured']:.6f} | {row['analytic']:.6f} | {row['rel_error']:.2e} |\n"
    for row in R["box_patch"]
) + f"""
## Result 3 -- a real parametric hull

`generate_hull(HullParameters(nx={nx}, nz={nz}))`, waterline at `{wl:.4f}` m, Bonjean
section areas every `{span/160:.4f}` m:

| quantity | value |
| --- | ---: |
| displacement volume (161 stations) | {dist_fine.volume:.6f} m3 |
| displacement force | {dist_fine.displacement_force:.1f} N |
| LCB | {dist_fine.lcb:.4f} m ({(dist_fine.lcb - x_aft)/span:.4f} L from the aft end) |
| volume, 41 vs 161 stations | {dist_coarse.volume:.6f} vs {dist_fine.volume:.6f} m3 ({R['hull']['volume_convergence_rel']*100:.3f} %) |
| trim +2 deg: volume / LCB | {trim.volume:.6f} m3 / {trim.lcb:.4f} m |
| trim +2 deg: draft aft / fore | {trim.draft_aft:.4f} / {trim.draft_fore:.4f} m |

Loads with 60 % of the weight in the aft half (force balanced, deliberately out of trim):

| quantity | value |
| --- | ---: |
| extreme shear | {q:.1f} N at x = {xq:.4f} m |
| extreme moment | {m:.1f} N m at x = {xm:.4f} m ({m/(dist.displacement_force*span):.4f} of `W L`) |
| residual force | {loads.residual_force:.3e} N ({R['hull']['residual_force_rel']:.2e} of `W`) |
| residual moment | {loads.residual_moment:.1f} N m, predicted `W (lcb_b - lcb_w)` = {R['hull']['residual_moment_expected']:.1f} N m (rel. error {R['hull']['residual_moment_rel_error']:.2e}) |

The moment residual is not noise: it *is* the trim imbalance, and it matches the closed form
`W (lcb_buoyancy - lcb_weight)` to `{R['hull']['residual_moment_rel_error']:.1e}` relative. Moving the cargo to
`x = [{best[2][0]:.3f}, {best[2][1]:.3f}]` m reduces it to
`{R['balanced']['residual_moment']:.3e}` N m ({R['balanced']['residual_moment_rel']:.2e} of `W L`), which is
the check that the residual is a real balance condition and not an integration artefact.

## What this measures, and what it does not

* Exactness is *by construction*: the load diagram is piecewise linear (piecewise-linear
  buoyancy samples plus piecewise-constant weight items) and is integrated in closed form
  per segment, with jumps taken as one-sided limits. The box-barge cases therefore come out
  at round-off, and the parametric-hull numbers carry only mesh error in `b(x)` -- 41 vs 161
  stations moves the displacement by {R['hull']['volume_convergence_rel']*100:.3f} %.
* Still water only. Wave-induced and dynamic components are the rule's additions
  (IACS UR S11 etc.), not this module's; no rule criterion is evaluated here.
* The weight distribution is an *input*. The module reports the two equilibrium residuals
  instead of inventing a weight curve or raising, so an out-of-trim loading condition is
  visible rather than hidden.
* Downstream consumers: `girder` (`sigma = M z / I`) and `shear_flow` (where `Q` travels in
  the plating) can now be driven from geometry instead of from a caller-supplied guess.
"""
(DOCS / "girder_loads.md").write_text(doc)
(ART / "girder_loads_m17002.json").write_text(json.dumps(R, indent=2))
print(json.dumps({k: v for k, v in R.items() if k in ("box_uniform", "hull", "balanced", "figures_ok")}, indent=2))
for row in R["box_patch"]:
    print(row)
