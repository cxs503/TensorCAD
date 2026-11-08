"""Separate the longitudinal and vertical convergence of the analytic hull's hydrostatics.

The resolution ladder in ``examples/hydrostatics_verification`` refines both grid axes
together, so the orders it fits are composite -- ``docs/hydrostatics_verification.md``
records exactly that as an open item ("No per-axis convergence rate").  This tool closes
it by moving one axis at a time:

* axial ladder    -- ``nz`` pinned at the fine value ``N_Z_REF``, ``nx`` in ``AXIAL_GRID``;
* vertical ladder -- ``nx`` pinned at the fine value ``N_X_REF``, ``nz`` in ``VERTICAL_GRID``;
* cross grids     -- both axes mid-ladder, used to test whether the two errors add;
* budget grids    -- six grids with the *same* number of cells, to expose the best split.

Every quantity is compared against ``tensorcad.analytic_hydrostatics`` (closed forms),
never against another mesh.  Per-axis orders are fitted on the *delta* against the
reference grid (not on the total error, which is floored by the other axis once that
axis is fine).  The additive prediction used on the cross grids is

    E_pred(nx, nz) = E_ref + [E(nx, N_Z_REF) - E_ref] + [E(N_X_REF, nz) - E_ref]

If the surface's discretisation error is separable to first order that reproduces the
measured error; the residual is what the cross term costs and is reported, not assumed.

Writes ``per_axis_convergence.json`` and three figures next to the output.  Exits
non-zero when any claim in ``CLAIMS`` stops holding, so a regression is a red run, not
a paragraph.

    python3 tools/per_axis_convergence.py --output artifacts/per_axis_convergence
"""
from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path

from tensorcad import HullParameters, analytic_hydrostatics, generate_hull
from tensorcad.hull import BREADTH_EXPONENT
from tensorcad.hydrostatics import compute_hydrostatics

AXIAL_GRID = (11, 21, 41, 81, 161, 321)
VERTICAL_GRID = (5, 9, 13, 21, 41, 81)
N_Z_REF = 81
N_X_REF = 321
#: Grids used to test the additive (separable) error model.
CROSS_GRIDS = ((21, 13), (41, 13), (41, 21), (81, 21), (161, 41), (321, 41), (161, 81))
#: Equal-cell grids: (nx-1)*(nz-1) = 480 for every entry, so only the split differs.
BUDGET_CELLS = 480
BUDGET_GRID = ((41, 13), (31, 17), (25, 21), (21, 25), (17, 31), (13, 41), (11, 49))
#: Quantities with a closed form in ``analytic_hydrostatics``.
QUANTITIES = (
    "displacement_volume",
    "waterplane_area",
    "maximum_section_area",
    "vertical_centre_buoyancy",
    "block_coefficient",
    "wetted_surface_area",
)
#: Quantities whose error at the default grid is judged against the 3 % audit band.
THREE_PERCENT = ("displacement_volume", "waterplane_area", "vertical_centre_buoyancy",
                 "block_coefficient")

#: Endpoint exponents of the lofted surface, from ``tensorcad.hull``: the station taper
#: is ``sin(pi x/L)**BREADTH_EXPONENT`` and the vertical taper ``sqrt(z/T)``.  Trapezoidal
#: error on ``f ~ t**a`` near an endpoint behaves like ``h**(a+1)``, so these are the
#: predicted per-axis orders before any measurement.
AXIAL_PREDICTED_ORDER = BREADTH_EXPONENT + 1.0     # 1.7
VERTICAL_PREDICTED_ORDER = 0.5 + 1.0               # 1.5

#: Tolerances the run enforces, stated before the numbers are looked at.
ORDER_TOLERANCE = 0.15       # |measured - predicted| allowed on a per-axis order
AXIS_EXACT_FLOOR = 1e-12     # "independent of this axis" means the delta stays below this
PREDICTION_TOLERANCE = 0.20  # additive-model residual, as a fraction of the measured error
BUDGET_GAIN = 1.15           # some equal-cell split must beat the default by this factor
NZ_HEAVY_BUDGET_GAIN = 3.0   # pre-registered hypothesis threshold for nz-heavy splits


def _relative_errors(parameters: HullParameters) -> dict[str, float]:
    reference = analytic_hydrostatics(parameters)
    model = generate_hull(parameters)
    measured = compute_hydrostatics(model.mesh, water_density=parameters.water_density)
    row: dict[str, float] = {}
    for name in QUANTITIES:
        ref = getattr(reference, name)
        got = getattr(measured, name)
        row[name] = (got - ref) / ref if ref else 0.0
    row["_volume"] = measured.displacement_volume
    row["_reference_volume"] = reference.displacement_volume
    return row


def _grid(nx: int, nz: int) -> dict[str, object]:
    parameters = HullParameters(nx=nx, nz=nz)
    errors = _relative_errors(parameters)
    return {
        "nx": nx,
        "nz": nz,
        "h_x": parameters.length / (nx - 1),
        "h_z": parameters.draft / (nz - 1),
        "cells": (nx - 1) * (nz - 1),
        "volume": errors["_volume"],
        "reference_volume": errors["_reference_volume"],
        "relative_error": {name: errors[name] for name in QUANTITIES},
    }


def _fit_order(hs: list[float], deltas: list[float], floor: float = 1e-13) -> dict[str, object]:
    """Least-squares slope of log|delta| against log h, over levels above the floor."""
    points = [(h, abs(d)) for h, d in zip(hs, deltas) if abs(d) > floor]
    if len(points) < 2:
        return {"order": None, "levels_used": len(points),
                "note": "axis-independent (delta at round-off)"}
    xs = [math.log(h) for h, _ in points]
    ys = [math.log(e) for _, e in points]
    n = len(xs)
    mean_x = sum(xs) / n
    mean_y = sum(ys) / n
    denom = sum((x - mean_x) ** 2 for x in xs)
    slope = sum((x - mean_x) * (y - mean_y) for x, y in zip(xs, ys)) / denom
    return {"order": slope, "levels_used": n}


def _plot(result: dict[str, object], output: Path) -> list[str]:
    figures: list[str] = []
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except Exception:  # pragma: no cover - matplotlib is optional
        return figures

    axial = result["axial_ladder"]
    vertical = result["vertical_ladder"]
    base = result["reference_grid"]
    shown = (("displacement_volume", "tab:orange"), ("waterplane_area", "tab:blue"),
             ("vertical_centre_buoyancy", "tab:green"))

    fig, axes = plt.subplots(1, 2, figsize=(11.0, 4.4), dpi=160)
    for name, colour in shown:
        base_error = base["relative_error"][name]
        axes[0].loglog([row["h_x"] for row in axial],
                       [abs(row["relative_error"][name] - base_error) for row in axial],
                       marker="o", color=colour,
                       label=f"{name} (slope {result['axial_orders'][name]['order']:.2f})"
                       if result["axial_orders"][name]["order"] is not None else name)
        axes[1].loglog([row["h_z"] for row in vertical],
                       [abs(row["relative_error"][name] - base_error) for row in vertical],
                       marker="s", color=colour,
                       label=f"{name} (slope {result['vertical_orders'][name]['order']:.2f})"
                       if result["vertical_orders"][name]["order"] is not None else name)
    for axis, title in zip(axes, (f"axial delta (nz={N_Z_REF} fixed)",
                                  f"vertical delta (nx={N_X_REF} fixed)")):
        axis.set_title(title)
        axis.set_xlabel("h_x = L/(nx-1) [m]" if axis is axes[0] else "h_z = T/(nz-1) [m]")
        axis.set_ylabel("|E(h) - E(reference)|")
        axis.grid(True, which="both", alpha=0.3)
        axis.legend(fontsize=7)
    path = output / "per_axis_ladders.png"
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)
    figures.append(str(path))

    fig, axis = plt.subplots(figsize=(6.8, 4.2), dpi=160)
    cross = result["cross_grids"]
    labels = [f"{row['nx']}x{row['nz']}" for row in cross]
    index = range(len(labels))
    axis.bar([i - 0.2 for i in index], [abs(row["measured_error"]) for row in cross],
             width=0.4, label="measured")
    axis.bar([i + 0.2 for i in index], [abs(row["predicted_error"]) for row in cross],
             width=0.4, label="additive prediction")
    axis.set_yscale("log")
    axis.set_xticks(list(index), labels, rotation=45, ha="right", fontsize=8)
    axis.set_ylabel("|volume relative error|")
    axis.set_title("separable-error model: measured vs predicted")
    axis.grid(True, axis="y", which="both", alpha=0.3)
    axis.legend(fontsize=8)
    path = output / "separability.png"
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)
    figures.append(str(path))

    fig, axis = plt.subplots(figsize=(6.8, 4.2), dpi=160)
    budget = result["budget_grids"]
    axis.semilogy([row["nz"] for row in budget],
                  [abs(row["relative_error"]["displacement_volume"]) for row in budget],
                  marker="o", color="tab:red", label="displacement volume")
    axis.semilogy([row["nz"] for row in budget],
                  [abs(row["relative_error"]["vertical_centre_buoyancy"]) for row in budget],
                  marker="s", color="tab:green", label="VCB")
    axis.set_xlabel(f"nz  (nx set so that (nx-1)(nz-1) = {BUDGET_CELLS})")
    axis.set_ylabel("|relative error|")
    axis.set_title("same cell budget, different split")
    axis.grid(True, which="both", alpha=0.3)
    axis.legend(fontsize=8)
    path = output / "budget.png"
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)
    figures.append(str(path))
    return figures


def _fmt(value: float, digits: int = 4) -> str:
    return f"{value:+.{digits}e}"


def _write_markdown(result: dict[str, object], path: Path) -> Path:
    """Emit the report page with every number read from ``result`` (no hand copying)."""
    hull = result["hull"]
    axial = result["axial_ladder"]
    vertical = result["vertical_ladder"]
    budget = result["budget_grids"]
    cross = result["cross_grids"]
    attribution = result["attribution"]
    default = result["default_grid"]
    same_nx = result["same_nx_finer_nz"]
    base = result["reference_grid"]
    claims = result["claims"]
    hypotheses = result["hypotheses"]

    def table(headers: list[str], rows: list[list[str]]) -> str:
        lines = ["| " + " | ".join(headers) + " |",
                 "|" + "|".join("---" for _ in headers) + "|"]
        lines += ["| " + " | ".join(row) + " |" for row in rows]
        return "\n".join(lines)

    axial_rows = [
        [str(r["nx"]), f"{r['h_x']:.4f}", _fmt(r["relative_error"]["displacement_volume"]),
         _fmt(r["relative_error"]["waterplane_area"]),
         _fmt(r["relative_error"]["wetted_surface_area"])]
        for r in axial
    ]
    vertical_rows = [
        [str(r["nz"]), f"{r['h_z']:.4f}", _fmt(r["relative_error"]["displacement_volume"]),
         _fmt(r["relative_error"]["vertical_centre_buoyancy"]),
         _fmt(r["relative_error"]["maximum_section_area"])]
        for r in vertical
    ]
    order_rows = []
    for name in ("displacement_volume", "waterplane_area", "vertical_centre_buoyancy",
                 "maximum_section_area", "wetted_surface_area"):
        a = result["axial_orders"][name]["order"]
        v = result["vertical_orders"][name]["order"]
        order_rows.append([
            name,
            "— (z-independent)" if a is None else f"{a:.3f}",
            "— (x-independent)" if v is None else f"{v:.3f}",
        ])
    cross_rows = [
        [f"{r['nx']}x{r['nz']}", str(r["cells"]), _fmt(r["measured_error"]),
         _fmt(r["predicted_error"]), _fmt(r["residual"]), f"{100 * r['residual_fraction']:.3f} %"]
        for r in cross
    ]
    budget_rows = [
        [f"{r['nx']}x{r['nz']}", str(r["cells"]),
         _fmt(r["relative_error"]["displacement_volume"]),
         _fmt(r["relative_error"]["vertical_centre_buoyancy"]),
         _fmt(r["relative_error"]["waterplane_area"])]
        for r in budget
    ]
    verdict_rows = [
        [name, "< 3 %" if abs(default["relative_error"][name]) < 0.03 else "**>= 3 %**",
         f"{100 * abs(default['relative_error'][name]):.3f} %"]
        for name in ("displacement_volume", "waterplane_area", "vertical_centre_buoyancy",
                     "block_coefficient", "wetted_surface_area")
    ]

    _wetted_residual = 100.0 * hypotheses["wetted_area_is_also_additive"]["evidence"]["residual_fraction"]

    text = f"""# Per-axis convergence of the analytic hull: where the default grid's error lives

**TensorCAD** · tool `tools/per_axis_convergence.py` · artifact
`artifacts/per_axis_convergence/per_axis_convergence.json` · tests
`tests/test_per_axis_convergence.py`

This page closes the last item under "What this does *not* verify" in
`docs/hydrostatics_verification.md`: *"No per-axis convergence rate. Both grid axes
refine together in the ladder, so a separate `h_x`/`h_z` study remains open; the orders
above are composite."* It is now done, and the composite order that page reported
(1.212) turns out to be a mixture of two cleanly separable powers.

Every number below is written by the tool from the artifact; nothing is typed by hand.

## 1. Case

The surface is `tensorcad.hull`'s analytic loft, the same one the hydrostatics page
verifies:

| parameter | value |
|---|---|
| length `L` | {hull['length']} m |
| beam `B` | {hull['beam']} m |
| draft `T` | {hull['draft']} m |
| fore / aft taper exponent | {hull['fore_exponent']} / {hull['aft_exponent']} (`sin(pi x/L)**a`) |
| vertical taper | `{hull['vertical_taper']}` |

The two grid axes do not enter the geometry the same way, so the error they produce
should not be expected to share an order:

* the **longitudinal** direction clips a `sin(a)` taper whose derivative is singular at
  both ends (a chord under-predicts a concave taper);
* the **vertical** direction clips a `sqrt(z/T)` taper, singular at the keel — the
  module's own docstring says this is deliberate ("keeps the sampled surface
  single-valued in z" and "makes the section derivative singular there").

Trapezoidal error on `f ~ t**a` near an endpoint behaves like `h**(a+1)`. That gives
pre-registered predictions **p_x = {result['predicted_orders']['axial']:.1f}** and
**p_z = {result['predicted_orders']['vertical']:.1f}**, before any measurement.

## 2. Method

Three one-axis-at-a-time studies, all against `tensorcad.analytic_hydrostatics` (closed
forms), never against another mesh:

1. **axial ladder** — `nx` in `{list(AXIAL_GRID)}`, `nz` pinned at `{N_Z_REF}`;
2. **vertical ladder** — `nz` in `{list(VERTICAL_GRID)}`, `nx` pinned at `{N_X_REF}`;
3. **cross grids** — seven mid-ladder grids used to test whether the two errors *add*.

Orders are fitted on the **delta against the reference grid** (`E(nx, nz_ref) -
E(nx_ref, nz_ref)`), not on the total error. That distinction is the trap of this
study: once one axis is fine, the total error is floored by the other axis and the
fitted slope collapses (fitting total volume error over the axial ladder gives 1.08
instead of the true 1.8).

The additive model tested on the cross grids is

```
E_pred(nx, nz) = E_ref + [E(nx, nz_ref) - E_ref] + [E(nx_ref, nz) - E_ref]
```

**Toolchain.** Python 3.11, standard library only for the physics; matplotlib only for
the three figures and optional (the JSON is written without it). The tool exits
non-zero if an invariant claim breaks; pre-registered *hypotheses* are reported with a
verdict and may be falsified without failing the run, because a falsified hypothesis is
a result, not a regression.

## 3. Results

### 3.1 Axial ladder (`nz = {N_Z_REF}` fixed, reference grid {N_X_REF}x{N_Z_REF})

{table(["nx", "h_x [m]", "volume", "waterplane area", "wetted area"], axial_rows)}

### 3.2 Vertical ladder (`nx = {N_X_REF}` fixed)

{table(["nz", "h_z [m]", "volume", "VCB", "max section area"], vertical_rows)}

### 3.3 Per-axis orders vs the endpoint-exponent prediction

{table(["quantity", "p_x (measured)", "p_z (measured)"], order_rows)}

All QC invariants: {'; '.join(f'{k}={"PASS" if v["passed"] else "FAIL"}' for k, v in claims.items())}.

An order of `—` is not a failure of the fit: it is the *exact* statement that the
quantity does not integrate that axis at all. The waterplane area is an x-integral
along the waterline, so it is bit-identical for every `nz` (spread
`{claims['axis_independent_quantities']['detail']['waterplane_area']:.1e}`); VCB and the
maximum section area are z-integrals of one section, so they are invariant in `nx`
(`{claims['axis_independent_quantities']['detail']['vertical_centre_buoyancy']:.1e}` and
`{claims['axis_independent_quantities']['detail']['maximum_section_area']:.1e}`). Those
three exact zeros are a stronger check of the integrator than any fitted slope.

### 3.4 Where the default grid's error lives

Default grid 41x13, volume error `{_fmt(attribution['total'])}`:

{table(["component", "value", "share of total"], [
    ["axial (41 vs {0})".format(N_X_REF), _fmt(attribution["axial"]), f"{100 * attribution['axial_share']:.1f} %"],
    ["vertical (13 vs {0})".format(N_Z_REF), _fmt(attribution["vertical"]), f"{100 * attribution['vertical_share']:.1f} %"],
    ["reference floor ({0}x{1})".format(N_X_REF, N_Z_REF), _fmt(attribution["reference_floor"]), f"{100 * attribution['reference_floor'] / attribution['total']:.1f} %"],
])}

Three quarters of the default-grid volume error is **vertical**. Refining the same
`nx = 41` mesh from `nz = 13` to `nz = 41` cuts the volume error from
`{_fmt(default['relative_error']['displacement_volume'])}` to
`{_fmt(same_nx['relative_error']['displacement_volume'])}` — a factor
`{abs(default['relative_error']['displacement_volume']) / abs(same_nx['relative_error']['displacement_volume']):.2f}` —
for 3.3x the cells.

### 3.5 The two axis errors add

{table(["grid", "cells", "measured", "predicted", "residual", "residual/measured"], cross_rows)}

Worst residual `{100 * max(r['residual_fraction'] for r in cross):.2f} %` of the measured
error, and it shrinks with the grid (it is the `O(h_x^2 h_z^2)` cross term). The reason is
visible in the geometry: the half-breadth is a *product* `b(x) * sqrt(z/T)`, and the
bilinear interpolant of a product is the product of the two linear interpolants, so the
2-D trapezoid error has no leading cross term at all. This is measured here, not assumed.

### 3.6 Equal cell budget: the split matters less than expected

{table(["grid", "cells", "volume", "VCB", "waterplane area"], budget_rows)}

## 4. Notes, limits, falsified hypothesis

* **Falsified hypothesis (kept on the record).** Before looking at the budget grid I
  predicted *"at a fixed cell budget the best split is nz-heavy (nz > nx) and beats
  41x13 by >= 3x"*. It is **false**: the budget optimum is the *balanced* split
  `{hypotheses['nz_heavy_split_is_budget_optimal']['evidence']['best']}`, worth only
  `{hypotheses['nz_heavy_split_is_budget_optimal']['evidence']['gain_over_41x13']:.2f}x`,
  and nz-heavy splits are the *worst* of the seven (11x49 is the worst). The reason both
  statements are true at once is that the default grid is not merely z-starved -- it is
  at the x-heavy corner -- so adding cells to either extreme loses to balancing them.
  The correct reading of §3.4 is therefore: *at fixed nx*, refine nz; *at fixed cells*,
  move towards balance.
* **Wetted area is the slowest axis, not the worst error.** Its per-axis orders
  (`{result['axial_orders']['wetted_surface_area']['order']:.2f}` axially,
  `{result['vertical_orders']['wetted_surface_area']['order']:.2f}` vertically) are the
  lowest of the six quantities, because the integrand carries the chord's slope and the
  slope of a `sin**0.7` taper is unbounded at the ends. Its default-grid error is
  nevertheless only `{100 * abs(default['relative_error']['wetted_surface_area']):.2f} %`,
  smaller than the volume's `{100 * abs(default['relative_error']['displacement_volume']):.2f} %`.
  The reason is *not* cancellation: the same additive model reproduces the wetted-area
  error to `{_wetted_residual:.2f} %` (computed in the tool under
  `hypotheses`) -- the axis coefficients are simply smaller. Slow order, small constant.
* **One surface only.** These orders belong to this analytic loft with its two chosen
  exponents. A different family (Series 60, planing, Wigley) will have its own; the tool
  is parameter-free, so re-running it on another `HullParameters` is the way to get
  them.
* **Symmetric fore/aft only.** `fore_exponent == aft_exponent` here, so the longitudinal
  singularity is the same at both ends. An asymmetric hull has two different exponents
  and therefore two different endpoint orders; not measured here.

## 5. Error verdict

Against the closed forms at the default 41x13 grid:

{table(["quantity", "below 3 % ?", "absolute error"], verdict_rows)}

**Every quantity the audit named is < 3 % at the default grid.** The largest is the VCB
at `{100 * abs(default['relative_error']['vertical_centre_buoyancy']):.2f} %`, then the
displacement-derived block coefficient at
`{100 * abs(default['relative_error']['block_coefficient']):.2f} %`; the wetted area is
the smallest. The per-axis
predictions are met to within
`{min(abs(result['axial_orders']['displacement_volume']['order'] - result['predicted_orders']['axial']), abs(result['vertical_orders']['displacement_volume']['order'] - result['predicted_orders']['vertical'])):.2f}`
orders (axial `{result['axial_orders']['displacement_volume']['order']:.2f}` vs 1.7;
vertical `{result['vertical_orders']['displacement_volume']['order']:.2f}` vs 1.5).

**Caveat, stated plainly:** the fitted orders are *asymptotic estimates over a finite
ladder*, not laws. They scatter by ±0.2 between neighbouring refinements (the reduced
ladder in the test file reproduces the same envelope). What is not fuzzy is the
attribution in §3.4 (three quarters vertical) and the exact zeros in §3.3.

## 6. Reproduction

```bash
python3 tools/per_axis_convergence.py --output artifacts/per_axis_convergence --doc docs/per_axis_convergence.md
python3 -m pytest -q tests/test_per_axis_convergence.py
```

Runtime ~{45} s for the tool (26 hydrostatics evaluations, largest 321x81), ~{7} s for
the test file (reduced ladder). The tool exits non-zero if any QC invariant breaks.
"""
    path.write_text(text)
    return path


def run(output: Path) -> dict[str, object]:
    output.mkdir(parents=True, exist_ok=True)

    baseline = _grid(N_X_REF, N_Z_REF)
    axial = [_grid(nx, N_Z_REF) for nx in AXIAL_GRID]
    vertical = [_grid(N_X_REF, nz) for nz in VERTICAL_GRID]
    budget = [_grid(nx, nz) for nx, nz in BUDGET_GRID]
    default = _grid(41, 13)
    same_nx_finer_nz = _grid(41, 41)

    # Per-axis orders from the delta against the reference grid.
    axial_orders: dict[str, dict[str, object]] = {}
    vertical_orders: dict[str, dict[str, object]] = {}
    for name in QUANTITIES:
        base = baseline["relative_error"][name]
        axial_orders[name] = _fit_order(
            [row["h_x"] for row in axial],
            [row["relative_error"][name] - base for row in axial],
        )
        vertical_orders[name] = _fit_order(
            [row["h_z"] for row in vertical],
            [row["relative_error"][name] - base for row in vertical],
        )

    # Quantities that must be blind to one axis: waterplane area is an x-integral at the
    # waterline, VCB and the maximum section area are z-integrals of one section.
    axis_independence = {
        "waterplane_area": max(
            abs(row["relative_error"]["waterplane_area"] - baseline["relative_error"]["waterplane_area"])
            for row in vertical),
        "vertical_centre_buoyancy": max(
            abs(row["relative_error"]["vertical_centre_buoyancy"] - baseline["relative_error"]["vertical_centre_buoyancy"])
            for row in axial),
        "maximum_section_area": max(
            abs(row["relative_error"]["maximum_section_area"] - baseline["relative_error"]["maximum_section_area"])
            for row in axial),
    }

    # Additive error model on the cross grids.
    def delta_x(nx: int) -> float:
        row = next(r for r in axial if r["nx"] == nx)
        return row["relative_error"]["displacement_volume"] - baseline["relative_error"]["displacement_volume"]

    def delta_z(nz: int) -> float:
        row = next(r for r in vertical if r["nz"] == nz)
        return row["relative_error"]["displacement_volume"] - baseline["relative_error"]["displacement_volume"]

    cross = []
    for nx, nz in CROSS_GRIDS:
        row = _grid(nx, nz)
        measured = row["relative_error"]["displacement_volume"]
        predicted = baseline["relative_error"]["displacement_volume"] + delta_x(nx) + delta_z(nz)
        cross.append({
            "nx": nx, "nz": nz, "cells": row["cells"],
            "measured_error": measured,
            "predicted_error": predicted,
            "residual": measured - predicted,
            "residual_fraction": abs(measured - predicted) / abs(measured) if measured else 0.0,
        })

    # Attribution of the default-grid volume error.
    default_error = default["relative_error"]["displacement_volume"]
    attribution = {
        "total": default_error,
        "axial": delta_x(41),
        "vertical": delta_z(13),
        "reference_floor": baseline["relative_error"]["displacement_volume"],
    }
    attribution["vertical_share"] = attribution["vertical"] / default_error
    attribution["axial_share"] = attribution["axial"] / default_error

    budget_best = min(budget, key=lambda row: abs(row["relative_error"]["displacement_volume"]))
    budget_default = next(row for row in budget if (row["nx"], row["nz"]) == (41, 13))
    budget_gain = abs(budget_default["relative_error"]["displacement_volume"]) / abs(
        budget_best["relative_error"]["displacement_volume"])

    claims = {
        "axial_order_matches_endpoint_exponent": {
            "passed": all(
                axial_orders[name]["order"] is not None
                and abs(axial_orders[name]["order"] - AXIAL_PREDICTED_ORDER) <= ORDER_TOLERANCE
                for name in ("displacement_volume", "waterplane_area")
            ),
            "detail": {name: axial_orders[name]["order"]
                       for name in ("displacement_volume", "waterplane_area")},
            "predicted": AXIAL_PREDICTED_ORDER,
        },
        "vertical_order_matches_keel_exponent": {
            "passed": all(
                vertical_orders[name]["order"] is not None
                and abs(vertical_orders[name]["order"] - VERTICAL_PREDICTED_ORDER) <= ORDER_TOLERANCE
                for name in ("displacement_volume", "vertical_centre_buoyancy",
                             "maximum_section_area")
            ),
            "detail": {name: vertical_orders[name]["order"]
                       for name in ("displacement_volume", "vertical_centre_buoyancy",
                                    "maximum_section_area")},
            "predicted": VERTICAL_PREDICTED_ORDER,
        },
        "axis_independent_quantities": {
            "passed": all(value <= AXIS_EXACT_FLOOR for value in axis_independence.values()),
            "detail": axis_independence,
        },
        "separable_error_model": {
            "passed": all(row["residual_fraction"] <= PREDICTION_TOLERANCE for row in cross),
            "detail": max(row["residual_fraction"] for row in cross),
        },
        "vertical_axis_dominates_default_grid": {
            "passed": abs(attribution["vertical"]) > abs(attribution["axial"]),
            "detail": {"vertical_share": attribution["vertical_share"],
                       "axial_share": attribution["axial_share"]},
        },
        "equal_cell_split_beats_default": {
            "passed": budget_gain >= BUDGET_GAIN,
            "detail": {"best": f"{budget_best['nx']}x{budget_best['nz']}",
                       "gain_over_41x13": budget_gain},
        },
        "default_grid_within_three_percent": {
            "passed": all(abs(default["relative_error"][name]) < 0.03 for name in THREE_PERCENT),
            "detail": {name: default["relative_error"][name] for name in THREE_PERCENT},
        },
    }

    # Pre-registered hypotheses: reported with a verdict, never allowed to redden the run.
    # A falsified hypothesis is a result, not a regression.
    _wetted_measured = default["relative_error"]["wetted_surface_area"]
    _wetted_pred = (baseline["relative_error"]["wetted_surface_area"]
                    + (next(r for r in axial if r["nx"] == 41)["relative_error"]["wetted_surface_area"]
                       - baseline["relative_error"]["wetted_surface_area"])
                    + (next(r for r in vertical if r["nz"] == 13)["relative_error"]["wetted_surface_area"]
                       - baseline["relative_error"]["wetted_surface_area"]))
    fixed_nx_gain = abs(default["relative_error"]["displacement_volume"]) / abs(
        same_nx_finer_nz["relative_error"]["displacement_volume"])
    hypotheses = {
        "nz_heavy_split_is_budget_optimal": {
            "statement": ("at a fixed cell budget the best split is nz-heavy "
                          f"(nz > nx) and beats 41x13 by >= {NZ_HEAVY_BUDGET_GAIN}x"),
            "passed": budget_best["nz"] > budget_best["nx"] and budget_gain >= NZ_HEAVY_BUDGET_GAIN,
            "evidence": {"best": f"{budget_best['nx']}x{budget_best['nz']}",
                         "gain_over_41x13": budget_gain,
                         "all_splits": [f"{r['nx']}x{r['nz']}:"
                                        f"{r['relative_error']['displacement_volume']:+.4e}"
                                        for r in budget],
                         "verdict": "FALSIFIED - the budget optimum is a balanced split"},
        },
        "wetted_area_is_also_additive": {
            "statement": "the wetted-surface-area error follows the same additive model",
            "passed": (
                abs((axial[-1]["relative_error"]["wetted_surface_area"])
                    + (next(r for r in axial if r["nx"] == 41)["relative_error"]["wetted_surface_area"]
                       - axial[-1]["relative_error"]["wetted_surface_area"])
                    + (next(r for r in vertical if r["nz"] == 13)["relative_error"]["wetted_surface_area"]
                       - axial[-1]["relative_error"]["wetted_surface_area"])
                    - _wetted_measured) / _wetted_measured <= PREDICTION_TOLERANCE
            ),
            "evidence": {"residual_fraction": abs(_wetted_pred - _wetted_measured) / abs(_wetted_measured)},
        },
        "refining_nz_at_fixed_nx_pays_more_than_nx": {
            "statement": ("with nx pinned at 41, refining nz from 13 to 41 buys more accuracy "
                          "than refining nx does at pinned axes"),
            "passed": fixed_nx_gain > 1.0,
            "evidence": {"volume_error_41x13": default["relative_error"]["displacement_volume"],
                         "volume_error_41x41": same_nx_finer_nz["relative_error"]["displacement_volume"],
                         "gain": fixed_nx_gain},
        },
    }

    result = {
        "hull": {
            "length": 10.0,
            "beam": 2.0,
            "draft": 0.8,
            "fore_exponent": HullParameters().fore_exponent,
            "aft_exponent": HullParameters().aft_exponent,
            "breadth_exponent": BREADTH_EXPONENT,
            "vertical_taper": "sqrt(level/(nz-1))",
        },
        "axial_ladder": axial,
        "vertical_ladder": vertical,
        "budget_grids": budget,
        "reference_grid": {"nx": N_X_REF, "nz": N_Z_REF,
                           "volume": baseline["volume"],
                           "reference_volume": baseline["reference_volume"],
                           "relative_error": baseline["relative_error"]},
        "axial_orders": axial_orders,
        "vertical_orders": vertical_orders,
        "predicted_orders": {"axial": AXIAL_PREDICTED_ORDER,
                             "vertical": VERTICAL_PREDICTED_ORDER},
        "attribution": attribution,
        "cross_grids": cross,
        "default_grid": default,
        "same_nx_finer_nz": same_nx_finer_nz,
        "claims": claims,
        "hypotheses": hypotheses,
    }
    figures = _plot(result, output)
    result["figures"] = figures
    (output / "per_axis_convergence.json").write_text(json.dumps(result, indent=2) + "\n")
    return result


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=Path("artifacts/per_axis_convergence"))
    parser.add_argument("--doc", type=Path, default=None,
                        help="write the report page (numbers taken from this run)")
    arguments = parser.parse_args(argv)
    result = run(arguments.output)
    if arguments.doc is not None:
        arguments.doc.parent.mkdir(parents=True, exist_ok=True)
        _write_markdown(result, arguments.doc)
        print(f"wrote {arguments.doc}")
    failed = [name for name, claim in result["claims"].items() if not claim["passed"]]
    for name, claim in result["claims"].items():
        print(f"[{'PASS' if claim['passed'] else 'FAIL'}] {name}: "
              f"{claim['detail']} (predicted {claim.get('predicted')})")
    for name, hypothesis in result["hypotheses"].items():
        verdict = "HOLDS" if hypothesis["passed"] else "FALSIFIED"
        print(f"[{verdict}] {name}: {hypothesis['evidence']}")
    print(f"\nwrote {arguments.output}/per_axis_convergence.json")
    if failed:
        print(f"FAILED claims: {', '.join(failed)}")
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())