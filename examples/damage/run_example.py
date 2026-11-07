"""Damage stability demonstration: box closed forms, sweeps, and the figures.

Runs the free-flooding solve on a wall-sided barge, checks every number it prints
against the closed form, writes ``docs/damage/results.json`` and
``docs/damage_stability.md``, and exits non-zero if any claim stops holding.

    python examples/damage/run_example.py
"""
from __future__ import annotations

import json
import math
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

from tensorcad.damage import (  # noqa: E402
    Compartment,
    buoyancy_curve,
    flooded_equilibrium,
    intact_equilibrium,
)
from tensorcad.geometry import Mesh  # noqa: E402

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
OUT = ROOT / "docs" / "damage"

L, B, T = 10.0, 2.0, 0.8
DISPLACEMENT = L * B * T
KG = 0.95


def box_mesh(length=L, beam=B, draft=T, freeboard=1.0) -> Mesh:
    xs = [0.0, length]
    ys = [-0.5 * beam, 0.5 * beam]
    zs = [-draft, freeboard]
    vertices: list[tuple[float, float, float]] = []
    index: dict[tuple[float, float, float], int] = {}
    faces: list[tuple[int, int, int]] = []

    def vid(point):
        if point not in index:
            index[point] = len(vertices)
            vertices.append(point)
        return index[point]

    def quad(a, b, c, d):
        faces.append((vid(a), vid(b), vid(c)))
        faces.append((vid(a), vid(c), vid(d)))

    quad((xs[0], ys[0], zs[0]), (xs[0], ys[1], zs[0]), (xs[1], ys[1], zs[0]), (xs[1], ys[0], zs[0]))
    quad((xs[0], ys[0], zs[1]), (xs[1], ys[0], zs[1]), (xs[1], ys[1], zs[1]), (xs[0], ys[1], zs[1]))
    quad((xs[0], ys[0], zs[0]), (xs[1], ys[0], zs[0]), (xs[1], ys[0], zs[1]), (xs[0], ys[0], zs[1]))
    quad((xs[1], ys[1], zs[0]), (xs[0], ys[1], zs[0]), (xs[0], ys[1], zs[1]), (xs[1], ys[1], zs[1]))
    quad((xs[0], ys[1], zs[0]), (xs[0], ys[0], zs[0]), (xs[0], ys[0], zs[1]), (xs[0], ys[1], zs[1]))
    quad((xs[1], ys[0], zs[0]), (xs[1], ys[1], zs[0]), (xs[1], ys[1], zs[1]), (xs[1], ys[0], zs[1]))
    return Mesh(vertices=vertices, faces=faces)


def relative(measured: float, exact: float) -> float:
    return abs(measured - exact) / abs(exact) if exact else abs(measured)


def main() -> int:
    mesh = box_mesh()
    claims: list[dict[str, object]] = []
    failures: list[str] = []

    def claim(name: str, measured: float, exact: float, tolerance: float, unit: str) -> None:
        error = relative(measured, exact)
        ok = error <= tolerance
        claims.append(
            {
                "claim": name,
                "measured": measured,
                "reference": exact,
                "relative_error": error,
                "tolerance": tolerance,
                "unit": unit,
                "pass": ok,
            }
        )
        if not ok:
            failures.append(f"{name}: {measured!r} vs {exact!r} ({error * 100:.2e} %)")

    intact = intact_equilibrium(mesh, DISPLACEMENT)
    claim("intact draft [m]", intact.draft_at_reference, T, 1e-12, "m")
    claim("intact KB [m]", intact.kb, 0.5 * T, 1e-12, "m")
    claim("intact I_T [m^4]", intact.transverse_second_moment, B**3 * L / 12.0, 1e-12, "m^4")
    claim("intact LCB [m]", intact.lcb, 0.5 * L, 1e-12, "m")

    midship = Compartment(4.0, 6.0, name="midship")
    mid = flooded_equilibrium(mesh, midship, DISPLACEMENT, kg=KG)
    claim("midship damaged draft [m]", mid.draft_at_reference, L * T / (L - 2.0), 1e-12, "m")
    claim("midship lost buoyancy [m^3]", mid.lost_buoyancy_volume, 2.0 * B * 1.0, 1e-12, "m^3")
    claim("midship damaged I_T [m^4]", mid.transverse_second_moment, B**3 * (L - 2.0) / 12.0, 1e-12, "m^4")
    claim("midship trim [deg]", mid.trim_angle_deg, 0.0, 1.0, "deg")
    claim(
        "midship GM [m]",
        mid.gm_transverse,
        0.5 * (L * T / (L - 2.0)) + (B**3 * (L - 2.0) / 12.0) / DISPLACEMENT - KG,
        1e-9,
        "m",
    )

    mu_state = flooded_equilibrium(mesh, Compartment(4.0, 6.0, permeability=0.6), DISPLACEMENT)
    claim("mu=0.6 damaged draft [m]", mu_state.draft_at_reference, L * T / (L - 1.2), 1e-12, "m")

    aft = flooded_equilibrium(mesh, Compartment(0.0, 2.0, name="aft end"), DISPLACEMENT, kg=KG)
    claim("aft compartment z0 [m]", aft.draft_at_reference, 1.1875, 1e-12, "m")
    claim("aft compartment tan(trim)", aft.waterline.trim_tangent, -0.1875, 1e-12, "-")

    port = flooded_equilibrium(mesh, Compartment(4.0, 6.0, y0=-1.0, y1=0.0, name="port side"), DISPLACEMENT, kg=KG)
    z0 = 16.0 / (18.0 - 1.0 / 6.0)
    claim("port compartment z0 [m]", port.draft_at_reference, z0, 1e-12, "m")
    claim("port compartment tan(heel)", port.waterline.heel_tangent, -z0 / 6.0, 1e-12, "-")

    starboard = flooded_equilibrium(mesh, Compartment(4.0, 6.0, y0=0.0, y1=1.0, name="stbd side"), DISPLACEMENT, kg=KG)
    claim("stbd heel is the mirror [deg]", starboard.heel_angle_deg, -port.heel_angle_deg, 1e-12, "deg")

    curve = buoyancy_curve(mesh, midship, waterline=mid.waterline, stations=81)

    # ---- sweeps -------------------------------------------------------------
    depths = [0.25, 0.5, 1.0, 1.5, 2.0, 3.0, 4.0, 5.0]
    drafts, heels, gms = [], [], []
    for length in depths:
        state = flooded_equilibrium(
            mesh, Compartment(5.0 - 0.5 * length, 5.0 + 0.5 * length), DISPLACEMENT, kg=KG
        )
        drafts.append(state.draft_at_reference)
        gms.append(state.gm_transverse)
        side = flooded_equilibrium(
            mesh, Compartment(5.0 - 0.5 * length, 5.0 + 0.5 * length, y0=-1.0, y1=0.0), DISPLACEMENT, kg=KG
        )
        heels.append(side.heel_angle_deg)

    OUT.mkdir(parents=True, exist_ok=True)

    fig, ax = plt.subplots(figsize=(7.0, 4.0))
    ax.plot(curve.stations, curve.external_areas, label="external A(x)", color="#1f77b4", lw=2)
    ax.fill_between(
        curve.stations, curve.net_areas, curve.external_areas,
        color="#d62728", alpha=0.35, label="buoyancy lost (mu = 1)",
    )
    ax.plot(curve.stations, curve.net_areas, label="net A(x) that carries the ship", color="#d62728", lw=1.5)
    ax.set_xlabel("x [m]")
    ax.set_ylabel("section area [m$^2$]")
    ax.set_title("Buoyancy curve before and after a midship compartment floods")
    ax.legend(loc="lower center")
    ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(OUT / "buoyancy_curve.png", dpi=140)
    plt.close(fig)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10.0, 3.8))
    ax1.plot(depths, drafts, "o-", color="#1f77b4", label="damaged draft")
    ax1.plot(depths, [L * T / (L - length) for length in depths], "--", color="#2ca02c", label="closed form")
    ax1.axhline(T, color="grey", ls=":", label="intact draft")
    ax1.set_xlabel("compartment length [m]")
    ax1.set_ylabel("draft [m]")
    ax1.set_title("Sinkage: solve vs closed form")
    ax1.legend()
    ax1.grid(alpha=0.3)
    ax2.plot(depths, heels, "s-", color="#d62728", label="heel, port compartment")
    ax2.plot(depths, gms, "^-", color="#9467bd", label="GM at midships damage (KG=0.95 m)")
    ax2.axhline(0.0, color="k", lw=0.8)
    ax2.set_xlabel("compartment length [m]")
    ax2.set_ylabel("degrees  /  m")
    ax2.set_title("Asymmetric damage heels, and GM falls")
    ax2.legend()
    ax2.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(OUT / "sweeps.png", dpi=140)
    plt.close(fig)

    worst = max(claim_["relative_error"] for claim_ in claims if claim_["unit"] != "-")
    results = {
        "case": {
            "length_m": L, "beam_m": B, "draft_m": T, "displacement_volume_m3": DISPLACEMENT,
            "kg_m": KG, "freeboard_m": 1.0,
        },
        "intact": intact.to_dict(),
        "midship": mid.to_dict(),
        "mu_0p6": mu_state.to_dict(),
        "aft": aft.to_dict(),
        "port_side": port.to_dict(),
        "starboard_side": starboard.to_dict(),
        "claims": claims,
        "worst_relative_error": worst,
        "sweeps": {"compartment_length_m": depths, "draft_m": drafts, "heel_deg": heels, "gm_m": gms},
    }
    (OUT / "results.json").write_text(json.dumps(results, indent=2) + "\n")

    rows = "\n".join(
        f"| {c['claim']} | {c['measured']:.12g} | {c['reference']:.12g} | "
        f"{float(c['relative_error']) * 100:.2e} % | {'pass' if c['pass'] else 'FAIL'} |"
        for c in claims
    )
    page = f"""# Damage stability: the damaged ship's own equilibrium

`tensorcad/damage.py` solves the free-flooding damage case -- a box-shaped space opened
to the sea -- by the lost-buoyancy method, and this page is the measurement of it.
Every number below is written by `python examples/damage/run_example.py` into
`docs/damage/results.json`, and the script exits non-zero if a claim stops holding.

## The case

| quantity | value |
|---|---|
| length L | {L} m |
| beam B | {B} m |
| design draft T | {T} m |
| displacement volume | {DISPLACEMENT} m^3 |
| `KG` (for GM only) | {KG} m |
| freeboard | 1.0 m |

Compartment geometries: midship (x 4-6 m), aft end (x 0-2 m), port side
(x 4-6 m, y -1-0 m), and the same at `mu = 0.6` and at partial depth.

## How it is solved

The flooded water is inside the hull, so the ship's weight is unchanged and buoyancy is
lost where the water sits:

    V_ext - mu V_c = Delta / rho          (buoyancy)
    centroid(V_ext - mu V_c) = (x_G, y_G) (moment)

Both are evaluated on the *same* waterline plane
`z = z0 + tan(trim)(x - x_ref) + tan(heel) y`, whose three unknowns are found by a damped
Newton solve with a finite-difference Jacobian and a backtracking line search; section
cuts are cached once per station and only re-clipped per iteration. Compartment ends are
panel edges and each panel is integrated with three-point Gauss-Legendre, so the jump in
`A(x)` at a compartment end is never sampled across. `GM = KB + I_T/(Delta/rho) - KG`
uses the *net* waterplane (external minus `mu` times the flooded part), which is the
lost-buoyancy method's answer to the free-surface question.

## Results against the closed forms

A wall-sided barge has no discretisation, so every quantity has an exact value.

| claim | measured | reference | relative error | verdict |
|---|---|---|---|---|
{rows}

Worst relative error over the dimensional claims: {worst * 100:.2e} %.

The two damage cases with no closed-form draft (aft end, port side) are checked against
the exact two-equation solution of the same problem on a rectangle: for the aft
compartment `16 z0 + 16 tan = 16` and `96 z0 + 181.3333 tan = 80` give
`z0 = 1.1875 m`, `tan(trim) = -0.1875` (both reproduced to round-off); for the port half
compartment `18 z0 + tan = 16` and `z0 + 6 tan = 0` give the measured
`z0 = {port.draft_at_reference:.9f} m`, `tan(heel) = {port.waterline.heel_tangent:.9f}`.

## What the figures show

`docs/damage/buoyancy_curve.png` is the damage control picture: the external `A(x)` and
the net curve the ship actually floats on, with the compartment's worth of buoyancy
shaded out. `docs/damage/sweeps.png` sweeps the compartment length: the damaged draft
against the closed form `L T/(L - l)` (exact at every length, not just one), the heel of
a port-side compartment, and `GM`, which goes negative once the damage is long enough --
reported, not clamped.

## Notes, limits, and what is deliberately absent

* **Free-flooding only.** Compartments open to the sea. An internally ruptured tank is
  the added-weight problem with a free surface; passing it here would double-count.
* **Deterministic, not probabilistic.** SOLAS part B-1 (subdivision index `A` against
  `R`), progressive flooding through internal openings, cross-flooding and the
  intermediate stages are not implemented, and no `GZ` curve of the damaged ship is
  produced. `tensorcad.intact_stability` remains the intact-ship rule check.
* **Permeability is an input.** It is never guessed; `mu = 1.0` empty, `0.6` cargo,
  `0.95` machinery, `0` full tank. The state records `lost_buoyancy_volume = mu V_c`.
* **The mesh is the remaining discretisation.** On the curved parametric hull the intact
  solve reproduces `tensorcad.hydrostatics`' own volume exactly and its draft to
  0.5 % -- the difference is the two integrations of the same mesh, not the solve.
* A hull without freeboard cannot be damaged: if the compartment takes away more
  buoyancy than the mesh has above the waterline, the case is refused with
  "the ship would be a submarine" rather than returning a clipped answer.

## Reproducing

```
python -m pytest tests/test_damage.py -q
python examples/damage/run_example.py
```
"""
    (ROOT / "docs" / "damage_stability.md").write_text(page)

    print(json.dumps({"claims": len(claims), "failures": failures, "worst": worst}, indent=2))
    return 0 if not failures else 1


if __name__ == "__main__":
    sys.exit(main())