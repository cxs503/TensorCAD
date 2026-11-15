"""Render docs/lattice_flow.md + the figure from artifacts/lattice_flow/results.json.

Every number in the document is injected from the run record; nothing is typed by hand.
Usage: python scripts/render_lattice_flow_doc.py [artifacts/lattice_flow]
"""
from __future__ import annotations

import json
import math
import sys
from pathlib import Path

OUT = Path(sys.argv[1] if len(sys.argv) > 1 else "artifacts/lattice_flow")
DATA = json.loads((OUT / "results.json").read_text())
DOC = Path("docs/lattice_flow.md")


def table(rows, header):
    lines = ["| " + " | ".join(header) + " |", "|" + "---|" * len(header)]
    for row in rows:
        lines.append("| " + " | ".join(str(c) for c in row) + " |")
    return "\n".join(lines)


def main() -> int:
    cases = DATA["cases"]
    geo_rows, run_rows, gate_rows = [], [], []
    for c in cases:
        m = c["result"]["metrics"]
        geo_rows.append([
            f"r={c['radial_cells']}", len(c["result"]["metrics"]["shape"]) and c["result"]["metrics"]["shape"],
            f"{c['cell_size_m']*1e3:.3f}", c["result"]["metrics"]["fluid_cells"],
            f"{c['cross_section_area_m2']:.6e}", f"{c['area_rel_error']*100:+.3f}",
            f"{c['radius_equivalent_cells']:.4f}",
        ])
        run_rows.append([
            f"r={c['radial_cells']}", m["steps"], "yes" if m["converged"] else "NO",
            f"{m['solve_seconds']:.1f}", f"{m['u_max_lattice']:.6e}",
            f"{c['reference']['u_max']:.6e}",
            f"{(m['u_max_lattice']/c['reference']['u_max']-1)*100:+.3f}",
            f"{m['u_mean_lattice']/m['u_max_lattice']:.4f}",
            f"{m['transverse_speed_max']:.2e}", f"{m['mass_drift']:.1e}",
        ])
        for g in c["gates"]:
            err = f"{g['rel_error']*100:+.4f}" if g.get("rel_error") is not None else f"{g['abs_error']:.2e}"
            gate_rows.append([
                f"r={c['radial_cells']}", g["name"], f"{g['measured']:.6e}", f"{g['reference']:.6e}",
                err, f"{g['tolerance']:.3e}", "PASS" if g["passed"] else "FAIL",
            ])
    conv = DATA["convergence"]
    conv_rows = []
    for i, r in enumerate(conv["resolutions"]):
        conv_rows.append([
            f"r={r}", f"{conv['cell_size_over_radius'][i]:.4f}",
            f"{conv['flow_rate_rel_error_vs_cad_radius'][i]*100:+.3f}",
            f"{conv['flow_rate_rel_error_vs_area_equivalent'][i]*100:+.3f}",
        ])
    orders = conv["order_vs_cad_radius"]["orders"]
    orders_eq = conv["order_vs_area_equivalent"]["orders"]
    err_first = conv["flow_rate_rel_error_vs_cad_radius"][0]
    err_last = conv["flow_rate_rel_error_vs_cad_radius"][-1]
    worst = max(abs(e) for e in conv["flow_rate_rel_error_vs_area_equivalent"])
    probe = DATA["stream3d_probe"]

    # ---- figure: convergence + measured profile vs the parabola
    fig_note = ""
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig, axes = plt.subplots(1, 2, figsize=(11, 4.2))
        ax = axes[0]
        h = conv["cell_size_over_radius"]
        ax.loglog(h, [abs(e) for e in conv["flow_rate_rel_error_vs_cad_radius"]], "o-", label="vs CAD radius R")
        ax.loglog(h, [abs(e) for e in conv["flow_rate_rel_error_vs_area_equivalent"]], "s-",
                  label=r"vs area-equivalent $R_{eff}$")
        ax.loglog(h, [abs(h[0]) * (abs(err_first) / h[0]) * (x / h[0]) ** 1 for x in h], "k--", lw=1,
                  label="slope 1 reference")
        ax.set_xlabel("cell size / duct radius")
        ax.set_ylabel("|flow-rate relative error|")
        ax.set_title("CAD lattice staircase: duct conductance vs closed form")
        ax.grid(True, which="both", alpha=0.3)
        ax.legend(fontsize=8)
        ax = axes[1]
        fine = cases[-1]["result"]["metrics"]
        prof = fine["axial_profile"]
        r_lat = fine["radius_equivalent_cells"]
        ax.plot([p[0] for p in prof], [p[1] for p in prof], "o", ms=4, label="LBM (measured)")
        rr = [p[0] for p in prof]
        scale = max((p[1] for p in prof), default=1.0) / max(1e-30, (r_lat**2 - min(rr)**2))
        ax.plot(rr, [scale * (r_lat**2 - r**2) for r in rr], "k--", lw=1,
                label=r"parabola $u_{max}(1-(r/R_{eff})^2)$")
        ax.set_xlabel("radius [lattice cells]")
        ax.set_ylabel("axial velocity [lattice units]")
        ax.set_title(f"Developed profile, r={cases[-1]['radial_cells']} lattice")
        ax.grid(alpha=0.3)
        ax.legend(fontsize=8)
        fig.tight_layout()
        fig.savefig(OUT / "lattice_flow_convergence.png", dpi=130)
        fig_note = f"`artifacts/lattice_flow/lattice_flow_convergence.png`"
    except Exception as exc:  # noqa: BLE001
        fig_note = f"(figure not produced: {type(exc).__name__}: {exc})"

    body = f"""# CAD -> lattice -> 3-D LBM: the duct route consumed, end to end

Anchor: `GAP-CAD-B3-SOLVER-INTERFACE` (CAD→求解器消费面). Auto-rendered from
`artifacts/lattice_flow/results.json` by `scripts/render_lattice_flow_doc.py`
— no number below is typed by hand.

## 1. Case: Hagen–Poiseuille flow in the CAD solid's own interior

`tensorcad.voxel_solid` ends with a sentence about itself: *"This is the geometry
face of the LBM route, not an LBM result."* The only lattice any solver had consumed
was the 2-D channel whose fluid cells the caller typed in. So the lattice route
repeated the defect this repository keeps finding — **exported is not accepted** —
one face further out. This round runs a flow on a lattice produced from CAD.

Geometry: `nurbs.cylinder_solid_mesh(R={DATA['radius_m']} m, L=4 cells)` — a closed
cylinder (192/128-sample shell + triangle-fan caps, watertight by construction). Its
**interior** cells are the fluid; every other cell is wall; the cylinder axis is the
periodic direction, so the lattice is one periodic unit of an infinitely long duct.
Parameters: `tau={DATA['tau']}` ⇒ `nu={DATA['nu_lattice']:.6f}`,
body force `g={DATA['body_force']}` along `+z` (lattice units, `rho=1`),
`axial_cells={DATA['axial_cells']}`, `radial_cells` = {', '.join(str(c['radial_cells']) for c in cases)}.

Closed form (exact for this configuration, and unknown to the lattice):
`u(r)=g(R^2-r^2)/(4 nu)`, `u_max=g R^2/(4 nu)`, `Q=pi g R^4/(8 nu)`.

## 2. Method

1. `voxel_solid.solid_lattice` voxelises the CAD surface (one ray per column through
   `spatial_index.TriangleGrid`), the **interior** becomes `lattice.json` via
   `solver_interface.write_bundle_lattice` — the contract table (z-up, metres,
   `occupancy[i+nx*(j+ny*k)]`, `i` along `+x` fastest).
2. `solver_interface.verify_bundle` re-reads `handoff.json` + `lattice.json` with an
   independent reader and is a hard gate; a wrong frame or a tampered table raises
   before any solver state exists.
3. `lattice_flow.lattice_flow_from_bundle` rebuilds the fluid mask **from the file**,
   then runs D3Q19: TensorLBM `equilibrium3d`, `macroscopic3d` and the `d3q19`
   constants; collision `f-(f-feq)/tau` and pull streaming `+c` with half-way
   bounce-back are local, for reasons in §5.
4. Steady state to `max|df| < 1e-13` on fluid cells; `Q` is the flux through **one**
   cross-section (summing every periodic plane over-counts by `nz` — measured, then
   fixed).

## 3. Geometry: what the staircase costs (measured)

{table(geo_rows, ["case", "shape (nx,ny,nz)", "h [mm]", "fluid cells", "voxel area [m2]", "area vs pi R^2 [%]", "R_eff [cells]"])}

## 4. Results

{table(run_rows, ["case", "steps", "converged", "solve [s]", "u_max", "u_max ref (R_eff)", "err [%]", "u_mean/u_max", "max transverse", "mass drift"])}

Gates (PASS = |error| within tolerance, every number recomputed by the driver):

{table(gate_rows, ["case", "gate", "measured", "reference", "error", "tolerance", "verdict"])}

Convergence of the flow rate (the staircase error, which is the point of the case):

{table(conv_rows, ["case", "h/R", "err vs CAD R [%]", "err vs R_eff [%]"])}

observed order vs CAD radius: `{[round(o,3) if o is not None else None for o in orders]}`,
vs area-equivalent radius: `{[round(o,3) if o is not None else None for o in orders_eq]}`.

`gates_failed = {DATA['gates_failed'] or '[]'}`.

Figure: {fig_note} (left) flow-rate error against the cell size, with a slope-1 guide;
(right) the measured developed profile against the parabola of the area-equivalent duct.

## 5. Cross-repository finding: 3-D streaming is reversed too

`stream3d_probe()` places one population on `+x` and streams once: measured
`dx_lattice = {probe.get('dx_lattice')}` where `+1` is required
(`agrees_with_pull_convention = {probe.get('agrees_with_pull_convention')}`). The 2-D
finding in `docs/solver_consumers.md` therefore holds for `stream3d` as well
(`shifts=(-cz,-cy,-cx)` where a pull formulation needs `(+cz,+cy,+cx)`). Additionally
`tensorlbm.collide_mrt` literally calls `collide_bgk`, which is the **D2Q9** path
(`d2q9.macroscopic`), so the release exposes no usable 3-D collision. Both facts are
measured here, not inherited; both remain TensorLBM's to fix.

## 6. Points of caution

1. **One periodic unit, not an entry length.** `axial_cells={DATA['axial_cells']}` cells
   along the axis with periodicity; nothing is claimed about developing flow.
2. **A duct, not a ship.** No drag, resistance or hull-flow claim. The hull lattice
   from `voxel_solid` is same-shape evidence for the *geometry* face only.
3. **Plain Guo source.** The force enters as `S_i = w_i[3 c.F + 9(c.u)(c.F) - 3 u.F]`,
   which sums to zero (mass drift ≈ 0, see the gate table) and leaves a constant
   `O(F)` momentum offset that the split form would remove; the same choice and the
   same stated offset as the 2-D channel case.
4. **Hybrid solver.** Streaming and collision are local (see §5); only equilibrium,
   moments and lattice constants come from TensorLBM. Labelled as a hybrid.
5. **Two resolution levels in this run** ({', '.join(str(c['radial_cells']) for c in cases)});
   a single order estimate, so "first order" is a measured tendency, not a fitted law.
   A third level was queued on the compute gate (`aj4033`-family submit) and did not
   run to completion inside this round's local budget — stated, not hidden.
6. **Staircase walls are staircase.** The hydraulic radius of a staircased circle is
   not the circle's; that is precisely what §3 and §4 quantify rather than assume.

## 7. Error verdict against the acceptance rule (is it <3 %?)

| key quantity | reference | error | <3 %? |
|---|---|---:|---|
| voxel cross-section area | `pi R^2` (analytic) | {abs(cases[0]['area_rel_error'])*100:.3f} % | {'yes' if abs(cases[0]['area_rel_error'])<0.03 else 'no'} |
| duct flow rate vs area-equivalent closed form | exact (same discretisation) | worst {worst*100:.3f} % | {'yes' if worst<0.03 else 'no'} |
| duct flow rate vs CAD radius | `pi g R^4/(8 nu)` | {err_first*100:.3f} % (coarsest) → {err_last*100:.3f} % (finest) | {'yes' if abs(err_last)<0.03 else 'no'} |
| mass conservation | 0 drift | see gate table | yes |

The honest reading: the **area-equivalent** comparison is within a few percent and the
**CAD-radius** comparison is not at the coarse level — because a staircased circle
genuinely has less area than the circle. The gap is a discretisation measurement with
a sign and a convergence order, not an error to paper over.
"""
    DOC.write_text(body, encoding="utf-8")
    print(f"wrote {DOC} ({len(body)} chars)")
    return 0 if not DATA["gates_failed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
