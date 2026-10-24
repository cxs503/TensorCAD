"""A design loop: target in, verified geometry and a solver manifest out.

This is the first example in the repository where the *input* is a requirement rather
than a parameter set, and the output is a design someone could act on:

    displacement + a shape ratio + a wetted-area budget
        -> damped Newton over (length, beam, draft)
        -> every iteration measured, validated and hashed
        -> the final hull written as OBJ, gmsh 2.2 and a solver manifest

The three things this script refuses to do, and says so in its own output:

* it does not run a solver. The manifest carries ``solver.executed_here: false`` and the
  reference quantities a flow or structure code needs; the solve belongs to TensorFEM /
  TensorLBM;
* it does not report a converged design as usable without the validation layer agreeing.
  One of the cases below converges to its target exactly and is still refused, because
  its triangles fail the mesh-quality checks;
* it does not hide the discretisation. The displacement is measured on a 41x13 grid, and
  the script re-measures the final design on a doubled grid so the reader can see how much
  of the answer is the grid.

Exit code is non-zero if any acceptance check fails.
"""

from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from tensorcad.adapters.fem import fem_surface, write_fem_json, write_gmsh
from tensorcad.design import DesignProblem, DesignTargets, write_design_json, write_solver_manifest
from tensorcad.export import write_obj
from tensorcad.hull import HullParameters, build_mesh
from tensorcad.hydrostatics import compute_hydrostatics

#: The hull the brief is built from. Its measured values become the targets, so the
#: loop has a known answer to find -- a design exercise with a checkable solution.
BRIEF_HULL = HullParameters(length=12.0, beam=2.4, draft=0.9)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--output", default="", help="write the manifest, OBJ, msh and JSON here")
    args = parser.parse_args(argv)
    output = Path(args.output) if args.output else None
    failures: list[str] = []

    def expect(label: str, condition: bool, detail: str = "") -> None:
        print(f"  [{'PASS' if condition else 'FAIL'}] {label}{': ' + detail if detail else ''}")
        if not condition:
            failures.append(label)

    def measure(parameters: HullParameters):
        mesh, _ = build_mesh(parameters)
        return mesh, compute_hydrostatics(mesh, 0.0, water_density=parameters.water_density)

    print("1. the brief")
    _, brief = measure(BRIEF_HULL)
    print(
        f"     reference hull L = {BRIEF_HULL.length} m, B = {BRIEF_HULL.beam} m, "
        f"T = {BRIEF_HULL.draft} m -> displacement {brief.displacement_volume:.6f} m^3, "
        f"wetted area {brief.wetted_surface_area:.6f} m^2"
    )
    targets = DesignTargets(
        displacement_m3=brief.displacement_volume,
        beam_to_draft=BRIEF_HULL.beam / BRIEF_HULL.draft,
        length_to_beam=BRIEF_HULL.length / BRIEF_HULL.beam,
        wetted_area_m2_max=1.05 * brief.wetted_surface_area,
    )
    print(
        f"     targets: displacement {targets.displacement_m3:.6f} m^3, "
        f"B/T {targets.beam_to_draft:.4f}, L/B {targets.length_to_beam:.4f}, "
        f"wetted area <= {targets.wetted_area_m2_max:.6f} m^2 (a 5 % budget)"
    )
    print("     free variables: L, B, T -- three equalities, so the problem is square")

    print("2. the loop, from a point that is not the answer")
    problem = DesignProblem(
        base=HullParameters(length=10.0, beam=2.0, draft=0.8),
        targets=targets,
        variables=("length", "beam", "draft"),
    )
    started = time.perf_counter()
    solution = problem.solve()
    seconds = time.perf_counter() - started
    print(
        f"     {'it':>3s} {'L (m)':>10s} {'B (m)':>9s} {'T (m)':>8s} {'|r|':>10s} "
        f"{'V (m^3)':>11s} {'A (m^2)':>10s} {'digest':>9s} ok"
    )
    for record in solution.records:
        hydro = record.hydrostatics
        print(
            f"     {record.iteration:3d} {record.variables['length']:10.6f} "
            f"{record.variables['beam']:9.6f} {record.variables['draft']:8.6f} "
            f"{record.residual_norm:10.2e} {hydro.displacement_volume:11.6f} "
            f"{hydro.wetted_surface_area:10.6f} {record.digest[:8]:>9s} "
            f"{'y' if record.passed else 'n'}"
        )
    print(f"     status {solution.status!r} in {seconds:.3f} s, {len(solution.records)} evaluations")
    expect("the loop converges", solution.status == "converged")
    expect("the design passes every gate", solution.solved, f"failed: {solution.failed_checks}")
    final = solution.final
    for name in ("length", "beam", "draft"):
        recovered = getattr(final.parameters, name)
        expected = getattr(BRIEF_HULL, name)
        expect(
            f"{name} is recovered within 1e-7 relative",
            abs(recovered / expected - 1.0) < 1e-7,
            f"{recovered:.9f} vs {expected:.9f}",
        )
    expect(
        "the final displacement matches the target",
        abs(final.hydrostatics.displacement_volume / targets.displacement_m3 - 1.0) < 1e-9,
        f"{final.hydrostatics.displacement_volume:.9f} m^3",
    )
    expect(
        "the wetted-area budget is respected",
        final.hydrostatics.wetted_surface_area <= targets.wetted_area_m2_max,
        f"{final.hydrostatics.wetted_surface_area:.6f} <= {targets.wetted_area_m2_max:.6f} m^2",
    )
    expect(
        "the hull hands off watertight, genus 0",
        bool(final.handoff["watertight"]) and final.handoff["euler_characteristic"] == 2,
        f"{final.handoff['triangles']} triangles, digest {final.digest[:12]}",
    )

    print("3. is it quadratic, and does the Jacobian measure the geometry?")
    order = solution.convergence_order
    print(
        f"     observed order {order:.3f} from the residual norms "
        + ", ".join(f"{record.residual_norm:.1e}" for record in solution.records[-4:])
    )
    for variable, rows in solution.sensitivity().items():
        print(
            f"     d ln(measured)/d ln({variable}): "
            + ", ".join(f"{name} {value:+.6f}" for name, value in rows.items())
        )
    expect(
        "the contraction is quadratic to within 0.2",
        order is not None and 1.8 < order < 2.3,
        f"{order:.3f}" if order is not None else "no estimate",
    )
    sensitivity = solution.sensitivity()
    expect(
        "a shape ratio has sensitivity +1 and -1 to its own terms",
        abs(sensitivity["beam"]["beam_to_draft"] - 1.0) < 1e-9
        and abs(sensitivity["draft"]["beam_to_draft"] + 1.0) < 1e-9,
    )
    expect(
        "displacement scales linearly with beam and draft (sensitivity 1.0)",
        all(abs(sensitivity[v]["displacement_m3"] - 1.0) < 1e-6 for v in ("beam", "draft")),
    )

    print("4. replay and re-open")
    replay = DesignProblem(
        base=HullParameters(length=10.0, beam=2.0, draft=0.8),
        targets=targets,
        variables=("length", "beam", "draft"),
    ).solve()
    print(f"     replay digest {replay.final.digest[:16]} vs {final.digest[:16]}")
    expect("a replay is bit-identical", replay.final.digest == final.digest)
    heavier_targets = DesignTargets(
        displacement_m3=1.15 * brief.displacement_volume,
        beam_to_draft=BRIEF_HULL.beam / BRIEF_HULL.draft,
        length_to_beam=BRIEF_HULL.length / BRIEF_HULL.beam,
        wetted_area_m2_max=1.5 * brief.wetted_surface_area,
    )
    heavier_problem = DesignProblem(
        base=HullParameters(length=12.0, beam=2.0, draft=0.8),
        targets=heavier_targets,
        variables=("length", "beam", "draft"),
    )
    heavier = heavier_problem.solve(start=final.variables)
    print(
        f"     +15 % displacement from the previous answer: {len(heavier.records)} evaluations, "
        f"L = {heavier.final.parameters.length:.6f}, B = {heavier.final.parameters.beam:.6f}, "
        f"T = {heavier.final.parameters.draft:.6f}, V = "
        f"{heavier.final.hydrostatics.displacement_volume:.6f} m^3"
    )
    expect("re-opening a solved design converges", heavier.solved)
    expect(
        "the heavier design carries 15 % more displacement",
        abs(heavier.final.hydrostatics.displacement_volume / heavier_targets.displacement_m3 - 1.0) < 1e-9,
    )
    expect(
        "the heavier design is genuinely bigger",
        heavier.final.parameters.beam > final.parameters.beam,
        f"{heavier.final.parameters.beam:.6f} > {final.parameters.beam:.6f} m",
    )

    print("5. the case the loop must refuse")
    blocked = DesignProblem(
        base=HullParameters(length=10.0, beam=2.0, draft=0.8),
        targets=DesignTargets(displacement_m3=200.0, beam_to_draft=2.5),
    ).solve()
    print(
        f"     target 200 m^3 at L = 10 m: status {blocked.status!r}, |r| = "
        f"{blocked.final.residual_norm:.2e}, solved = {blocked.solved}"
    )
    print(f"     failed gates: {', '.join(blocked.failed_checks) or 'none'}")
    print(
        f"     L/B = {blocked.final.parameters.length / blocked.final.parameters.beam:.3f} "
        f"(the model's guidance is 2.5-8); the mesh is still watertight, which is why only "
        "the quality checks stop it"
    )
    expect("a converged-but-rejected design is not reported as solved", not blocked.solved)
    expect("the rejection names what failed", bool(blocked.failed_checks))
    expect(
        "the manifest marks it not ready for a solver",
        blocked.manifest()["solver"]["status"] == "not_ready",
    )

    print("6. how much of the answer is the grid")
    refined = HullParameters(
        length=final.parameters.length,
        beam=final.parameters.beam,
        draft=final.parameters.draft,
        nx=2 * final.parameters.nx - 1,
        nz=2 * final.parameters.nz - 1,
    )
    _, refined_hydro = measure(refined)
    change = refined_hydro.displacement_volume / final.hydrostatics.displacement_volume - 1.0
    print(
        f"     (41, 13) -> ({refined.nx}, {refined.nz}): displacement "
        f"{final.hydrostatics.displacement_volume:.6f} -> "
        f"{refined_hydro.displacement_volume:.6f} m^3 ({100 * change:+.3f} %), "
        f"wetted area {final.hydrostatics.wetted_surface_area:.6f} -> "
        f"{refined_hydro.wetted_surface_area:.6f} m^2"
    )
    print(
        "     the loop converges on the grid it measures; the record states that grid, and "
        "this is the number to quote when someone asks how exact the design is"
    )
    expect(
        "the design is grid-independent to within 1 %",
        abs(change) < 0.01,
        f"{100 * change:+.3f} %",
    )

    print("7. the handoff")
    manifest = solution.manifest()
    hydro = final.hydrostatics
    print(f"     schema {manifest['schema']}, status {manifest['status']}, solved {manifest['solved']}")
    print(
        f"     reference quantities: waterline z = {hydro.waterline_z:+.3f} m, "
        f"rho = {hydro.water_density:.1f} kg/m^3, S = {hydro.wetted_surface_area:.6f} m^2, "
        f"Cb = {hydro.block_coefficient:.6f}, Cp = {hydro.prismatic_coefficient:.6f}, "
        f"LCB = {hydro.longitudinal_centre_buoyancy:.6f} m, VCB = {hydro.vertical_centre_buoyancy:.6f} m"
    )
    print(
        f"     solver: {manifest['solver']['status']}, executed_here "
        f"{manifest['solver']['executed_here']} -- {manifest['solver']['reason']}"
    )
    print(f"     limits: {'; '.join(manifest['limits'])}")
    expect("the manifest says the solver was not run", manifest["solver"]["executed_here"] is False)
    expect(
        "the manifest carries the geometry handoff",
        manifest["handoff"]["schema"] == "tensorcad.geometry.handoff/1"
        and manifest["handoff"]["mesh_sha256"] == final.digest,
    )

    if output is not None:
        output.mkdir(parents=True, exist_ok=True)
        write_obj(final.mesh, output / "design_hull.obj")
        surface = fem_surface(final.mesh, units=final.mesh.units, coordinate_frame=final.mesh.coordinate_frame)
        write_gmsh(surface, output / "design_hull.msh")
        write_fem_json(surface, output / "design_hull.json")
        write_solver_manifest(solution, output / "solver_manifest.json")
        write_design_json(solution, output / "design_path.json")
        summary = {
            "schema": "tensorcad.design.example/1",
            "brief": targets.__dict__,
            "converged_in": len(solution.records),
            "convergence_order": solution.convergence_order,
            "sensitivity": solution.sensitivity(),
            "final": final.to_dict(),
            "grid_refinement": {
                "resolution": [final.parameters.nx, final.parameters.nz],
                "refined_resolution": [refined.nx, refined.nz],
                "displacement_change": change,
            },
            "rejected_design": {
                "targets": [200.0, 2.5],
                "status": blocked.status,
                "failed": list(blocked.failed_checks),
            },
            "passed": not failures,
        }
        (output / "design_loop.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
        print(f"     wrote {output}/design_hull.obj, .msh, .json, solver_manifest.json, design_path.json")

    print()
    if failures:
        print(f"FAILED ({len(failures)}): " + "; ".join(failures))
        return 1
    print("all design-loop checks passed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
