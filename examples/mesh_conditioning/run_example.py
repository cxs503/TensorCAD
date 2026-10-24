"""Mesh conditioning on the AFF-8 model -- the numbers, including the bad ones.

The claim this example tests is narrow and checkable: marching-tet output is watertight
but full of slivers, and smoothing plus a projection back onto the analytic field can
remove the slivers *without moving the shape*. "Without moving the shape" is measured
against the analytic wetted volume and area of the AFF-8 model, not against the input
mesh -- the input mesh is the thing being repaired.

What the script does, in order:

1. Samples SUBOFF with sail and fins into a closed mesh and reports its element quality
   (the input has slivers; that is the point).
2. Runs Taubin smoothing. The lam/mu pair is a band-pass, so unlike a plain Laplacian it
   does not shrink the body: the volume and area drift are printed, not assumed.
3. Projects every vertex back onto the fields with a bounded, backtracking Newton step.
   The step limit is per vertex (a fraction of that vertex's own mean edge length), which
   is what stopped the projection from walking crease vertices off the surface.
4. Decimates the conditioned surface to several face budgets and reports quality, drift
   and time for each.
5. Runs the *other* order -- decimate first, then condition, the ``condition_mesh`` path
   used by the test suite -- and prints its numbers beside the first, because at the same
   face count the decimate-first result is measurably worse. Both are reported; only one
   is recommended.
6. Hands the final surface to the FEM exporter (gmsh 2.2 + JSON) and prints the diff
   between the input and output handoff records.

Exit code is non-zero if any acceptance check fails: watertight, edge-manifold, Euler
characteristic unchanged, no zero-area faces, volume drift within 0.5 %, area drift
within 1 %, and no sliver left. Element quality is reported and *not* gated where the
measurement says it must degrade: decimating a 0.008 degree needle by a factor of two
cannot produce a 3 degree triangle out of the needle's own neighbourhood, and pretending
otherwise is how a "quality gate" becomes a rubber stamp. The one number that is gated
is the sliver share, because a sliver is a topological defect (a degenerate fold), not a
taste judgement about angles.
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
from tensorcad.export import write_obj
from tensorcad.geometry_field import to_mesh
from tensorcad.mesh_conditioning import (
    condition_mesh,
    decimate,
    measure_mesh,
    project_to_field,
    taubin_smooth,
)
from tensorcad.suboff_appendages import SuboffFinsField, SuboffFullField, SuboffHullField, SuboffSailField

#: AFF-8 reference volume and wetted area from the analytic fields (docs/suboff_appendages.md).
REFERENCE_VOLUME_M3 = 0.70754064
REFERENCE_AREA_M2 = 6.350646


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--resolution", default="150,18,26", help="cells per axis for the input sample")
    parser.add_argument("--targets", default="45000,30000,20000", help="face budgets for decimation")
    parser.add_argument("--output", default="", help="write OBJ/msh/JSON here")
    args = parser.parse_args(argv)

    resolution = tuple(int(part) for part in args.resolution.split(","))
    targets = [int(part) for part in args.targets.split(",")]
    output = Path(args.output) if args.output else None
    failures: list[str] = []

    def expect(label: str, condition: bool, detail: str = "") -> None:
        print(f"  [{'PASS' if condition else 'FAIL'}] {label}{': ' + detail if detail else ''}")
        if not condition:
            failures.append(label)

    def drift(measured: float, reference: float) -> float:
        return 100.0 * (measured / reference - 1.0)

    def row(label: str, mesh, seconds: float | None = None) -> None:
        q = measure_mesh(mesh)
        timing = f"{seconds:6.1f}" if seconds is not None else "     -"
        print(
            f"     {label:>16s} {q.triangles:8d} {q.min_angle_degrees:9.4f} "
            f"{q.max_aspect_ratio:9.1f} {100.0 * q.sliver_share:8.4f} "
            f"{drift(q.signed_volume, REFERENCE_VOLUME_M3):+8.3f} "
            f"{drift(q.surface_area, REFERENCE_AREA_M2):+8.3f} {timing}"
        )

    print("1. input: marching tets on the analytic AFF-8 field")
    started = time.perf_counter()
    sampled = to_mesh(SuboffFullField(), resolution, margin=2)
    seconds_sample = time.perf_counter() - started
    incoming = sampled.mesh
    print(
        f"     {len(incoming.vertices)} vertices, {len(incoming.faces)} triangles in "
        f"{seconds_sample:.2f} s; sampler checks passed: {sampled.passed}"
    )
    print(
        f"     {'stage':>16s} {'tris':>8s} {'min ang':>9s} {'max ar':>9s} {'slivers%':>8s} "
        f"{'dvol %':>8s} {'darea %':>8s} {'s':>6s}"
    )
    row("input", incoming)
    before = measure_mesh(incoming)
    expect("the input is watertight", before.watertight)
    expect(
        "the input is the sliver-heavy mesh this example exists for",
        before.sliver_share > 0.005,
        f"{100 * before.sliver_share:.3f} % of triangles are slivers",
    )

    print("2. Taubin smoothing (band-pass: lambda 0.5, mu -0.53)")
    started = time.perf_counter()
    smoothed, smoothing_move = taubin_smooth(incoming, iterations=30)
    seconds_smooth = time.perf_counter() - started
    row("smooth", smoothed, seconds_smooth)
    after_smooth = measure_mesh(smoothed)
    print(
        f"     max vertex displacement {1e3 * smoothing_move:.2f} mm; min angle improved "
        f"{before.min_angle_degrees:.4f} deg -> {after_smooth.min_angle_degrees:.4f} deg "
        f"(x{after_smooth.min_angle_degrees / before.min_angle_degrees:.0f})"
    )
    expect("smoothing keeps the mesh watertight", after_smooth.watertight)
    expect(
        "smoothing removes the slivers",
        after_smooth.sliver_share == 0.0,
        f"{100 * after_smooth.sliver_share:.4f} % left",
    )
    expect(
        "smoothing alone keeps the volume within 0.5 %",
        abs(drift(after_smooth.signed_volume, REFERENCE_VOLUME_M3)) < 0.5,
        f"{drift(after_smooth.signed_volume, REFERENCE_VOLUME_M3):+.3f} %",
    )

    print("3. projection back onto the analytic fields")
    started = time.perf_counter()
    projected, projection_move = project_to_field(smoothed, SuboffFullField(), iterations=30)
    seconds_project = time.perf_counter() - started
    row("project", projected, seconds_project)
    after_project = measure_mesh(projected)
    print(f"     max vertex displacement {1e3 * projection_move:.3f} mm")
    residual = SuboffFullField().signed_distance(projected.vertices)
    worst = max(abs(value) for value in residual)
    print(f"     worst |f| over {len(residual)} vertices: {worst * 1e3:.4f} mm")
    expect("projection keeps the mesh watertight", after_project.watertight)
    expect(
        "the conditioned surface sits on the analytic field",
        worst < 1e-3,
        f"worst |f| {worst * 1e3:.4f} mm",
    )
    expect(
        "the conditioned volume is within 0.5 % of the analytic value",
        abs(drift(after_project.signed_volume, REFERENCE_VOLUME_M3)) < 0.5,
        f"{drift(after_project.signed_volume, REFERENCE_VOLUME_M3):+.3f} %",
    )
    expect(
        "the conditioned area is within 1 % of the analytic value",
        abs(drift(after_project.surface_area, REFERENCE_AREA_M2)) < 1.0,
        f"{drift(after_project.surface_area, REFERENCE_AREA_M2):+.3f} %",
    )

    print("4. decimating the conditioned surface (recommended order)")
    for target in targets:
        started = time.perf_counter()
        reduced, stats = decimate(projected, target_faces=target)
        seconds = time.perf_counter() - started
        row(f"-> {target}", reduced, seconds)
        measured = measure_mesh(reduced)
        print(
            f"                       {stats['collapsed']} collapses, "
            f"{stats['rejected_link']} rejected by the link condition, "
            f"{stats['rejected_flip']} rejected for flipping a triangle"
        )
        expect(f"{target} faces stay watertight", measured.watertight)
        expect(
            f"{target} faces keep every sliver gone",
            measured.sliver_share == 0.0,
            f"{100 * measured.sliver_share:.4f} %",
        )
        expect(
            f"{target} faces keep the volume within 0.5 %",
            abs(drift(measured.signed_volume, REFERENCE_VOLUME_M3)) < 0.5,
            f"{drift(measured.signed_volume, REFERENCE_VOLUME_M3):+.3f} %",
        )
        if output is not None:
            write_obj(reduced, output / f"aff8_{target}.obj")

    print("5. the same budget in the other order (decimate first, then condition)")
    print(
        f"     {'stage':>16s} {'tris':>8s} {'min ang':>9s} {'max ar':>9s} {'slivers%':>8s} "
        f"{'dvol %':>8s} {'darea %':>8s} {'s':>6s}"
    )
    for target in targets[:2]:
        started = time.perf_counter()
        conditioned = condition_mesh(
            incoming,
            field=SuboffFullField(),
            target_faces=target,
            smooth_iterations=20,
            projection_iterations=2,
        )
        seconds = time.perf_counter() - started
        row(f"condition {target}", conditioned.mesh, seconds)
        print(
            f"                       {conditioned.collapsed_edges} collapses, "
            f"max smoothing move {1e3 * conditioned.max_smoothing_move:.2f} mm, "
            f"max projection move {1e3 * conditioned.max_projection_move:.3f} mm; "
            f"gates passed: {conditioned.passed}"
        )
        expect(
            f"condition_mesh({target}) keeps the declared tolerances",
            conditioned.passed,
        )

    print("6. handoff: the record a solver reads, before and after")
    final_surface = fem_surface(
        projected,
        components={"hull": SuboffHullField(), "sail": SuboffSailField(), "fins": SuboffFinsField()},
    )
    incoming_record = sampled.handoff()
    outgoing = final_surface.to_dict()
    print(f"     {'field':>22s} {'input':>22s} {'output':>22s}")
    for key in (
        "units",
        "coordinate_frame",
        "watertight",
        "edge_manifold",
        "euler_characteristic",
        "triangles",
        "signed_volume_m3",
        "surface_area_m2",
    ):
        print(f"     {key:>22s} {incoming_record[key]!s:>22s} {outgoing[key]!s:>22s}")
    print(f"     {'mesh digest':>22s} {incoming_record['mesh_sha256'][:16]:>22s} {final_surface.digest[:16]:>22s}")
    print(
        "     regions: "
        + ", ".join(
            f"{stats.name} {stats.triangles} tris / {stats.area_m2:.4f} m^2"
            for stats in final_surface.regions
        )
    )
    print(
        f"     every centroid is within {1e3 * max(stats.max_assignment_distance_m for stats in final_surface.regions):.2f} mm "
        "of the analytic component it was assigned to"
    )
    expect(
        "the output keeps the Euler characteristic",
        outgoing["euler_characteristic"] == incoming_record["euler_characteristic"],
    )
    expect(
        "the output is watertight",
        outgoing["watertight"] and outgoing["edge_manifold"],
    )
    expect(
        "the digest changes, because the surface changed",
        incoming_record["mesh_sha256"] != final_surface.digest,
    )
    expect(
        "all three analytic components come back as patches",
        all(stats.triangles > 0 for stats in final_surface.regions),
    )

    if output is not None:
        output.mkdir(parents=True, exist_ok=True)
        write_obj(projected, output / "aff8_conditioned.obj")
        write_gmsh(final_surface, output / "aff8_conditioned.msh")
        write_fem_json(final_surface, output / "aff8_conditioned.json")
        summary = {
            "schema": "tensorcad.mesh.conditioning/1",
            "resolution": list(resolution),
            "reference_volume_m3": REFERENCE_VOLUME_M3,
            "reference_area_m2": REFERENCE_AREA_M2,
            "incoming": before.to_dict(),
            "smoothed": after_smooth.to_dict(),
            "projected": after_project.to_dict(),
            "max_smoothing_move_m": smoothing_move,
            "max_projection_move_m": projection_move,
            "fem": final_surface.to_dict(),
            "passed": not failures,
        }
        (output / "conditioning.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
        print(f"     wrote {output}/aff8_conditioned.obj, .msh, .json and conditioning.json")

    print()
    if failures:
        print(f"FAILED ({len(failures)}): " + "; ".join(failures))
        return 1
    print("all conditioning checks passed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
