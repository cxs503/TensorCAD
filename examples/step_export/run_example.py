"""STEP export (optional OCP) and the TensorFEM handoff, measured end to end.

Two deliveries, one script, because they answer two different questions that are easy
to confuse:

* **The handoff** (always runs, no extra required): the AFF-8 surface is written as a
  gmsh 2.2 file with one physical tag per analytic component and as a JSON manifest,
  both re-read and compared against the mesh they came from. This is what a
  finite-element solver consumes; nothing here needs OpenCASCADE.
* **The STEP file** (optional, needs ``pip install "tensorcad[ocp]"``): the same surface
  is sewn into a closed shell, turned into a solid, checked by OpenCASCADE's own
  validity analyser, written as STEP, read back, and the volume compared against the
  mesh's signed volume. The script prints the numbers it measured and, when the wheels
  are absent, says so and scores only the handoff part.

What STEP export here is *not*: a B-rep model of the hull. Sewing a marching-tet surface
produces one planar face per triangle -- the script prints that ratio so nobody has to
guess. A designer who wants fillets and analytic surfaces needs a feature model, not a
tessellation; a solver that wants a closed, valid, measured solid gets exactly that.

Acceptance checks (the script exits non-zero if any fails):

1. the mesh is watertight and edge-manifold before anything is exported;
2. the gmsh round trip preserves coordinates exactly and areas to rounding;
3. the JSON round trip preserves the digest;
4. the patch areas sum to the measured wetted area, and all three components are present;
5. with OCP: the shell has no free or multiple edges, closes, the solid passes
   ``BRepCheck_Analyzer``, the kernel volume equals the mesh volume to 1e-9 relative,
   and a re-import returns that volume.
"""

from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from tensorcad.adapters.fem import fem_surface, read_fem_json, read_gmsh, write_fem_json, write_gmsh
from tensorcad.adapters.ocp import have_ocp, read_step, write_step
from tensorcad.geometry_field import to_mesh
from tensorcad.mesh import total_signed_volume
from tensorcad.suboff_appendages import SuboffFinsField, SuboffFullField, SuboffHullField, SuboffSailField


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--resolution", default="72,9,13", help="cells per axis for the sampled surface")
    parser.add_argument("--output", default="", help="write .msh/.json/.step here")
    parser.add_argument("--length-scale", type=float, default=1000.0, help="STEP units (mm per model unit)")
    args = parser.parse_args(argv)

    resolution = tuple(int(part) for part in args.resolution.split(","))
    output = Path(args.output) if args.output else Path(".")
    output.mkdir(parents=True, exist_ok=True)
    failures: list[str] = []

    def expect(label: str, condition: bool, detail: str = "") -> None:
        print(f"  [{'PASS' if condition else 'FAIL'}] {label}{': ' + detail if detail else ''}")
        if not condition:
            failures.append(label)

    print("1. the surface, and the patches it splits into")
    started = time.perf_counter()
    sampled = to_mesh(SuboffFullField(), resolution, margin=2)
    print(
        f"     {len(sampled.mesh.vertices)} vertices, {len(sampled.mesh.faces)} triangles, "
        f"{time.perf_counter() - started:.2f} s; watertight: {sampled.report.watertight}"
    )
    started = time.perf_counter()
    surface = fem_surface(
        sampled,
        components={"hull": SuboffHullField(), "sail": SuboffSailField(), "fins": SuboffFinsField()},
    )
    print(f"     labelled by nearest analytic component in {time.perf_counter() - started:.2f} s")
    print(f"     {'patch':>8s} {'tris':>8s} {'area m^2':>10s} {'share %':>8s} {'max dist mm':>12s}")
    for stats in surface.regions:
        print(
            f"     {stats.name:>8s} {stats.triangles:8d} {stats.area_m2:10.5f} "
            f"{100.0 * stats.area_m2 / surface.area_m2:8.3f} "
            f"{1e3 * stats.max_assignment_distance_m:12.2f}"
        )
    print(
        f"     volume {surface.volume_m3:.6f} m^3, wetted area {surface.area_m2:.6f} m^2, "
        f"Euler characteristic {surface.euler_characteristic}"
    )
    expect("the surface is watertight and edge-manifold", surface.watertight and surface.edge_manifold)
    expect(
        "all three analytic components appear as patches",
        all(stats.triangles > 0 for stats in surface.regions),
    )
    expect(
        "the patch areas sum to the measured area",
        abs(sum(stats.area_m2 for stats in surface.regions) - surface.area_m2) < 1e-9,
    )
    expect(
        "the hull dominates the wetted area, as it must",
        100.0 * surface.regions[0].area_m2 / surface.area_m2 > 80.0,
        f"{100.0 * surface.regions[0].area_m2 / surface.area_m2:.2f} %",
    )

    print("2. gmsh 2.2 round trip (no optional dependency)")
    gmsh_path = write_gmsh(surface, output / "aff8_surface.msh")
    parsed = read_gmsh(gmsh_path)
    print(
        f"     {gmsh_path.stat().st_size / 1e6:.2f} MB, {parsed['nodes']} nodes, "
        f"{parsed['triangles']} triangles, physical names "
        + ", ".join(f"{tag}:{name}" for tag, name in parsed["physical_names"].items())
    )
    worst_node = max(
        abs(a - b) for node_id, node in enumerate(surface.nodes, start=1)
        for a, b in zip(parsed["node_table"][node_id], node, strict=True)
    )
    worst_area = max(
        abs(parsed["area_by_region_m2"][stats.physical_id] - stats.area_m2) for stats in surface.regions
    )
    print(f"     worst node displacement {worst_node:.1e} m, worst per-patch area difference {worst_area:.1e} m^2")
    expect("the node table survives the round trip exactly", worst_node == 0.0)
    expect("every patch area survives the round trip", worst_area < 1e-9, f"{worst_area:.1e} m^2")

    print("3. JSON round trip")
    json_path = write_fem_json(surface, output / "aff8_surface.json")
    restored = read_fem_json(json_path)
    print(f"     {json_path.stat().st_size / 1e6:.2f} MB, digest {restored.digest[:16]}...")
    expect("the JSON digest matches the surface", restored.digest == surface.digest)
    expect("the JSON keeps the source digest", restored.source_digest == surface.source_digest)
    expect("the JSON keeps the topology flags", restored.watertight == surface.watertight)

    print("4. STEP export" + ("" if have_ocp() else " (skipped: OCP is not installed)"))
    if not have_ocp():
        print('     install the extra to run this part: pip install "tensorcad[ocp]"')
        print("     nothing below is simulated; the handoff above is the deliverable without OCP")
    else:
        step_path = output / "aff8_solid.step"
        started = time.perf_counter()
        record = write_step(sampled.mesh, step_path, length_scale=args.length_scale, round_trip=True)
        seconds = time.perf_counter() - started
        mesh_volume = total_signed_volume(sampled.mesh)
        print(
            f"     sewn in {record.seconds_sew:.1f} s, written in {record.seconds_write:.1f} s "
            f"({seconds:.1f} s total); {record.file_bytes / 1e6:.2f} MB"
        )
        print(
            f"     free edges {record.free_edges}, multiple edges {record.multiple_edges}, "
            f"shells {record.shells}, solids {record.solids}, "
            f"shell closed {record.shell_closed}, solid valid {record.solid_valid}"
        )
        print(
            f"     kernel volume {record.volume_m3:.9f} m^3 vs mesh signed volume "
            f"{mesh_volume:.9f} m^3 (rel {record.volume_m3 / mesh_volume - 1.0:+.2e})"
        )
        print(
            f"     re-imported volume {record.round_trip_volume_m3:.9f} m^3 "
            f"(rel {record.round_trip_volume_error:+.2e})"
        )
        imported = read_step(step_path, length_scale=args.length_scale)
        print(
            f"     STEP contains {imported.solids} solid(s), {imported.faces} faces, "
            f"valid {imported.valid}; one planar face per triangle, as sewing implies"
        )
        expect("the sewn shell has no free edges", record.free_edges == 0)
        expect("the sewn shell has no multiple edges", record.multiple_edges == 0)
        expect("the shell is closed", record.shell_closed)
        expect("OpenCASCADE accepts the solid", record.solid_valid)
        expect(
            "the kernel volume equals the mesh volume",
            abs(record.volume_m3 / mesh_volume - 1.0) < 1e-9,
            f"{record.volume_m3 / mesh_volume - 1.0:+.2e}",
        )
        expect(
            "a re-import returns the same volume",
            abs(record.round_trip_volume_error or 0.0) < 1e-9,
            f"{record.round_trip_volume_error:+.2e}",
        )
        expect(
            "every closed shell became a solid in the file",
            imported.solids == record.shells,
            f"{imported.solids} solid(s) for {record.shells} shell(s)",
        )
        expect(
            "the STEP face count is the triangle count, i.e. this is a faceted solid",
            imported.faces == len(sampled.mesh.faces),
            f"{imported.faces} faces for {len(sampled.mesh.faces)} triangles",
        )
        if record.shells > 1:
            print(
                "     note: the input mesh has more than one closed island (a marching-tet "
                "artefact at this resolution); each island is a solid, none is dropped"
            )

    print()
    if failures:
        print(f"FAILED ({len(failures)}): " + "; ".join(failures))
        return 1
    print(f"all checks passed; artifacts in {output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
