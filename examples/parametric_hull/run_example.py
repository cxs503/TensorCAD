#!/usr/bin/env python3
"""Run the first TensorCAD geometry benchmark and *verify* it numerically.

This is not a screenshot test. It generates a parametric hull, checks that the
exported surface is a closed oriented solid, compares the mesh-integrated
displacement volume against the analytic value of the generating functions, runs
a resolution study, and exits non-zero if any claim fails. The printed table is
measured output from the code in this repository.

Usage::

    python examples/parametric_hull/run_example.py [--output results/hull-demo]
"""
from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
if str(REPO_ROOT) not in sys.path:  # allow running from a checkout without installing
    sys.path.insert(0, str(REPO_ROOT))

import itertools

from tensorcad import HullParameters, build_mesh, compute_hydrostatics, export_hull, generate_hull
from tensorcad import mesh as meshlib
from tensorcad.hull import TAPER_INTEGRAL, analytic_volume


def analytic_waterplane_area(length: float, beam: float) -> float:
    return beam * length * TAPER_INTEGRAL


def resolution_study(length: float, beam: float, draft: float, grids: list[tuple[int, int]]) -> list[dict[str, float]]:
    rows: list[dict[str, float]] = []
    for nx, nz in grids:
        model_mesh, _offsets = build_mesh(HullParameters(length, beam, draft, nx, nz))
        hydro = compute_hydrostatics(model_mesh, model_mesh.bounding_box().maximum[2], closed=True)
        rows.append(
            {
                "nx": nx,
                "nz": nz,
                "volume": hydro.displacement_volume,
                "volume_error_percent": 100.0 * (hydro.displacement_volume / analytic_volume(length, beam, draft) - 1.0),
                "waterplane_area": hydro.waterplane_area,
                "waterplane_error_percent": 100.0
                * (hydro.waterplane_area / analytic_waterplane_area(length, beam) - 1.0),
                "wetted_area": hydro.wetted_surface_area,
            }
        )
    return rows


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--length", type=float, default=10.0)
    parser.add_argument("--beam", type=float, default=2.0)
    parser.add_argument("--draft", type=float, default=0.8)
    parser.add_argument("--output", default="results/hull-demo", help="directory for exported artifacts")
    args = parser.parse_args(argv)

    failures: list[str] = []

    def expect(label: str, condition: bool, detail: str = "") -> None:
        status = "PASS" if condition else "FAIL"
        print(f"  [{status}] {label}{': ' + detail if detail else ''}")
        if not condition:
            failures.append(label)

    parameters = HullParameters(args.length, args.beam, args.draft, 41, 13)
    print(f"1. generating hull for {json.dumps(parameters.as_dict(), sort_keys=True)}")
    mesh, offsets = build_mesh(parameters)

    box = mesh.bounding_box()
    print("2. closed-solid checks")
    expect("watertight (every edge used by exactly two faces)", meshlib.is_watertight(mesh))
    expect("edge-manifold", meshlib.is_edge_manifold(mesh))
    expect("no zero-area faces", not meshlib.zero_area_faces(mesh))
    expect("no boundary edges", not meshlib.boundary_edges(mesh))
    expect(
        "outward orientation (signed volume > 0)",
        meshlib.total_signed_volume(mesh) > 0,
        f"volume = {meshlib.total_signed_volume(mesh):.6f} m^3",
    )
    euler = len(mesh.vertices) - len(meshlib.edge_face_counts(mesh)) + len(mesh.faces)
    expect("Euler characteristic V-E+F = 2", euler == 2, f"V={len(mesh.vertices)} F={len(mesh.faces)} -> {euler}")

    print("3. dimensional fidelity")
    expect("length matches request", math.isclose(box.sizes[0], args.length, rel_tol=1e-12), f"{box.sizes[0]:.9g}")
    expect("draft matches request", math.isclose(box.sizes[2], args.draft, rel_tol=1e-12), f"{box.sizes[2]:.9g}")
    expect(
        "beam is realised on an odd station grid",
        math.isclose(box.sizes[1], args.beam, rel_tol=1e-12),
        f"{box.sizes[1]:.9g}",
    )
    expect("waterline plane is exactly z=0", box.maximum[2] == 0.0)
    expect(
        "offsets: zero breadth at both ends, one row per station",
        len(offsets) == parameters.nx and offsets[0][1] == 0.0 and offsets[-1][1] == 0.0,
    )

    print("4. hydrostatics vs analytic value of the sampled surface")
    hydro = compute_hydrostatics(mesh, box.maximum[2], closed=True)
    exact = analytic_volume(args.length, args.beam, args.draft)
    print(f"     mesh volume       = {hydro.displacement_volume:.6f} m^3")
    print(f"     analytic volume     = {exact:.6f} m^3")
    print(f"     wetted surface    = {hydro.wetted_surface_area:.6f} m^2")
    print(f"     waterplane area   = {hydro.waterplane_area:.6f} m^2")
    print(f"     LCB / LCF / VCB   = {hydro.longitudinal_centre_buoyancy:.4f} / "
          f"{hydro.longitudinal_centre_flotation:.4f} / {hydro.vertical_centre_buoyancy:.4f} m")
    print(f"     Cb / Cp           = {hydro.block_coefficient:.4f} / {hydro.prismatic_coefficient:.4f}")
    expect("volume within 1% of analytic", abs(hydro.displacement_volume / exact - 1.0) < 0.01)
    expect("LCB on the centre of flotation (form is symmetric)",
           math.isclose(hydro.longitudinal_centre_buoyancy, hydro.longitudinal_centre_flotation, rel_tol=1e-9))
    expect("Cb consistent with volume", math.isclose(hydro.block_coefficient, exact / (args.length * args.beam * args.draft), rel_tol=2e-2))

    print("5. resolution study (error must shrink as the grid refines)")
    rows = resolution_study(args.length, args.beam, args.draft, [(11, 5), (21, 9), (41, 13), (81, 21), (161, 41)])
    print("       nx   nz     volume [m^3]    V err %     Awp [m^2]    Awp err %")
    for row in rows:
        print(f"     {row['nx']:4d} {row['nz']:3d}   {row['volume']:11.6f}  {row['volume_error_percent']:+7.3f}   "
              f"{row['waterplane_area']:10.6f}  {row['waterplane_error_percent']:+7.3f}")
    errors = [abs(row["volume_error_percent"]) for row in rows]
    expect("volume error decreases monotonically with refinement", all(a > b for a, b in itertools.pairwise(errors)))
    expect("finest grid within 0.2% of analytic", errors[-1] < 0.2)

    print("6. export artifacts")
    model = generate_hull(parameters)
    artifacts = export_hull(model, args.output, convergence=True)
    for name, path in artifacts.items():
        size = Path(path).stat().st_size
        expect(f"{name} artifact written ({size} bytes)", size > 0)
    written = json.loads(Path(artifacts["metadata"]).read_text())
    expect("metadata records a measured convergence estimate",
           written.get("discretisation", {}).get("status") == "measured")

    print()
    if failures:
        print(f"FAILED {len(failures)} check(s): {', '.join(failures)}")
        return 1
    print("all checks passed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
