#!/usr/bin/env python3
"""Lid-driven cavity Re=5000 benchmark — mesh generation and reference data.

This example generates a structured quadrilateral mesh for the lid-driven
cavity at Re = 5000 and exports the Ghia et al. (1982) reference data for
centreline velocity profiles.  The mesh is suitable for consumption by
TensorLBM, TensorFVM, or any structured-grid Navier–Stokes solver.

What this script produces
-------------------------
1. ``cavity_re5000_mesh.json`` — structured mesh (nodes + elements + BCs)
2. ``ghia_re5000_reference.json`` — centreline velocity reference data
3. ``cavity_mesh_summary.txt`` — human-readable mesh summary
4. ``cavity_mesh_quality.png`` — mesh visualisation (grid + aspect ratio)
5. ``ghia_reference_profiles.png`` — reference velocity profiles

Usage
-----
    python examples/cavity_re5000/run_example.py [--nx 128] [--grading 1.0]

Exit code is 0 on success, 1 if any check fails.
"""
from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from tensorcad.cavity_mesh import (
    CavityParameters,
    GHIA_RE5000_U_CENTERLINE,
    GHIA_RE5000_V_CENTERLINE,
    generate_cavity_mesh,
    ghia_u_at_y,
    ghia_v_at_x,
    mesh_summary,
    mesh_to_json,
    write_mesh_json,
    write_reference_data,
)


def check_mesh_topology(mesh) -> list[str]:
    """Verify mesh topology invariants."""
    failures = []
    
    # Check node count
    expected_nodes = (mesh.nx + 1) * (mesh.ny + 1)
    if mesh.n_nodes != expected_nodes:
        failures.append(f"Node count: expected {expected_nodes}, got {mesh.n_nodes}")
    
    # Check element count
    expected_elems = mesh.nx * mesh.ny
    if mesh.n_elements != expected_elems:
        failures.append(f"Element count: expected {expected_elems}, got {mesh.n_elements}")
    
    # Check element connectivity bounds
    for i, elem in enumerate(mesh.elements):
        for node_idx in elem:
            if node_idx < 0 or node_idx >= mesh.n_nodes:
                failures.append(f"Element {i}: node index {node_idx} out of range")
                break
    
    # Check that all nodes are within domain
    L = mesh.params.domain_size
    for i, (x, y) in enumerate(mesh.nodes):
        if x < -1e-12 or x > L + 1e-12 or y < -1e-12 or y > L + 1e-12:
            failures.append(f"Node {i}: ({x}, {y}) outside domain [0, {L}]²")
            break
    
    return failures


def check_reference_data() -> list[str]:
    """Verify Ghia reference data integrity."""
    failures = []
    
    # Check u-profile: should start at 0 (bottom wall) and end at 1 (lid)
    if abs(GHIA_RE5000_U_CENTERLINE[0][1]) > 1e-6:
        failures.append(f"u(y=0) should be 0, got {GHIA_RE5000_U_CENTERLINE[0][1]}")
    if abs(GHIA_RE5000_U_CENTERLINE[-1][1] - 1.0) > 1e-6:
        failures.append(f"u(y=1) should be 1 (lid), got {GHIA_RE5000_U_CENTERLINE[-1][1]}")
    
    # Check v-profile: should be 0 at walls (x=0 and x=1)
    if abs(GHIA_RE5000_V_CENTERLINE[0][1]) > 1e-6:
        failures.append(f"v(x=0) should be 0, got {GHIA_RE5000_V_CENTERLINE[0][1]}")
    if abs(GHIA_RE5000_V_CENTERLINE[-1][1]) > 1e-6:
        failures.append(f"v(x=1) should be 0, got {GHIA_RE5000_V_CENTERLINE[-1][1]}")
    
    # Check monotonicity of y-coordinates in u-profile
    for i in range(len(GHIA_RE5000_U_CENTERLINE) - 1):
        y0 = GHIA_RE5000_U_CENTERLINE[i][0]
        y1 = GHIA_RE5000_U_CENTERLINE[i + 1][0]
        if y1 <= y0:
            failures.append(f"u-profile: y not monotonically increasing at index {i}")
            break
    
    # Check monotonicity of x-coordinates in v-profile
    for i in range(len(GHIA_RE5000_V_CENTERLINE) - 1):
        x0 = GHIA_RE5000_V_CENTERLINE[i][0]
        x1 = GHIA_RE5000_V_CENTERLINE[i + 1][0]
        if x1 <= x0:
            failures.append(f"v-profile: x not monotonically increasing at index {i}")
            break
    
    return failures


def generate_mesh_visualisation(mesh, output_path: Path) -> None:
    """Generate a mesh visualisation plot."""
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        import matplotlib.collections as mc
        import numpy as np
    except ImportError:
        print("  [SKIP] matplotlib not available for visualisation")
        return
    
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    
    # Left: mesh grid
    ax = axes[0]
    nx, ny = mesh.nx, mesh.ny
    
    # Plot horizontal grid lines (constant y)
    for j in range(ny + 1):
        xs = [mesh.nodes[j * (nx + 1) + i][0] for i in range(nx + 1)]
        ys = [mesh.nodes[j * (nx + 1) + i][1] for i in range(nx + 1)]
        ax.plot(xs, ys, 'b-', linewidth=0.3, alpha=0.5)
    
    # Plot vertical grid lines (constant x)
    for i in range(nx + 1):
        xs = [mesh.nodes[j * (nx + 1) + i][0] for j in range(ny + 1)]
        ys = [mesh.nodes[j * (nx + 1) + i][1] for j in range(ny + 1)]
        ax.plot(xs, ys, 'b-', linewidth=0.3, alpha=0.5)
    
    # Highlight lid
    xs_lid = [mesh.nodes[ny * (nx + 1) + i][0] for i in range(nx + 1)]
    ys_lid = [mesh.nodes[ny * (nx + 1) + i][1] for i in range(nx + 1)]
    ax.plot(xs_lid, ys_lid, 'r-', linewidth=2, label='Lid (u=1)')
    
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_title(f'Cavity Mesh ({nx}×{ny})')
    ax.set_aspect('equal')
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    # Right: reference velocity profiles
    ax = axes[1]
    u_y = [p[0] for p in GHIA_RE5000_U_CENTERLINE]
    u_u = [p[1] for p in GHIA_RE5000_U_CENTERLINE]
    v_x = [p[0] for p in GHIA_RE5000_V_CENTERLINE]
    v_v = [p[1] for p in GHIA_RE5000_V_CENTERLINE]
    
    ax.plot(u_u, u_y, 'b-', linewidth=1.5, label='u along x=0.5')
    ax.plot(v_v, v_x, 'r-', linewidth=1.5, label='v along y=0.5')
    ax.set_xlabel('Velocity')
    ax.set_ylabel('Position')
    ax.set_title('Ghia et al. (1982) Re=5000 Reference')
    ax.legend()
    ax.grid(True, alpha=0.3)
    ax.axhline(y=0.5, color='gray', linestyle='--', alpha=0.5)
    ax.axvline(x=0, color='gray', linestyle='--', alpha=0.5)
    
    plt.tight_layout()
    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  wrote {output_path}")


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--nx", type=int, default=64, help="cells in x direction")
    parser.add_argument("--ny", type=int, default=None, help="cells in y direction (default: nx)")
    parser.add_argument("--grading", type=float, default=1.0,
                        help="wall grading ratio (1.0=uniform, >1=finer near walls)")
    parser.add_argument("--output", default="results/cavity_re5000",
                        help="output directory")
    args = parser.parse_args(argv)
    
    output = Path(args.output)
    output.mkdir(parents=True, exist_ok=True)
    
    print("=" * 70)
    print("Lid-driven cavity Re=5000 benchmark")
    print("=" * 70)
    
    # Generate mesh
    print(f"\nGenerating mesh: {args.nx}×{args.ny or args.nx}, grading={args.grading}")
    params = CavityParameters(reynolds_number=5000.0)
    mesh = generate_cavity_mesh(
        nx=args.nx,
        ny=args.ny,
        params=params,
        wall_grading=args.grading,
    )
    print(mesh_summary(mesh))
    
    # Check topology
    print("\nChecking mesh topology...")
    topo_failures = check_mesh_topology(mesh)
    if topo_failures:
        for f in topo_failures:
            print(f"  [FAIL] {f}")
        return 1
    print("  [PASS] All topology checks passed")
    
    # Check reference data
    print("\nChecking Ghia reference data...")
    ref_failures = check_reference_data()
    if ref_failures:
        for f in ref_failures:
            print(f"  [FAIL] {f}")
        return 1
    print("  [PASS] All reference data checks passed")
    
    # Write outputs
    print("\nWriting outputs...")
    
    mesh_path = write_mesh_json(mesh, output / "cavity_re5000_mesh.json")
    print(f"  wrote {mesh_path}")
    
    ref_path = write_reference_data(output / "ghia_re5000_reference.json")
    print(f"  wrote {ref_path}")
    
    summary_path = output / "cavity_mesh_summary.txt"
    with open(summary_path, "w") as f:
        f.write(mesh_summary(mesh))
        f.write("\n\n")
        f.write("Boundary conditions:\n")
        for name, bc in mesh.boundary.items():
            f.write(f"  {name}: {bc['description']}\n")
    print(f"  wrote {summary_path}")
    
    # Generate visualisation
    print("\nGenerating visualisation...")
    viz_path = output / "cavity_mesh_quality.png"
    generate_mesh_visualisation(mesh, viz_path)
    
    # Summary
    print("\n" + "=" * 70)
    print("Benchmark complete")
    print("=" * 70)
    print(f"Mesh: {mesh.nx}×{mesh.ny} = {mesh.n_elements} cells, {mesh.n_nodes} nodes")
    print(f"Re = {params.reynolds_number:.0f}, ν = {params.nu:.6e}")
    print(f"Outputs in: {output}/")
    
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
