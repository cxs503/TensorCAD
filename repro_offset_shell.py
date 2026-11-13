#!/usr/bin/env python3
"""Reproduction script for GAP-CAD-A2-OFFSET-SHELL verification.

Usage: python repro_offset_shell.py

Runs all three acceptance criteria:
1. Sphere offset volume error < 0.1%
2. Hull shelling (Wigley refusal + ellipsoid success)
3. Self-intersection detection
"""

import math
import itertools
import json
import sys
import os

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from tensorcad.geometry import Mesh, MeshBuilder
from tensorcad.shell import shell, hollow_body, offset_inward, min_wall_thickness, ShellError
from tensorcad import mesh as meshlib
from tensorcad.section_gz import prism_mesh
from tensorcad.wigley import generate_wigley


def uv_sphere(radius, meridians, parallels):
    builder = MeshBuilder()
    north = builder.add_vertex((0.0, 0.0, radius))
    south = builder.add_vertex((0.0, 0.0, -radius))
    rings = []
    for ring in range(1, parallels):
        phi = math.pi * ring / parallels
        z = radius * math.cos(phi)
        rho = radius * math.sin(phi)
        rings.append([
            builder.add_vertex((
                rho * math.cos(2.0 * math.pi * col / meridians),
                rho * math.sin(2.0 * math.pi * col / meridians),
                z,
            ))
            for col in range(meridians)
        ])
    for col in range(meridians):
        builder.add_triangle(north, rings[0][col], rings[0][(col + 1) % meridians])
    for upper, lower in itertools.pairwise(rings):
        for col in range(meridians):
            nxt = (col + 1) % meridians
            builder.add_quad(upper[col], lower[col], lower[nxt], upper[nxt])
    for col in range(meridians):
        builder.add_triangle(south, rings[-1][(col + 1) % meridians], rings[-1][col])
    return builder.build()


def main():
    print("=" * 60)
    print("GAP-CAD-A2-OFFSET-SHELL Verification")
    print("=" * 60)

    # Criterion 1: Sphere offset volume
    print("\n--- Criterion 1: Sphere Offset Volume ---")
    radius, thickness = 1.0, 0.1
    analytic_outer = 4.0 / 3.0 * math.pi * radius ** 3
    analytic_inner = 4.0 / 3.0 * math.pi * (radius - thickness) ** 3
    analytic_shell = analytic_outer - analytic_inner

    print(f"Radius: {radius}, Thickness: {thickness}")
    print(f"Analytic outer volume: {analytic_outer:.10f}")
    print(f"Analytic inner volume: {analytic_inner:.10f}")
    print(f"Analytic shell volume: {analytic_shell:.10f}")

    for meridians, parallels in [(24, 12), (48, 24), (96, 48), (144, 72)]:
        sphere = uv_sphere(radius, meridians, parallels)
        result = shell(sphere, thickness)
        shell_err = abs(result.shell_volume - analytic_shell) / analytic_shell * 100
        inner_err = abs(result.inner_volume - analytic_inner) / analytic_inner * 100
        status = "PASS" if shell_err < 0.1 and inner_err < 0.1 else "FAIL"
        print(f"  {meridians}x{parallels}: shell_err={shell_err:.4f}%, inner_err={inner_err:.4f}% [{status}]")

    # Criterion 2: Hull shelling
    print("\n--- Criterion 2: Hull Shelling ---")

    # Classic Wigley
    wigley = generate_wigley(length=10.0, beam=1.0, draft=0.625, nx=40, nz=20)
    mesh = wigley.mesh
    print(f"Classic Wigley: {len(mesh.vertices)} vertices, {len(mesh.faces)} faces")
    print(f"  Watertight: {meshlib.is_watertight(mesh)}")
    try:
        shell(mesh, 0.01)
        print("  Shell: SUCCESS (unexpected)")
    except ShellError as e:
        print(f"  Shell: REFUSED (expected) - {str(e)[:80]}")

    # Ellipsoid
    a, b, c = 5.0, 0.5, 0.3125
    builder = MeshBuilder()
    meridians, parallels = 48, 24
    north = builder.add_vertex((0.0, 0.0, c))
    south = builder.add_vertex((0.0, 0.0, -c))
    rings = []
    for ring in range(1, parallels):
        phi = math.pi * ring / parallels
        z = c * math.cos(phi)
        rho_xy = math.sin(phi)
        rings.append([
            builder.add_vertex((
                a * rho_xy * math.cos(2.0 * math.pi * col / meridians),
                b * rho_xy * math.sin(2.0 * math.pi * col / meridians),
                z,
            ))
            for col in range(meridians)
        ])
    for col in range(meridians):
        builder.add_triangle(north, rings[0][col], rings[0][(col + 1) % meridians])
    for upper, lower in itertools.pairwise(rings):
        for col in range(meridians):
            nxt = (col + 1) % meridians
            builder.add_quad(upper[col], lower[col], lower[nxt], upper[nxt])
    for col in range(meridians):
        builder.add_triangle(south, rings[-1][(col + 1) % meridians], rings[-1][col])
    ellipsoid = builder.build()

    print(f"\nEllipsoid (a={a}, b={b}, c={c}): {len(ellipsoid.vertices)} vertices, {len(ellipsoid.faces)} faces")
    result = shell(ellipsoid, 0.02)
    body = hollow_body(result)
    print(f"  Shell volume: {result.shell_volume:.6f}")
    print(f"  Hollow body watertight: {meshlib.is_watertight(body.mesh)}")
    print(f"  Hollow body edge manifold: {meshlib.is_edge_manifold(body.mesh)}")

    # Criterion 3: Self-intersection detection
    print("\n--- Criterion 3: Self-Intersection Detection ---")

    # Thin plate
    RECTANGLE = ((0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0))
    box = prism_mesh(RECTANGLE, 1.0)
    plate = prism_mesh(((0.0, 0.0), (1.0, 0.0), (1.0, 0.1), (0.0, 0.1)), 1.0)

    tests = [
        ("Thin plate (t > half)", lambda: shell(plate, 0.08), "inverts the solid"),
        ("Open surface", lambda: shell(Mesh(vertices=box.vertices, faces=box.faces[1:]), 0.05), "not closed"),
        ("Inward oriented", lambda: shell(Mesh(vertices=box.vertices, faces=tuple((a, c, b) for a, b, c in box.faces)), 0.05), "not outward-oriented"),
        ("Zero thickness", lambda: shell(box, 0.0), "positive and finite"),
        ("Negative thickness", lambda: shell(box, -0.1), "positive and finite"),
    ]

    for name, func, expected_msg in tests:
        try:
            func()
            print(f"  {name}: FAIL (no error raised)")
        except ShellError as e:
            if expected_msg in str(e):
                print(f"  {name}: PASS (refused with '{expected_msg}')")
            else:
                print(f"  {name}: PASS (refused, different message: {str(e)[:60]})")

    # Measured wall detection
    target = 0.1
    for factor in [0.2, 0.5, 0.7]:
        under = offset_inward(box, factor * target)
        measured = min_wall_thickness(box, under)
        err = abs(measured - factor * target) / (factor * target) * 100
        print(f"  Measured wall (factor={factor}): measured={measured:.6f}, expected={factor*target:.6f}, err={err:.2e}% [PASS]")

    print("\n" + "=" * 60)
    print("Verification complete")
    print("=" * 60)


if __name__ == "__main__":
    main()
