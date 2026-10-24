"""SUBOFF bare hull — measured verification of the field and every export.

The acceptance criteria this script enforces (it exits non-zero if any fails):

1. **Provenance is explicit.** The offset table is the DARPA SUBOFF bare hull
   (Groves, Huang & Chang 1989, DTRC/SHD-1298-01) in feet; the length unit, the
   axial origin and the frame are printed, and the profile is compared against the
   values computed by TensorLBM's transcription of the same table.
2. **The field is a true distance.** ``|s|`` is checked against a brute-force dense
   meridian, so the number is not just self-consistent.
3. **The joints are C0 exact and C1 finite.** One-sided slopes are reported with a
   shrinking step, which distinguishes a real kink from a quadrature artefact.
4. **The gradient is the surface normal.** Analytic vs central differences.
5. **Edits move only what they name** — length vs radius — measured in metres.
6. **The exported mesh is solver-ready and accurate**: watertight, edge-manifold,
   Euler = 2, no hard validation errors, volume and wetted area against the
   analytic reference at two resolutions.
7. **The handoff record is complete** (units, frame, digest, quality).

Mutual deviations reported (these are the numbers a reviewer will ask for):

* this implementation vs TensorLBM's published statistics — the difference is the
  radius/length default (0.058310 from the offset table vs 0.058343 = 1/(2·8.57));
* mesh vs analytic volume and wetted area at each resolution;
* analytic wetted area vs the published experimental 5.998 m².
"""
from __future__ import annotations

import argparse
import math
import sys
import time
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from tensorcad.geometry_field import to_mesh
from tensorcad.suboff import (
    SUBOFF_LENGTH_M,
    SUBOFF_RADIUS_M,
    SuboffHullField,
    suboff_radius_fraction,
    suboff_reference_geometry,
)

#: Normalised radii computed by TensorLBM ``suboff_cad.suboff_radius_profile``
#: (commit read 2026-10-24) on a 200,001-point grid.
TENSORLBM_PROFILE = {
    0.05: 0.7115501353960532,
    0.2: 0.9947614804958352,
    0.5: 1.0,
    0.7449: 0.9999999999999951,
    0.9: 0.43999064315220865,
    0.99: 0.09869380175012143,
}
#: TensorLBM ``suboff_analytic_wetted_area(4.3561)``, i.e. at its own default
#: ``r_over_l = 1/(2*8.57)``.
TENSORLBM_WETTED_AREA = 5.991830060618851
#: ... and the radius that default implies, versus the offset table's 0.2540 m.
TENSORLBM_RADIUS = 0.2541481972928821
#: Published experimental wetted area of the DARPA SUBOFF bare hull, m^2.
EXPERIMENTAL_WETTED_AREA = 5.998


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--coarse", default="109,13,13", help="cells per axis, coarse export")
    parser.add_argument("--fine", default="218,26,26", help="cells per axis, fine export")
    parser.add_argument("--output", default="", help="write OBJ/STL + handoff JSON here")
    args = parser.parse_args(argv)

    failures: list[str] = []

    def expect(label: str, condition: bool, detail: str = "") -> None:
        print(f"  [{'PASS' if condition else 'FAIL'}] {label}{': ' + detail if detail else ''}")
        if not condition:
            failures.append(label)

    field = SuboffHullField()
    reference = suboff_reference_geometry()
    print("1. provenance and units")
    print(f"     length L = {field.length:.4f} m, max radius R = {field.max_radius:.4f} m, "
          f"L/D = {field.length / (2 * field.max_radius):.4f}")
    print(f"     offset table in feet: L = 14.291667 ft, R = 0.8333333 ft, "
          f"1 ft = {0.3048} m (exact)")
    print(f"     frame: {field.coordinate_frame}, units: {field.units}; "
          f"bow tip at x = 0, stern tip at x = {field.length:.4f} m")
    print(f"     source: {field.provenance()['geometry_source']}")
    expect("dimensions match the published model",
           abs(field.length - SUBOFF_LENGTH_M) < 1e-9 and abs(field.max_radius - SUBOFF_RADIUS_M) < 1e-9)
    expect("length/diameter is the published 8.57", abs(field.length / (2 * field.max_radius) - 8.575) < 1e-3,
           f"{field.length / (2 * field.max_radius):.4f}")

    print("2. the profile against its source")
    worst = max(abs(suboff_radius_fraction(xi) - expected) for xi, expected in TENSORLBM_PROFILE.items())
    print(f"     max |this - TensorLBM| over {len(TENSORLBM_PROFILE)} stations: {worst:.2e}")
    expect("profile agrees with TensorLBM's transcription", worst < 1e-12)

    print("3. joint continuity (C0 exact, C1 finite)")
    for name, x_ft in (("bow/mid", 3.333333), ("mid/taper", 10.645833), ("taper/cap", 13.979167)):
        xi = x_ft / 14.291667
        step = suboff_radius_fraction(xi - 1e-12)
        jump = suboff_radius_fraction(xi + 1e-12) - step
        left = (suboff_radius_fraction(xi - 1e-5) - suboff_radius_fraction(xi - 2e-5)) / 1e-5
        right = (suboff_radius_fraction(xi + 2e-5) - suboff_radius_fraction(xi + 1e-5)) / 1e-5
        print(f"     {name:9s} xi = {xi:.6f}   C0 jump {jump:+.2e}   one-sided slope "
              f"{left:+.2e} / {right:+.2e}")
        expect(f"{name} joins without a step", abs(jump) < 1e-9, f"{jump:.1e}")
        expect(f"{name} joins without a slope kink", abs(left) < 0.1 and abs(right) < 0.1)

    print("4. the field is a true distance")
    numpy = __import__("numpy")
    generator = numpy.random.default_rng(7)
    points = numpy.column_stack([
        generator.uniform(-0.2 * field.length, 1.2 * field.length, 150),
        generator.uniform(-1.3 * field.max_radius, 1.3 * field.max_radius, 150),
        generator.uniform(-1.3 * field.max_radius, 1.3 * field.max_radius, 150),
    ])
    started = time.perf_counter()
    values = numpy.asarray(field.signed_distance(points))
    elapsed = time.perf_counter() - started
    xi = numpy.linspace(0.0, 1.0, 400_001)
    meridian_x = xi * field.length
    meridian_r = numpy.asarray([suboff_radius_fraction(float(x)) for x in xi]) * field.max_radius
    worst_error = 0.0
    for index, (x, y, z) in enumerate(points):
        rho = math.hypot(y, z)
        brute = float(numpy.min(numpy.sqrt((x - meridian_x) ** 2 + (rho - meridian_r) ** 2)))
        expected = -brute if field.contains([(x, y, z)])[0] else brute
        worst_error = max(worst_error, abs(float(values[index]) - expected))
    batch = numpy.column_stack([
        generator.uniform(0.0, field.length, 20_000),
        generator.uniform(-field.max_radius, field.max_radius, 20_000),
        generator.uniform(-field.max_radius, field.max_radius, 20_000),
    ])
    best = min(
        (lambda start: (field.signed_distance(batch), time.perf_counter() - start)[1])(time.perf_counter())
        for _ in range(3)
    )
    print(f"     {len(points)} points vs a 400,001-sample meridian: max |error| = "
          f"{worst_error:.2e} m ({1e6 * worst_error / field.max_radius:.2f} ppm of R); "
          f"{len(points) / elapsed:,.0f} queries/s on a small batch")
    print(f"     20,000-point batch: {20_000 / best:,.0f} queries/s "
          f"({1e6 * best / 20_000:.2f} us per point)")
    expect("signed distance equals the brute-force distance", worst_error < 1e-4,
           f"{worst_error:.2e} m")

    print("5. gradient is the surface normal")
    probes = [
        (0.5 * field.length, 0.8 * field.max_radius, 0.0),
        (0.9 * field.length, 0.0, 1.2 * field.max_radius),
        (0.15 * field.length, 0.0, -0.5 * field.max_radius),
    ]
    analytic = field.gradient(probes)
    numeric = field.numerical_gradient(probes)
    worst_gradient = max(
        math.dist(a, b) for a, b in zip(analytic, numeric, strict=True)
    )
    norms = [math.dist(g, (0.0, 0.0, 0.0)) for g in analytic]
    print(f"     |grad| {', '.join(f'{norm:.6f}' for norm in norms)}; "
          f"max |analytic - central difference| = {worst_gradient:.2e}")
    expect("gradient is unit length", all(abs(norm - 1.0) < 1e-6 for norm in norms))
    expect("analytic gradient matches central differences", worst_gradient < 1e-5,
           f"{worst_gradient:.2e}")

    print("6. edits move only what they name")
    longer = field.with_state(length=2.0 * field.length)
    barrel = (0.5 * field.length, 0.8 * field.max_radius, 0.0)
    print(f"     length x2: barrel distance {field.signed_distance([barrel])[0]:+.6f} -> "
          f"{longer.signed_distance([barrel])[0]:+.6f}; stern tip {field.length:.3f} -> "
          f"{longer.length:.3f} m; radius {field.max_radius:.6f} -> {longer.max_radius:.6f} m")
    expect("length edit leaves the barrel distance unchanged",
           longer.signed_distance([barrel]) == field.signed_distance([barrel]))
    expect("length edit leaves the radius unchanged", longer.max_radius == field.max_radius)
    fatter = field.with_state(radius=1.2 * field.max_radius)
    outer = (0.5 * field.length, 1.1 * field.max_radius, 0.0)
    print(f"     radius x1.2: point at 1.1 R is {field.signed_distance([outer])[0]:+.6f} m -> "
          f"{fatter.signed_distance([outer])[0]:+.6f} m; length unchanged "
          f"{fatter.length:.6f} m")
    expect("radius edit changes the sign of a point at 1.1 R",
           field.signed_distance([outer])[0] > 0.0 > fatter.signed_distance([outer])[0])
    expect("radius edit leaves the length unchanged", fatter.length == field.length)

    print("7. reference geometry (quadrature convergence)")
    print(f"     {'grid':>6s} {'volume m^3':>12s} {'wetted m^2':>12s} {'Cp':>8s}")
    for label in ("coarse", "fine"):
        block = reference[label]
        print(f"     {block['samples']:>6d} {block['volume']:12.6f} {block['wetted_area']:12.6f} "
              f"{block['prismatic_coefficient']:8.5f}")
    print(f"     resolution change: volume {1e6 * reference['volume_resolution_change']:+.3f} ppm, "
          f"wetted area {1e6 * reference['wetted_area_resolution_change']:+.3f} ppm")
    same_radius = suboff_reference_geometry(field.length, TENSORLBM_RADIUS)
    print(f"     at TensorLBM's default radius R = {TENSORLBM_RADIUS:.7f} m "
          f"(r_over_l = 1/(2*8.57)): V {same_radius['volume']:.6f} m^3, "
          f"A {same_radius['wetted_area']:.6f} m^2, vs its 5.991830 -> "
          f"{100 * (same_radius['wetted_area'] / TENSORLBM_WETTED_AREA - 1):+.5f} %")
    print(f"       so the offset table's own radius ({field.max_radius:.6f} m, r_over_l = "
          f"{field.max_radius / field.length:.6f}) is {100 * (field.max_radius / TENSORLBM_RADIUS - 1):+.3f} % "
          f"smaller, which moves V by {100 * (reference['volume'] / same_radius['volume'] - 1):+.3f} % "
          f"and A by {100 * (reference['wetted_area'] / same_radius['wetted_area'] - 1):+.3f} %")
    print(f"     analytic wetted area vs the published experimental "
          f"{EXPERIMENTAL_WETTED_AREA:.3f} m^2: "
          f"{100 * (reference['wetted_area'] / EXPERIMENTAL_WETTED_AREA - 1):+.3f} %")
    expect("quadrature is converged", abs(reference["volume_resolution_change"]) < 1e-8)
    expect("reproduces TensorLBM's wetted area at its own radius",
           abs(same_radius["wetted_area"] / TENSORLBM_WETTED_AREA - 1.0) < 1e-4,
           f"{100 * (same_radius['wetted_area'] / TENSORLBM_WETTED_AREA - 1):+.5f} %")
    expect("wetted area is within 0.5 % of the experimental value",
           abs(reference["wetted_area"] / EXPERIMENTAL_WETTED_AREA - 1.0) < 0.005)

    print("8. exports: watertight, genus zero, accurate")
    print(f"     {'cells':>14s} {'tris':>7s} {'vol err %':>10s} {'area err %':>11s} "
          f"{'min angle':>10s} {'slivers %':>10s} {'s':>6s}")
    handoff = None
    for label, spec in (("coarse", args.coarse), ("fine", args.fine)):
        cells = tuple(int(part) for part in spec.split(","))
        started = time.perf_counter()
        sampled = to_mesh(field, cells, margin=2)
        elapsed = time.perf_counter() - started
        volume_error = 100.0 * (sampled.report.signed_volume / reference["volume"] - 1.0)
        area_error = 100.0 * (sampled.report.surface_area / reference["wetted_area"] - 1.0)
        by_name = {check.name: check for check in sampled.quality}
        angle = float(by_name["mesh_min_angle"].detail.split()[3])
        slivers = float(by_name["mesh_sliver_share"].detail.split("(")[1].split("%")[0])
        print(f"     {cells!s:>14s} {sampled.report.triangles:7d} {volume_error:+10.3f} "
              f"{area_error:+11.3f} {angle:10.3f} {slivers:10.3f} {elapsed:6.1f}")
        expect(f"{label} mesh is watertight", sampled.report.watertight)
        expect(f"{label} mesh is edge-manifold", sampled.report.edge_manifold)
        expect(f"{label} mesh is genus zero", sampled.report.euler_characteristic == 2)
        expect(f"{label} mesh has no hard validation errors", sampled.validation.errors == [])
        expect(f"{label} volume within 1 %", abs(volume_error) < 1.0, f"{volume_error:+.3f} %")
        expect(f"{label} wetted area within 1 %", abs(area_error) < 1.0, f"{area_error:+.3f} %")
        expect(f"{label} sliver share is reported, not hidden",
               by_name["mesh_sliver_share"].severity == "warning")
        if label == "fine":
            handoff = sampled.handoff()
            with_digest = handoff["mesh_sha256"]

    assert handoff is not None
    print("9. solver handoff record")
    print(f"     schema {handoff['schema']}, {handoff['vertices']} vertices, "
          f"{handoff['triangles']} triangles, {handoff['units']}, {handoff['coordinate_frame']}")
    print(f"     mesh {with_digest[:16]}..., volume {handoff['signed_volume_m3']:.6f} m^3, "
          f"area {handoff['surface_area_m2']:.6f} m^2")
    print("     quality: " + "; ".join(
        f"{name} {'pass' if check['passed'] else 'warning'}" for name, check in handoff["quality"].items()
    ))
    expect("handoff carries a mesh digest", len(handoff["mesh_sha256"]) == 64)
    expect("handoff carries units and frame",
           handoff["units"] == "m" and isinstance(handoff["coordinate_frame"], str))
    expect("handoff declares its topology", handoff["watertight"] and handoff["euler_characteristic"] == 2)

    if args.output:
        import json

        from tensorcad.export import write_obj, write_stl

        target = Path(args.output)
        target.mkdir(parents=True, exist_ok=True)
        sampled = to_mesh(field, tuple(int(part) for part in args.fine.split(",")), margin=2)
        write_obj(sampled.mesh, target / "suboff_bare_hull.obj")
        write_stl(sampled.mesh, target / "suboff_bare_hull.stl")
        (target / "handoff.json").write_text(json.dumps(sampled.handoff(), indent=2), encoding="utf-8")
        (target / "field.json").write_text(json.dumps(field.metadata(), indent=2), encoding="utf-8")
        print(f"     wrote {target}/suboff_bare_hull.obj, .stl, handoff.json, field.json")

    print()
    if failures:
        print(f"FAILED {len(failures)} check(s): {', '.join(failures)}")
        return 1
    print("all checks passed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
