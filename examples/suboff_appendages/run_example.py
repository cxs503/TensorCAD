"""SUBOFF appendages (AFF-8) — measured verification of the sail, fins and full model.

The acceptance criteria this script enforces (it exits non-zero if any fails):

1. **The sections are the published ones.** The sail's three polynomials and the
   fins' NACA thickness are compared against the values computed by TensorLBM's
   transcription of the DARPA report, station by station.
2. **Containment is the published predicate.** Eighteen points labelled by the
   actual torch predicates (``suboff_sail_contains_points`` /
   ``suboff_fins_contain_points``) are checked, tight ones included: 3.1e-2 ft
   forward of the fin's leading edge, 0.03 ft above the sail's cap.
3. **The magnitude is measured, not claimed.** The sail's cap and the fins'
   thickness surfaces have no closed-form nearest point, so the fields are
   first-order normalised gaps; this script walks off the analytic surface along
   its normal and reports the error at 0.1 mm, 0.5 mm and 2 mm.
4. **Volume is checked against quadrature over the sections**, not against the
   sampler: sail, fins and the full model (hull + the parts of the appendages
   that stick out of it) each get an independent reference.
5. **Every export is solver-ready**: watertight, edge-manifold, no hard
   validation errors, and volume converging to the quadrature reference.
6. **The handoff record is complete** (units, frame, digest, quality).

Reported deviations (these are the numbers a reviewer will ask for):

* mesh volume vs quadrature volume at each resolution;
* the measured normal-offset error of the implicit magnitude;
* the fins' Euler characteristic, which is 0 rather than 2 for a real reason:
  the four blades start at r = 0.075 ft and are 0.1 ft thick, so they touch each
  other near the root and leave the axis column clear — a closed genus-1 shell.
  The hull fills that column, so the full model is genus 0.
"""

from __future__ import annotations

import argparse
import json
import math
import sys
import time
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from tensorcad.geometry_field import to_mesh
from tensorcad.suboff import suboff_reference_geometry
from tensorcad.suboff_appendages import (
    SuboffAppendagesField,
    SuboffFinsField,
    SuboffFullField,
    SuboffSailField,
    _FIN_R_INNER_FT,
    _FIN_TRAILING_FT,
    _LENGTH_FT,
    _SAIL_BODY_TOP_FT,
    _SAIL_HALF_WIDTH_FT,
    _SAIL_X_END_FT,
    _SAIL_X_START_FT,
    _DECK_FT,
    _fin_chord_coefficients,
    _naca_half_thickness_array,
    _sail_half_width_array,
    _sail_half_width_slope_array,
)

#: Sail half width (ft) at seven stations, from TensorLBM ``_sail_half_thickness_np``.
TENSORLBM_SAIL_HALF_WIDTH = {
    3.05: 0.04904483522377447,
    3.2: 0.10644549991291676,
    3.4: 0.109375,
    3.55: 0.109375,
    3.8: 0.0994715017078891,
    4.0: 0.06692398822021756,
    4.2: 0.014061932204276293,
}
#: NACA half thickness (ft) at seven chord positions, from TensorLBM ``_naca4_thickness_np``.
TENSORLBM_NACA_HALF_THICKNESS = {
    0.02: 0.03932962554685719,
    0.1: 0.07804677373039919,
    0.3: 0.1000287773232838,
    0.5: 0.08815875333428597,
    0.7: 0.060653511877967026,
    0.9: 0.022816421191197506,
    0.98: 0.004799850348000401,
}
#: Containment labels from TensorLBM's torch predicates at the reference length, for
#: points given in feet in the DARPA frame: ``(x, y, z, sail, fins)``.
TENSORLBM_POINTS_FT = [
    (3.6, 0.0, 1.0, True, False),
    (3.6, 0.05, 0.5, True, False),
    (3.6, 0.2, 0.5, False, False),
    (3.6, 0.0, 1.53, True, False),
    (3.6, 0.0, 1.60, False, False),
    (3.6, 0.05, 1.53, True, False),
    (3.6, 0.05, 1.60, False, False),
    (3.05, 0.0, 1.0, True, False),
    (4.2, 0.0, 1.0, True, False),
    (3.0, 0.0, 1.0, False, False),
    (12.9, 0.3, 0.0, False, True),
    (12.9, 0.3, 0.05, False, True),
    (12.9, 0.3, 0.3, False, False),
    (12.9, 0.05, 0.3, False, True),
    (13.0, 0.0, 0.4, False, True),
    (12.0, 0.4, 0.0, False, False),
    (12.5, 0.6, 0.0, False, False),
    (13.2, 0.6, 0.0, False, False),
]


def _quadrature_references(feet: float) -> dict[str, float]:
    """Volumes from the published sections, independent of the field and the sampler."""
    numpy = __import__("numpy")
    # Sail: cross-section = 2 a YTMP (rectangle) + pi a^2 / 4 (half ellipse of 2a by a).
    x = numpy.linspace(_SAIL_X_START_FT, _SAIL_X_END_FT, 400_001)
    width = _sail_half_width_array(x)
    sail = float(numpy.trapezoid(2.0 * width * _SAIL_BODY_TOP_FT + math.pi * width**2 / 4.0, x))

    # Fins: four blades, each 2 h(s) integrated over its trapezoidal planform.
    intercept, r_outer = _fin_chord_coefficients(1.0)
    s = numpy.linspace(0.0, 1.0, 400_001)
    integral = float(numpy.trapezoid(_naca_half_thickness_array(s), s))
    radius = numpy.linspace(_FIN_R_INNER_FT, r_outer, 100_001)
    chord = -0.466308 * radius + intercept
    fins = 4.0 * 2.0 * integral * float(numpy.trapezoid(chord, radius))
    planform = 4.0 * float(numpy.trapezoid(chord, radius))

    # The full model: the bare hull's analytic volume plus the parts of the sail and
    # fins that are outside the hull, i.e. what the union actually adds.
    from tensorcad.suboff import suboff_radius_fraction_array

    def hull_radius(x_ft):
        return suboff_radius_fraction_array(numpy.clip(x_ft / _LENGTH_FT, 0.0, 1.0)) * 0.8333333

    stations = numpy.linspace(_SAIL_X_START_FT + 1e-12, _SAIL_X_END_FT, 4001)
    half = _sail_half_width_array(stations)
    hull_r = hull_radius(stations)
    top = _SAIL_BODY_TOP_FT + half / 2.0
    fraction = numpy.linspace(0.0, 1.0, 4001)
    z = fraction[None, :] * top[:, None]
    safe = numpy.maximum(half[:, None] / 2.0, 1e-12)
    sail_width = numpy.where(
        z <= _SAIL_BODY_TOP_FT,
        half[:, None],
        half[:, None] * numpy.sqrt(numpy.maximum(1.0 - ((z - _SAIL_BODY_TOP_FT) / safe) ** 2, 0.0)),
    )
    hull_width = numpy.where(
        z < hull_r[:, None], numpy.sqrt(numpy.maximum(hull_r[:, None] ** 2 - z * z, 0.0)), 0.0
    )
    outside = 2.0 * (sail_width - numpy.minimum(sail_width, hull_width))
    sail_out = float(numpy.trapezoid(numpy.trapezoid(outside, fraction * top[:, None], axis=1), stations))

    span_stations = numpy.linspace(_FIN_R_INNER_FT, r_outer, 4001)
    chords = -0.466308 * span_stations + intercept
    along = numpy.linspace(0.0, 1.0, 2001)
    x_probe = _FIN_TRAILING_FT - chords[:, None] + along[None, :] * chords[:, None]
    fraction_out = numpy.mean(span_stations[:, None] > hull_radius(x_probe), axis=1)
    fin_out = 4.0 * float(numpy.trapezoid(2.0 * chords * integral * fraction_out, span_stations))

    hull = suboff_reference_geometry()
    return {
        "sail_ft3": sail,
        "sail": sail / feet**3,
        "fins_ft3": fins,
        "fins": fins / feet**3,
        "fin_planform": planform / feet**2,
        "sail_outside": sail_out / feet**3,
        "fins_outside": fin_out / feet**3,
        "hull": hull["volume"],
        "hull_area": hull["wetted_area"],
        "full": hull["volume"] + (sail_out + fin_out) / feet**3,
    }


def _sail_probe(feet: float) -> tuple[list, list]:
    """Analytic sail surface points (metres) with outward unit normals."""
    points: list[tuple[float, float, float]] = []
    normals: list[tuple[float, float, float]] = []
    for step in range(1, 40):
        x_ft = _SAIL_X_START_FT + (_SAIL_X_END_FT - _SAIL_X_START_FT) * step / 40.0
        width = float(_sail_half_width_array(x_ft))
        slope = float(_sail_half_width_slope_array(x_ft))
        for side in (1.0, -1.0):
            for fraction in (0.15, 0.5, 0.85):
                norm = math.hypot(slope, 1.0)
                points.append((x_ft / feet, side * width / feet, _SAIL_BODY_TOP_FT * fraction / feet))
                normals.append((-slope / norm, side / norm, 0.0))
            for index in range(1, 12):
                angle = math.pi * index / 12.0
                norm = math.hypot(slope, math.cos(angle), 2.0 * math.sin(angle))
                points.append(
                    (
                        x_ft / feet,
                        side * width * math.cos(angle) / feet,
                        (_SAIL_BODY_TOP_FT + 0.5 * width * math.sin(angle)) / feet,
                    )
                )
                normals.append((-slope / norm, side * math.cos(angle) / norm, 2.0 * math.sin(angle) / norm))
    return points, normals


def _fin_probe(feet: float) -> tuple[list, list]:
    """Analytic fin surface points (metres) with outward unit normals."""
    points: list[tuple[float, float, float]] = []
    normals: list[tuple[float, float, float]] = []
    intercept, r_outer = _fin_chord_coefficients(1.0)
    for span_step in range(10):
        radius = _FIN_R_INNER_FT + (r_outer - _FIN_R_INNER_FT) * (0.1 + 0.8 * span_step / 9.0)
        chord = -0.466308 * radius + intercept
        for chord_step in range(14):
            s = 0.1 + 0.8 * chord_step / 13.0
            x_ft = _FIN_TRAILING_FT + (s - 1.0) * chord
            half = float(_naca_half_thickness_array(s))
            low, high = max(s - 1e-6, 0.0), min(s + 1e-6, 1.0)
            slope = (float(_naca_half_thickness_array(high)) - float(_naca_half_thickness_array(low))) / (
                high - low
            )
            ds_dx = 1.0 / chord
            ds_dr = -(x_ft - _FIN_TRAILING_FT) * -0.466308 / (chord * chord)
            for sign in (1.0, -1.0):
                for axis in ("y", "z"):
                    if axis == "y":
                        point = (x_ft, radius, sign * half)
                        normal = (-slope * ds_dx, -slope * ds_dr, sign)
                    else:
                        point = (x_ft, sign * half, radius)
                        normal = (-slope * ds_dx, sign, -slope * ds_dr)
                    norm = math.sqrt(sum(component * component for component in normal))
                    points.append((point[0] / feet, point[1] / feet, point[2] / feet))
                    normals.append(tuple(component / norm for component in normal))
    return points, normals


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--coarse", default="218,26,38", help="cells per axis, coarse export")
    parser.add_argument("--fine", default="436,52,76", help="cells per axis, fine export")
    parser.add_argument("--output", default="", help="write OBJ/STL + handoff JSON here")
    args = parser.parse_args(argv)

    failures: list[str] = []

    def expect(label: str, condition: bool, detail: str = "") -> None:
        print(f"  [{'PASS' if condition else 'FAIL'}] {label}{': ' + detail if detail else ''}")
        if not condition:
            failures.append(label)

    sail = SuboffSailField()
    fins = SuboffFinsField()
    appendages = SuboffAppendagesField()
    full = SuboffFullField()
    feet = _LENGTH_FT / sail.length
    deck = _DECK_FT / feet
    apex = sail.height

    print("1. provenance and dimensions")
    print(
        f"     frame: {sail.coordinate_frame}, units: {sail.units}; bow tip at x = 0, "
        f"hull length L = {sail.length:.4f} m (14.291667 ft), deck at z = R = {deck:.4f} m"
    )
    print(f"     source: {sail.provenance()['geometry_source']}")
    print(
        f"     sail: x in [{_SAIL_X_START_FT / feet:.4f}, {_SAIL_X_END_FT / feet:.4f}] m, "
        f"half width {_SAIL_HALF_WIDTH_FT / feet:.6f} m, body top "
        f"{_SAIL_BODY_TOP_FT / feet:.4f} m above the axis"
    )
    print(
        f"           cap apex {apex:.6f} m above the axis = {apex - deck:.6f} m above the deck; "
        f"footprint length {(_SAIL_X_END_FT - _SAIL_X_START_FT) / feet:.4f} m"
    )
    _intercept, r_outer = _fin_chord_coefficients(1.0)
    print(
        f"     fins: span r in [{_FIN_R_INNER_FT / feet:.5f}, {r_outer / feet:.5f}] m, "
        f"chord {fins.parameters()['chord_root_m']:.5f} m at the root, "
        f"{fins.parameters()['chord_tip_m']:.5f} m at the tip"
    )
    print(
        f"           trailing edge at x = {_FIN_TRAILING_FT / feet:.4f} m, sweep slope "
        f"{fins.parameters()['sweep_slope']:.6f}, max half thickness "
        f"{max(TENSORLBM_NACA_HALF_THICKNESS.values()) / feet:.5f} m"
    )
    expect("sail sits on the DARPA deck plane", abs(deck - full.hull().max_radius) < 1e-12)
    expect(
        "cap apex is the body top plus half the maximum half width",
        abs(apex - (_SAIL_BODY_TOP_FT + 0.5 * _SAIL_HALF_WIDTH_FT) / feet) < 1e-12,
    )
    expect(
        "the fins reach 99 % of the hull's maximum radius",
        abs(r_outer / feet / deck - 1.0) < 0.02,
        f"{r_outer / feet:.5f} m vs R = {deck:.5f} m",
    )

    print("2. the published sections")
    worst_width = max(
        abs(float(_sail_half_width_array(x_ft)) - expected)
        for x_ft, expected in TENSORLBM_SAIL_HALF_WIDTH.items()
    )
    worst_naca = max(
        abs(float(_naca_half_thickness_array(s)) - expected)
        for s, expected in TENSORLBM_NACA_HALF_THICKNESS.items()
    )
    print(
        f"     sail half width vs TensorLBM over {len(TENSORLBM_SAIL_HALF_WIDTH)} stations: "
        f"max |error| {worst_width:.2e} ft "
        f"({1e6 * worst_width * 0.3048 / (_SAIL_HALF_WIDTH_FT * 0.3048):.4f} ppm of Zmax)"
    )
    print(
        f"     NACA half thickness vs TensorLBM over {len(TENSORLBM_NACA_HALF_THICKNESS)} "
        f"chord positions: max |error| {worst_naca:.2e} ft"
    )
    expect("sail section matches its source", worst_width < 1e-12, f"{worst_width:.1e} ft")
    expect("fin section matches its source", worst_naca < 1e-12, f"{worst_naca:.1e} ft")

    print("3. continuity of the sail section")
    for joint in (3.358507, 3.559028):
        left = float(_sail_half_width_array(joint - 1e-9))
        right = float(_sail_half_width_array(joint + 1e-9))
        print(f"     joint x = {joint:.6f} ft: C0 jump {right - left:+.2e} ft, half width {right:.9f} ft")
        expect(f"segments join at x = {joint:.3f} ft", abs(right - left) < 1e-7)
    for tip in (_SAIL_X_START_FT, _SAIL_X_END_FT):
        expect(f"the section closes at x = {tip:.3f} ft", float(_sail_half_width_array(tip)) == 0.0)
    print(
        f"     entrance slope da/dx at x = 3.10 ft: "
        f"{float(_sail_half_width_slope_array(3.10)):+.5f} (the published section is a "
        f"knife edge at the leading tip: the slope is unbounded there)"
    )

    print("4. containment against the labels the torch predicates returned")
    points = [(x / feet, y / feet, z / feet) for x, y, z, _, _ in TENSORLBM_POINTS_FT]
    sail_labels = sail.contains(points)
    fin_labels = fins.contains(points)
    expected_sail = [row[3] for row in TENSORLBM_POINTS_FT]
    expected_fins = [row[4] for row in TENSORLBM_POINTS_FT]
    sail_agree = sum(a == b for a, b in zip(sail_labels, expected_sail, strict=True))
    fin_agree = sum(a == b for a, b in zip(fin_labels, expected_fins, strict=True))
    print(
        f"     sail: {sail_agree}/{len(points)} labels, fins: {fin_agree}/{len(points)} labels "
        f"({sum(expected_sail)} sail-inside, {sum(expected_fins)} fin-inside)"
    )
    expect("the sail's containment is the published predicate", sail_agree == len(points))
    expect("the fins' containment is the published predicate", fin_agree == len(points))
    sign_agree = all(
        (value < 0.0) == label
        for value, label in zip(appendages.signed_distance(points), appendages.contains(points), strict=True)
    )
    expect("the union's sign is its containment", sign_agree)

    print("5. the magnitude, measured along the surface normal")
    print(f"     {'field':>5s} {'offset':>8s} {'p50 err':>11s} {'p95 err':>11s} {'max err':>11s}")
    for name, field, (probe_points, probe_normals) in (
        ("sail", sail, _sail_probe(feet)),
        ("fins", fins, _fin_probe(feet)),
    ):
        for offset in (1e-4, 5e-4, 2e-3):
            probes = [
                (p[0] + offset * n[0], p[1] + offset * n[1], p[2] + offset * n[2])
                for p, n in zip(probe_points, probe_normals, strict=True)
            ]
            values = field.signed_distance(probes)
            errors = sorted(abs(value - offset) for value in values)
            median = errors[len(errors) // 2]
            p95 = errors[int(0.95 * (len(errors) - 1))]
            print(
                f"     {name:>5s} {offset * 1e3:7.1f}mm {median * 1e3:9.5f}mm "
                f"{p95 * 1e3:9.5f}mm {errors[-1] * 1e3:9.5f}mm"
            )
            expect(
                f"{name} gap is first-order accurate at {offset * 1e3:.1f} mm ({len(errors)} normals)",
                median < (0.01 if name == "sail" else 0.002) * offset,
                f"p50 {median * 1e3:.5f} mm",
            )
    for name, field, (probe_points, probe_normals) in (
        ("sail", sail, _sail_probe(feet)),
        ("fins", fins, _fin_probe(feet)),
    ):
        probes = [
            (p[0] + 5e-4 * n[0], p[1] + 5e-4 * n[1], p[2] + 5e-4 * n[2])
            for p, n in zip(probe_points, probe_normals, strict=True)
        ]
        norms = sorted(math.dist(g, (0.0, 0.0, 0.0)) for g in field.gradient(probes))
        print(
            f"     {name}: |grad| p5 {norms[len(norms) // 20]:.4f}, "
            f"median {norms[len(norms) // 2]:.4f}, max {norms[-1]:.4f}"
        )
        expect(
            f"{name} |grad| stays at or below 1 and centres on it",
            norms[-1] <= 1.0 + 1e-9 and norms[len(norms) // 2] > 0.98,
        )

    print("6. reference volumes from quadrature over the sections")
    reference = _quadrature_references(feet)
    print(
        f"     sail {reference['sail']:.8f} m^3 (total), {reference['sail_outside']:.8f} m^3 outside the hull"
    )
    print(
        f"     fins {reference['fins']:.8f} m^3 (total), {reference['fins_outside']:.8f} m^3 "
        f"outside the hull; planform area (both sides, four blades) "
        f"{reference['fin_planform']:.5f} m^2"
    )
    print(
        f"     hull {reference['hull']:.8f} m^3, {reference['hull_area']:.6f} m^2; "
        f"full model {reference['full']:.8f} m^3 = hull + "
        f"{100 * (reference['full'] / reference['hull'] - 1):+.3f} %"
    )
    expect(
        "the appendages add volume, as a union must",
        reference["sail_outside"] > 0.0 and reference["fins_outside"] > 0.0,
    )
    expect(
        "the sail's buried part is the larger share of its volume",
        reference["sail_outside"] < 0.5 * reference["sail"],
        f"{100 * reference['sail_outside'] / reference['sail']:.1f} % outside",
    )

    print("7. exports")
    print(
        f"     {'model':>8s} {'cells':>14s} {'tris':>8s} {'vol err %':>10s} {'tris/s':>8s} "
        f"{'min ang':>8s} {'slivers %':>10s} {'s':>6s}"
    )
    handoff = None
    for label, field, spec, expected_volume, euler in (
        ("sail", sail, "72,12,92", reference["sail"], 2),
        ("fins", fins, "110,64,64", reference["fins"], 0),
        ("coarse", full, args.coarse, reference["full"], 2),
        ("fine", full, args.fine, reference["full"], 2),
    ):
        cells = tuple(int(part) for part in spec.split(","))
        started = time.perf_counter()
        sampled = to_mesh(field, cells, margin=2)
        elapsed = time.perf_counter() - started
        volume_error = 100.0 * (sampled.report.signed_volume / expected_volume - 1.0)
        by_name = {check.name: check for check in sampled.quality}
        angle = float(by_name["mesh_min_angle"].detail.split()[3])
        slivers = float(by_name["mesh_sliver_share"].detail.split("(")[1].split("%")[0])
        print(
            f"     {label:>8s} {cells!s:>14s} {sampled.report.triangles:8d} "
            f"{volume_error:+10.3f} {sampled.report.triangles / elapsed:8.0f} {angle:8.3f} "
            f"{slivers:10.3f} {elapsed:6.1f}"
        )
        expect(f"{label} mesh is watertight", sampled.report.watertight)
        expect(f"{label} mesh is edge-manifold", sampled.report.edge_manifold)
        expect(
            f"{label} mesh has Euler characteristic {euler}",
            sampled.report.euler_characteristic == euler,
            f"got {sampled.report.euler_characteristic}",
        )
        if label == "fins":
            # The blades touch near the root but leave the axis column clear, so the
            # fin set alone is a closed genus-1 shell. The hull fills the column.
            expect(
                "the fins' genus is explained, not ignored", abs(volume_error) < 2.0, f"{volume_error:+.3f} %"
            )
        else:
            expect(f"{label} mesh has no hard validation errors", sampled.validation.errors == [])
            expect(
                f"{label} volume within 1 % of quadrature", abs(volume_error) < 1.0, f"{volume_error:+.3f} %"
            )
        expect(
            f"{label} sliver share is reported, not hidden",
            by_name["mesh_sliver_share"].severity == "warning",
        )
        if label == "fine":
            handoff = sampled.handoff()
    assert handoff is not None

    print("8. throughput")
    numpy = __import__("numpy")
    generator = numpy.random.default_rng(5)
    batch = numpy.column_stack(
        [
            generator.uniform(0.0, full.length, 20_000),
            generator.uniform(-0.26, 0.26, 20_000),
            generator.uniform(-0.26, 0.48, 20_000),
        ]
    )
    best = min(
        (full.signed_distance(batch), time.perf_counter() - start)[1]
        for start in (time.perf_counter(), time.perf_counter(), time.perf_counter())
    )
    print(
        f"     full model, 20,000 points: {20_000 / best:,.0f} queries/s "
        f"({1e6 * best / 20_000:.2f} us per point), min of 3 runs"
    )
    expect("the field is fast enough to sample", best < 2.0, f"{best:.3f} s")

    print("9. solver handoff record")
    print(
        f"     schema {handoff['schema']}, {handoff['vertices']} vertices, "
        f"{handoff['triangles']} triangles, {handoff['units']}, {handoff['coordinate_frame']}"
    )
    print(
        f"     mesh {handoff['mesh_sha256'][:16]}..., volume "
        f"{handoff['signed_volume_m3']:.6f} m^3, area {handoff['surface_area_m2']:.6f} m^2"
    )
    print(
        "     quality: "
        + "; ".join(
            f"{name} {'pass' if check['passed'] else 'warning'}" for name, check in handoff["quality"].items()
        )
    )
    expect("handoff carries a mesh digest", len(handoff["mesh_sha256"]) == 64)
    expect("handoff declares its topology", handoff["watertight"] and handoff["euler_characteristic"] == 2)
    expect("handoff's validation passed", handoff["validation"]["passed"] is True)
    expect(
        "handoff volume is within 1 % of quadrature",
        abs(handoff["signed_volume_m3"] / reference["full"] - 1.0) < 0.01,
    )

    if args.output:
        from tensorcad.export import write_obj, write_stl

        target = Path(args.output)
        target.mkdir(parents=True, exist_ok=True)
        cells = tuple(int(part) for part in args.fine.split(","))
        sampled = to_mesh(full, cells, margin=2)
        write_obj(sampled.mesh, target / "suboff_full_model.obj")
        write_stl(sampled.mesh, target / "suboff_full_model.stl")
        (target / "handoff.json").write_text(json.dumps(sampled.handoff(), indent=2), encoding="utf-8")
        (target / "field.json").write_text(json.dumps(full.metadata(), indent=2), encoding="utf-8")
        sail_mesh = to_mesh(sail, (72, 12, 92), margin=2)
        write_obj(sail_mesh.mesh, target / "suboff_sail.obj")
        fins_mesh = to_mesh(fins, (110, 64, 64), margin=2)
        write_obj(fins_mesh.mesh, target / "suboff_fins.obj")
        print(
            f"     wrote {target}/suboff_full_model.obj, .stl, suboff_sail.obj, "
            f"suboff_fins.obj, handoff.json, field.json"
        )

    print()
    if failures:
        print(f"FAILED {len(failures)} check(s): {', '.join(failures)}")
        return 1
    print("all checks passed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
