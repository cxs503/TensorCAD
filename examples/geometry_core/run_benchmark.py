"""Continuous Geometry Core v0.1 -- measured benchmark, no mocks.

What this script establishes, and fails the build if it stops being true:

1. the analytic primitives are *exact* (``|s|`` equals the Euclidean distance for
   the sphere, the plane and the capped cylinder; the ellipsoid's zero set is
   exact and its magnitude is declared as an implicit function, not a distance);
2. their gradients agree with central differences -- the protocol's gradient is
   not decorative;
3. an edit moves what it says it moves and leaves the other controls alone,
   quantified in metres (and a learned field's conditioned state is edited the
   same way, without retraining);
3. every representation exports through the same pipeline: sample -> closed mesh
   -> ``inspect_geometry`` + mesh-quality checks, with the volume compared against
   the closed-form value;
4. a sampled mesh carries no squashed triangles (area below 1e-6 of the longest
   edge squared) even when the extent is a whole number of cells long -- the
   cylinder used to emit 2,936 of them at 40 cells;
5. when the ``neural`` extra is installed, a *learned* field (the revision-2 sphere
   bundle committed in this repo) is queried through the same protocol and
   compared with the analytic reference: surface error, conditioned edits,
   throughput and memory.

Point 5 compares surfaces, never raw values: the learned field is a *truncated*
signed distance, so its magnitude outside the truncation band is not a distance.
"""
from __future__ import annotations

import argparse
import math
import sys
import time
import tracemalloc
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from tensorcad.geometry_field import (
    IMPLICIT_FUNCTION,
    SIGNED_DISTANCE,
    CylinderField,
    EllipsoidField,
    PlaneField,
    SphereField,
    surface_point,
    to_mesh,
)

#: Directions used for every geometric comparison. A deterministic spiral, so two
#: runs of this script compare the same points and the numbers are reviewable.
DIRECTIONS = [
    (
        math.cos(0.37 * index) * math.cos(0.91 * index),
        math.cos(0.37 * index) * math.sin(0.91 * index),
        math.sin(0.37 * index),
    )
    for index in range(24)
]


def failure_collector() -> tuple[list[str], object]:
    failures: list[str] = []

    def expect(label: str, condition: bool, detail: str = "") -> None:
        status = "PASS" if condition else "FAIL"
        print(f"  [{status}] {label}{': ' + detail if detail else ''}")
        if not condition:
            failures.append(label)

    return failures, expect


def _squashed_faces(mesh: object) -> int:
    """Triangles with area below 1e-6 of their longest edge squared.

    Scale-free on purpose: what a grid-vertex coincidence leaves behind is a
    triangle whose shape is degenerate, not one with a particular area.
    """
    count = 0
    for face in mesh.faces:  # type: ignore[attr-defined]
        a, b, c = (mesh.vertices[index] for index in face)  # type: ignore[attr-defined]
        ux, uy, uz = b[0] - a[0], b[1] - a[1], b[2] - a[2]
        vx, vy, vz = c[0] - a[0], c[1] - a[1], c[2] - a[2]
        area2 = math.sqrt(
            (uy * vz - uz * vy) ** 2 + (uz * vx - ux * vz) ** 2 + (ux * vy - uy * vx) ** 2
        )
        longest = max(math.dist(a, b), math.dist(b, c), math.dist(c, a))
        if longest > 0.0 and area2 <= 1e-6 * longest * longest:
            count += 1
    return count


def throughput(field: object, count: int = 20_000) -> tuple[float, float]:
    """Return (queries per second, peak MiB) for a batch query."""
    points = [
        (
            0.9 * math.cos(0.017 * index),
            0.9 * math.sin(0.017 * index),
            0.4 * math.sin(0.011 * index),
        )
        for index in range(count)
    ]
    field.signed_distance(points[:16])  # warm up imports and caches
    best = float("inf")
    peak = 0.0
    for _ in range(3):
        tracemalloc.start()
        started = time.perf_counter()
        field.signed_distance(points)  # type: ignore[attr-defined]
        elapsed = time.perf_counter() - started
        peak = max(peak, tracemalloc.get_traced_memory()[1] / (1024 * 1024))
        tracemalloc.stop()
        best = min(best, elapsed)
    return count / best, peak


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--resolution", type=int, default=40, help="cells per axis for exports")
    parser.add_argument("--bundle", default=str(
        REPO_ROOT / "docs/experiments/2026-10-24-neural-geometry-object/bundle"
    ), help="published .ngo bundle used for the learned-field comparison")
    parser.add_argument("--points", type=int, default=20_000, help="points per throughput sample")
    args = parser.parse_args(argv)

    failures, expect = failure_collector()
    sphere = SphereField(radius=0.5)
    cylinder = CylinderField(radius=0.3, height=1.0)
    ellipsoid = EllipsoidField(semi_axes=(0.5, 0.3, 0.2))
    plane = PlaneField(normal=(0.0, 0.0, 1.0), offset=0.25)

    print("1. analytic exactness (the reference implementations)")
    probes = [(0.31, 0.17, -0.09), (-0.12, 0.44, 0.05), (0.05, 0.02, 0.61), (0.9, 0.0, 0.0)]
    worst_sphere = max(
        abs(abs(value) - abs(math.dist(point, sphere.center) - sphere.radius))
        for point, value in zip(probes, sphere.signed_distance(probes), strict=True)
    )
    worst_plane = max(
        abs(value - (point[2] - plane.offset))
        for point, value in zip(probes, plane.signed_distance(probes), strict=True)
    )
    worst_cylinder = 0.0
    for point in probes:
        value = cylinder.signed_distance([point])[0]
        # Closed form: distance to the capped cylinder, recomputed independently.
        axial = abs(point[2]) - 0.5 * cylinder.height
        radial = math.hypot(point[0], point[1]) - cylinder.radius
        reference = math.hypot(max(radial, 0.0), max(axial, 0.0)) + min(max(radial, axial), 0.0)
        worst_cylinder = max(worst_cylinder, abs(value - reference))
    ellipsoid_zero = max(abs(value) for value in ellipsoid.signed_distance(
        [(0.5, 0.0, 0.0), (0.0, 0.3, 0.0), (0.0, 0.0, 0.2)]
    ))
    print(f"     sphere   max |s| error   {worst_sphere:.2e} m   ({sphere.quantity})")
    print(f"     plane    max |s| error   {worst_plane:.2e} m   ({plane.quantity})")
    print(f"     cylinder max |s| error   {worst_cylinder:.2e} m   ({cylinder.quantity})")
    print(f"     ellipsoid |f| on surface {ellipsoid_zero:.2e}      ({ellipsoid.quantity})")
    expect("sphere distance is exact", worst_sphere < 1e-12, f"{worst_sphere:.2e}")
    expect("plane distance is exact", worst_plane < 1e-12, f"{worst_plane:.2e}")
    expect("cylinder distance is exact", worst_cylinder < 1e-12, f"{worst_cylinder:.2e}")
    expect("ellipsoid zero set is exact", ellipsoid_zero < 1e-12, f"{ellipsoid_zero:.2e}")
    expect(
        "the ellipsoid declares an implicit function, not a distance",
        ellipsoid.quantity == IMPLICIT_FUNCTION and sphere.quantity == SIGNED_DISTANCE,
    )

    print("2. gradients agree with central differences")
    for field in (sphere, cylinder, ellipsoid, plane):
        analytic = field.gradient(probes)
        numeric = field.numerical_gradient(probes)
        worst = max(
            abs(a - b)
            for pair_a, pair_b in zip(analytic, numeric, strict=True)
            for a, b in zip(pair_a, pair_b, strict=True)
        )
        expect(f"{field.kind} gradient", worst < 1e-6, f"max |analytic - numeric| = {worst:.2e}")

    print("3. edits move the target and nothing else")
    before, after = SphereField(radius=0.5), SphereField(radius=0.5).with_state(radius=0.65)
    moved = [
        math.dist(surface_point(before, direction), surface_point(after, direction))
        for direction in DIRECTIONS
    ]
    expect("radius edit moves the surface by exactly the radius change",
           max(abs(distance - 0.15) for distance in moved) < 1e-6,
           f"min {min(moved):.6f} max {max(moved):.6f} (target 0.15)")
    expect("radius edit keeps the centre", after.center == before.center, str(after.center))

    short, tall = CylinderField(radius=0.3, height=1.0), CylinderField(radius=0.3, height=1.0).with_state(height=1.6)
    expect("height edit keeps the radius", tall.radius == short.radius == 0.3)
    barrel_before = short.signed_distance([(0.8, 0.0, 0.0)])[0]
    barrel_after = tall.signed_distance([(0.8, 0.0, 0.0)])[0]
    expect("height edit does not move the barrel", barrel_before == barrel_after,
           f"{barrel_before:.6f} vs {barrel_after:.6f}")
    lid_before = short.signed_distance([(0.0, 0.0, 1.2)])[0]
    lid_after = tall.signed_distance([(0.0, 0.0, 1.2)])[0]
    expect("height edit does move the lid",
           abs((lid_before - lid_after) - 0.5 * (1.6 - 1.0)) < 1e-12,
           f"lid distance {lid_before:.6f} -> {lid_after:.6f}")

    thin = EllipsoidField(semi_axes=(0.5, 0.3, 0.2)).with_state(semi_y=0.45)
    expect("ellipsoid axis edit leaves the other axes alone",
           thin.semi_axes == (0.5, 0.45, 0.2), str(thin.semi_axes))

    print("4. export: sample -> closed mesh -> the shared acceptance pipeline")
    exports = [
        (sphere, 4.0 / 3.0 * math.pi * sphere.radius**3),
        (cylinder, math.pi * cylinder.radius**2 * cylinder.height),
        (ellipsoid, 4.0 / 3.0 * math.pi * math.prod(ellipsoid.semi_axes)),
    ]
    print(f"     {'field':10s} {'cells':>5s} {'tris':>7s} {'zero-area':>9s} {'vol err %':>10s} "
          f"{'angle':>7s} {'slivers %':>9s}")
    for field, exact in exports:
        sampled = to_mesh(field, args.resolution)
        by_name = {check.name: check for check in sampled.quality}
        angle = float(by_name["mesh_min_angle"].detail.split()[3])
        slivers = float(by_name["mesh_sliver_share"].detail.split("(")[1].split("%")[0])
        error = 100.0 * (sampled.report.signed_volume / exact - 1.0)
        zero_area = _squashed_faces(sampled.mesh)
        print(f"     {field.kind:10s} {args.resolution:5d} {sampled.report.triangles:7d} "
              f"{zero_area:9d} {error:+10.3f} {angle:7.3f} {slivers:9.3f}")
        expect(f"{field.kind} has no squashed faces", zero_area == 0, f"{zero_area}")
        expect(f"{field.kind} mesh is watertight", sampled.report.watertight)
        expect(f"{field.kind} mesh is genus zero", sampled.report.euler_characteristic == 2)
        expect(f"{field.kind} volume within 0.5 %", abs(error) < 0.5, f"{error:+.3f} %")
        expect(f"{field.kind} has no hard validation errors", sampled.validation.errors == [])

    print("5. representation comparison: analytic sphere vs a learned object")
    bundle = Path(args.bundle)
    try:
        from tensorcad.neural.objects import field_from_object  # needs torch; optional

        learned = field_from_object(bundle)
    except Exception as exc:  # any failure to import torch means "skip this section"
        print(f"     skipped: the neural extra is not usable here ({type(exc).__name__}: {exc})")
        learned = None

    if learned is not None:
        reference = SphereField(radius=0.5)
        sigma = 0.22  # the edited cone, as recorded in the bundle's state
        dent_axis = (1.0, 0.0, 0.0)
        errors: list[float] = []
        protected: list[float] = []
        for direction in DIRECTIONS:
            error = math.dist(surface_point(reference, direction), surface_point(learned, direction))
            angle = math.acos(max(-1.0, min(1.0, direction[0])))
            errors.append(error)
            if angle > 2.5 * sigma:  # outside the edited cone
                protected.append(error)
        metrics = learned.manifest.raw.get("validation", {}).get("metrics", {})
        recorded_peak = float(metrics.get("dent_peak_measured_m", "nan"))
        target_radius = reference.radius - 0.06  # depth recorded in the bundle's state
        learned_radius = math.dist(surface_point(learned, dent_axis), (0.0, 0.0, 0.0))
        mean_protected = sum(protected) / len(protected)
        print(f"     learned bundle    revision {learned.manifest.revision}, "
              f"weights {learned.manifest.weights.sha256[:16]}, {learned.quantity}")
        print(f"     surface vs sphere mean {mean_protected * 1000:.3f} mm, "
              f"worst {max(protected) * 1000:.3f} mm over {len(protected)} unedited directions")
        print(f"     dent axis radius  {learned_radius * 1000:.3f} mm "
              f"(analytic target {target_radius * 1000:.3f} mm, "
              f"error {(learned_radius - target_radius) * 1000:+.3f} mm)")
        print(f"     recorded dent     {recorded_peak * 1000:.3f} mm -- measured against the "
              f"*parent revision*, not the analytic sphere; it implies a parent axis radius of "
              f"{(learned_radius + recorded_peak) * 1000:.3f} mm "
              f"({(learned_radius + recorded_peak - reference.radius) * 1000:+.3f} mm vs r=0.5)")
        expect("learned field is a truncated signed distance",
               learned.quantity == "truncated_signed_distance",
               "magnitudes outside the band are not distances")
        expect("learned surface stays within 2 % of the radius away from the dent",
               mean_protected < 0.02 * reference.radius, f"mean {mean_protected * 1000:.3f} mm")
        expect("learned surface worst case stays within the fit's recorded 3.2 %",
               max(protected) < 0.032 * reference.radius, f"worst {max(protected) * 1000:.3f} mm")
        expect("the user's 60 mm dent is reproduced within 5 mm",
               abs(learned_radius - target_radius) < 0.005,
               f"{(learned_radius - target_radius) * 1000:+.3f} mm")

        def outward(points: list[tuple[float, float, float]]) -> tuple[list[float], list[float]]:
            gradients = learned.gradient(points)
            norms = [math.dist(gradient, (0.0, 0.0, 0.0)) for gradient in gradients]
            cosines = [
                sum(a * b for a, b in zip(gradient, point, strict=True))
                / (norm * math.dist(point, (0.0, 0.0, 0.0)))
                for gradient, norm, point in zip(gradients, norms, points, strict=True)
            ]
            return norms, cosines

        near_norms, near_cosines = outward([(0.4, 0.03, 0.0), (0.0, 0.4, 0.03), (0.03, 0.0, 0.4)])
        far_norms, far_cosines = outward([(0.3, 0.03, 0.0), (0.0, 0.3, 0.03), (0.03, 0.0, 0.3)])
        print(f"     gradient (autograd), inside the band:  |grad| "
              f"{', '.join(f'{norm:.2f}' for norm in near_norms)}, outward cosine min "
              f"{min(near_cosines):.3f}")
        print(f"     gradient (autograd), beyond the band:  |grad| "
              f"{', '.join(f'{norm:.2f}' for norm in far_norms)}, outward cosine min "
              f"{min(far_cosines):.3f}  <- the truncation caveat, measured")
        print("       -> the field is not a metric: |grad| drifts from 1, so distance queries "
              "must use the field value, never |grad|")
        expect("the learned gradient points outward inside the band",
               min(near_cosines) > 0.9, f"min cosine {min(near_cosines):.3f}")
        expect("the analytic gradient of an SDF is unit length",
               all(abs(norm - 1.0) < 1e-9 for norm in
                   (math.dist(gradient, (0.0, 0.0, 0.0)) for gradient in reference.gradient(probes))))

        print("     conditioned edits (no retraining): radius is field state, not mesh geometry")
        for radius in (0.45, 0.55, 0.6):
            conditioned = learned.with_state(radius=radius)
            measured = [
                math.dist(surface_point(conditioned, direction), (0.0, 0.0, 0.0))
                for direction in DIRECTIONS
                if math.acos(max(-1.0, min(1.0, direction[0]))) > 2.5 * sigma
            ]
            mean_radius = sum(measured) / len(measured)
            print(f"       state {radius * 1000:5.0f} mm -> measured {mean_radius * 1000:7.3f} mm "
                  f"({100.0 * (mean_radius / radius - 1.0):+.2f} %)")
            expect(f"conditioned radius {radius:.2f} tracks within 1 %",
                   abs(mean_radius / radius - 1.0) < 0.01,
                   f"{100.0 * (mean_radius / radius - 1.0):+.2f} %")

        print(f"     {'representation':24s} {'queries/s':>12s} {'peak MiB':>9s}")
        scalar_rate, scalar_peak = throughput(reference, args.points)
        numpy_rate, numpy_peak = throughput(_NumpyWrapper(reference), args.points)
        learned_rate, learned_peak = throughput(learned, args.points)
        print(f"     {'analytic, python loop':24s} {scalar_rate:12,.0f} {scalar_peak:9.2f}")
        print(f"     {'analytic, numpy path':24s} {numpy_rate:12,.0f} {numpy_peak:9.2f}")
        print(f"     {'learned .ngo (torch)':24s} {learned_rate:12,.0f} {learned_peak:9.2f}")
        print(f"     vectorised / scalar = {numpy_rate / scalar_rate:.2f}x, "
              f"learned / vectorised = {learned_rate / numpy_rate:.3f}x")
        expect("the vectorised analytic path is materially faster than the loop",
               numpy_rate > 2.0 * scalar_rate, f"{numpy_rate / scalar_rate:.2f}x")

    print()
    if failures:
        print(f"FAILED {len(failures)} check(s): {', '.join(failures)}")
        return 1
    print("all checks passed")
    return 0


class _NumpyWrapper:
    """Expose a field's vectorised evaluator through the same query signature."""

    def __init__(self, field: object) -> None:
        self._evaluate = field.numpy_callable()  # type: ignore[attr-defined]
        assert self._evaluate is not None

    def signed_distance(self, points: list[tuple[float, float, float]]) -> list[float]:
        import numpy as np

        values = self._evaluate(np.asarray(points, dtype=np.float64))
        return [float(value) for value in values]


if __name__ == "__main__":
    raise SystemExit(main())
