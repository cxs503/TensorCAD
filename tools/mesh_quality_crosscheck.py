"""Differential cross-check of :mod:`tensorcad.geomquality.quality`.

What this tool asserts, in the order it prints:

1. **Closed form.**  Every measure on canonical elements equals its ideal/analytic
   value: equilateral (all normalised measures 1, condition number ``1/sqrt(3)``),
   right isosceles, 3-4-5, degenerate.  Tolerance ``1e-12``.
2. **Independent route.**  :mod:`tools.quality_oracle` recomputes every measure along a
   different arithmetic path (Heron instead of cross product, LAPACK SVD instead of the
   invariant identity, ``2r/R`` instead of ``8A^2/(s abc)``).  Agreement on every mesh,
   every measure: ``< 1e-6``.
3. **Exact arithmetic.**  On the rational meshes the oracle runs in
   :class:`fractions.Fraction` up to one square root, so the float library's rounding is
   bounded, not just self-consistent.  Tolerance ``1e-12``.
4. **Upstream.**  The real ``vtk.vtkMeshQuality`` (what ``meshmagick.Mesh.print_quality``
   calls) and ``gmsh.model.mesh.getElementQualities`` are run *if importable*; this host
   has neither, and the report records that with the reason instead of pretending.
5. **Dialect audit.**  The four in-house "aspect ratio" conventions are measured against
   the canonical definitions on the same meshes, so their disagreement is a number.
6. **Gate.**  The default CI gate plus a threshold sweep on every mesh.

Artifacts land in ``artifacts/mesh_quality_lib/``: ``crosscheck.json`` and ``gate.json``.
Exit code is 1 if any tolerance is violated, so the tool doubles as the gate itself.

Usage::

    python tools/mesh_quality_crosscheck.py [--json-dir artifacts/mesh_quality_lib]
"""
from __future__ import annotations

import argparse
import json
import math
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:  # allow `python tools/...` from a source checkout
    sys.path.insert(0, str(ROOT))
TOOLS = Path(__file__).resolve().parent
if str(TOOLS) not in sys.path:
    sys.path.insert(0, str(TOOLS))

from tensorcad import gate as gate_module  # noqa: E402
from tensorcad import validation as validation_module  # noqa: E402
from tensorcad.geometry import Mesh, MeshBuilder  # noqa: E402
from tensorcad.geomquality.quality import (  # noqa: E402
    DEFAULT_CRITERIA,
    IDEAL_VALUES,
    METRIC_COMPARISON,
    METRIC_DEFINITIONS,
    mesh_quality,
    quality_gate,
    triangle_quality,
)
from tensorcad.mesh_conditioning import measure_mesh as conditioning_measure  # noqa: E402
from tensorcad.parametric_hull import mesh_quality_report  # noqa: E402

import quality_oracle  # noqa: E402

#: Max |library - independent route| accepted for a measure on a mesh.
ROUTE_TOLERANCE = 1e-6
#: Max |library - exact arithmetic| accepted on a rational mesh.  Angles get their own
#: (looser, still physically meaningless) bound: every angle here goes through
#: ``acos``, whose condition number is ``1/sqrt(1-x^2)``, so a sliver with a 0.01 degree
#: angle amplifies a 1e-16 relative error in ``cos`` to ~1e-12 degrees.  That is a
#: property of the inverse cosine, not of the mesh.
EXACT_TOLERANCE = 1e-12
ANGLE_TOLERANCE = 1e-9


def _tolerance(metric: str, default: float) -> float:
    """Per-measure tolerance: the angle bound is the ``acos`` bound, not an excuse."""
    return ANGLE_TOLERANCE if metric.endswith("_deg") else default
#: Max deviation from a closed-form ideal value on a canonical element.
ANALYTIC_TOLERANCE = 1e-12


def _canonical_triangles() -> dict[str, tuple[tuple[float, float, float], ...]]:
    return {
        "equilateral": ((0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (0.5, math.sqrt(3) / 2.0, 0.0)),
        "right_isosceles": ((0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (0.0, 1.0, 0.0)),
        "three_four_five": ((0.0, 0.0, 0.0), (3.0, 0.0, 0.0), (3.0, 4.0, 0.0)),
        "sliver_1e-4": ((0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (0.5, 1e-4, 0.0)),
        "degenerate": ((0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (2.0, 0.0, 0.0)),
    }


def _box_mesh() -> Mesh:
    """Unit cube, integer coordinates, 12 triangles -- the exact-arithmetic mesh."""
    builder = MeshBuilder()
    corners = [
        (0.0, 0.0, 0.0),
        (1.0, 0.0, 0.0),
        (1.0, 1.0, 0.0),
        (0.0, 1.0, 0.0),
        (0.0, 0.0, 1.0),
        (1.0, 0.0, 1.0),
        (1.0, 1.0, 1.0),
        (0.0, 1.0, 1.0),
    ]
    for index in corners:
        builder.add_vertex(index)
    quads = [
        (0, 3, 2, 1),
        (4, 5, 6, 7),
        (0, 1, 5, 4),
        (1, 2, 6, 5),
        (2, 3, 7, 6),
        (3, 0, 4, 7),
    ]
    for a, b, c, d in quads:
        builder.add_triangle(a, b, c)
        builder.add_triangle(a, c, d)
    return builder.build()


def _sliver_patch() -> Mesh:
    """A rational patch with a good element, a mild sliver and a full sliver.

    Coordinates are integers or halves, so the exact route can run on it; two of the
    three triangles must fail the default gate.
    """
    builder = MeshBuilder()
    points = [
        (0.0, 0.0, 0.0),
        (4.0, 0.0, 0.0),
        (4.0, 3.0, 0.0),  # 3-4-5 right triangle: fine
        (8.0, 0.0, 0.0),
        (8.0, 0.5, 0.0),  # base 4, height 0.5 : mild sliver
        (12.0, 0.0, 0.0),
        (12.0, 0.01, 0.0),  # base 4, height 0.01 : hard sliver
        (16.0, 0.0, 0.0),
    ]
    for point in points:
        builder.add_vertex(point)
    builder.add_triangle(0, 1, 2)
    builder.add_triangle(3, 4, 5)
    builder.add_triangle(5, 6, 7)
    return builder.build()


def _hull_meshes() -> dict[str, Mesh]:
    """Three real hull tessellations from this repository's own generators."""
    from tensorcad.geometry_field import to_mesh
    from tensorcad.series60 import generate_series60
    from tensorcad.suboff_appendages import SuboffFullField
    from tensorcad.wigley import generate_wigley

    meshes: dict[str, Mesh] = {}
    meshes["wigley"] = generate_wigley().mesh
    meshes["series60"] = generate_series60().mesh
    started = time.time()
    sampled = to_mesh(SuboffFullField(), (48, 8, 12), margin=2)
    meshes["suboff_coarse"] = sampled.mesh
    from tensorcad.mesh_conditioning import taubin_smooth

    smoothed, _drift = taubin_smooth(sampled.mesh, iterations=30)
    meshes["suboff_taubin30"] = smoothed
    return meshes


def _probe_upstream() -> dict[str, dict[str, object]]:
    """Which real oracle could be run on this host, and why not when it cannot."""
    probe: dict[str, dict[str, object]] = {}
    for name, module in (("vtk", "vtkmodules.all"), ("gmsh", "gmsh"), ("meshmagick", "meshmagick")):
        target = module if name == "vtk" else name
        try:
            __import__(target)
            probe[name] = {"available": True, "reason": ""}
        except Exception as error:  # noqa: BLE001 - the reason is the payload
            probe[name] = {"available": False, "reason": f"{type(error).__name__}: {error}"}
    probe["meshmagick"]["note"] = (
        "meshmagick.Mesh.print_quality is the meshmagick quality entry point; it needs "
        "VTK (meshmagick/mesh.py:459-540 wraps vtk.vtkMeshQuality), so it is only "
        "runnable when vtk is available."
    )
    probe["gmsh"]["note"] = (
        "the Gmsh measures compared here (gamma, minSICN family, min/max edge) are read "
        "back through gmsh.model.mesh.getElementQualities on a .msh file."
    )
    return probe


def _vtk_compare(meshes: dict[str, Mesh]) -> dict[str, object]:
    """Run vtkMeshQuality on the meshes when VTK exists and diff it against the library."""
    try:
        import vtkmodules.all as vtk  # type: ignore[import-not-found]
    except Exception as error:  # noqa: BLE001
        return {"ran": False, "reason": f"{type(error).__name__}: {error}"}
    measures = {
        "edge_ratio": ("SetTriangleQualityMeasureToEdgeRatio", "max"),
        "aspect_ratio": ("SetTriangleQualityMeasureToAspectRatio", "max"),
        "aspect_frobenius": ("SetTriangleQualityMeasureToAspectFrobenius", "max"),
        "min_angle_deg": ("SetTriangleQualityMeasureToMinAngle", "min"),
    }
    out: dict[str, object] = {"ran": True, "per_mesh": {}}
    for name, mesh in meshes.items():
        points = vtk.vtkPoints()
        for vertex in mesh.vertices:
            points.InsertNextPoint(*vertex)
        cells = vtk.vtkCellArray()
        for face in mesh.faces:
            triangle = vtk.vtkTriangle()
            for position, index in enumerate(face):
                triangle.GetPointIds().SetId(position, index)
            cells.InsertNextCell(triangle)
        polydata = vtk.vtkPolyData()
        polydata.SetPoints(points)
        polydata.SetPolys(cells)
        quality = vtk.vtkMeshQuality()
        quality.SetInputData(polydata)
        report = mesh_quality(mesh)
        per_measure: dict[str, object] = {}
        for metric, (setter, sense) in measures.items():
            getattr(quality, setter)()
            quality.Update()
            array = quality.GetOutput().GetFieldData().GetArray("Mesh Triangle Quality")
            if array is None:
                per_measure[metric] = {"ran": False, "reason": "no quality array"}
                continue
            count = array.GetNumberOfTuples()
            values = [array.GetTuple1(i) for i in range(count)]
            reference = min(values) if sense == "min" else max(values)
            ours = (
                report.summary[metric].minimum if sense == "min" else report.summary[metric].maximum
            )
            per_measure[metric] = {
                "vtk": reference,
                "library": ours,
                "difference": abs(reference - ours),
            }
        out["per_mesh"][name] = per_measure  # type: ignore[index]
    return out


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument(
        "--json-dir", default=str(ROOT / "artifacts" / "mesh_quality_lib"), help="artifact directory"
    )
    parser.add_argument("--quiet", action="store_true", help="only the verdict line")
    args = parser.parse_args(argv)
    outdir = Path(args.json_dir)
    outdir.mkdir(parents=True, exist_ok=True)

    failures: list[str] = []

    # 1. closed form ---------------------------------------------------------------
    analytic: dict[str, dict[str, object]] = {}
    for name, triangle in _canonical_triangles().items():
        measured = triangle_quality(*triangle)
        oracle = quality_oracle.oracle_metrics(*triangle)
        # The float route uses Heron's formula, which loses digits exactly where a sliver
        # is: on a height-1e-4 triangle it reports a 1e-8 relative area error and the
        # resulting aspect-ratio error is 1e-5.  Where the coordinates are rational the
        # exact route settles which of the two is wrong, so it is the reference there and
        # the float route is kept as a reported diagnostic.
        deviations = {
            metric: abs(measured.value(metric) - oracle[metric])
            for metric in quality_oracle.METRICS
            if math.isfinite(oracle[metric]) and math.isfinite(measured.value(metric))
        }
        from fractions import Fraction as _Fraction

        rational = all(
            isinstance(coordinate, (int, float)) and float(_Fraction(coordinate)) == coordinate
            for point in triangle
            for coordinate in point
        )
        exact_reference = (
            quality_oracle.oracle_metrics_exact(*triangle) if rational else None
        )
        if exact_reference is not None:
            exact_deviations = {
                metric: abs(measured.value(metric) - exact_reference[metric])
                for metric in quality_oracle.METRICS
                if metric in exact_reference and math.isfinite(measured.value(metric))
            }
            worst = max(exact_deviations.values(), default=0.0)
            offending = {
                metric: value
                for metric, value in exact_deviations.items()
                if value > _tolerance(metric, EXACT_TOLERANCE)
            }
            if offending:
                failures.append(f"analytic-exact/{name}: {offending}")
        else:
            worst = max(deviations.values(), default=0.0)
            if worst > ROUTE_TOLERANCE:
                failures.append(f"analytic/{name}: deviation {worst}")
        analytic[name] = {
            "measured": measured.to_dict(),
            "independent_float": {k: v for k, v in oracle.items() if k != "degenerate"},
            "independent_exact": exact_reference,
            "max_abs_deviation_vs_exact": worst,
            "max_abs_deviation_vs_float_route": max(deviations.values(), default=0.0),
        }
    ideal = triangle_quality(*_canonical_triangles()["equilateral"])
    ideal_check = {
        metric: {"expected": expected, "got": ideal.value(metric), "difference": abs(ideal.value(metric) - expected)}
        for metric, expected in IDEAL_VALUES.items()
    }
    if max(item["difference"] for item in ideal_check.values()) > ANALYTIC_TOLERANCE:
        failures.append("ideal-values: equilateral triangle does not score the ideal values")

    # 2/3. differential on real and synthetic meshes --------------------------------
    meshes = {**_hull_meshes(), "box_unit": _box_mesh(), "sliver_patch": _sliver_patch()}
    differential: dict[str, object] = {}
    timings: dict[str, float] = {}
    for name, mesh in meshes.items():
        started = time.time()
        report = mesh_quality(mesh)
        timings[name] = time.time() - started
        oracle_faces = [
            quality_oracle.oracle_metrics(
                mesh.vertices[face[0]], mesh.vertices[face[1]], mesh.vertices[face[2]]
            )
            for face in mesh.faces
        ]
        per_metric: dict[str, dict[str, float]] = {}
        for metric in METRIC_DEFINITIONS:
            comparison = METRIC_COMPARISON[metric]
            quantity = report.summary[metric]
            reference_values = [item[metric] for item in oracle_faces]
            reference = (
                min(reference_values) if comparison == "min" else max(reference_values)
            )
            ours = quantity.minimum if comparison == "min" else quantity.maximum
            difference = abs(reference - ours)
            per_metric[metric] = {"library": ours, "independent": reference, "difference": difference}
            if math.isfinite(difference) and difference > ROUTE_TOLERANCE:
                failures.append(f"{name}/{metric}: {difference:.3e} > {ROUTE_TOLERANCE:g}")
        differential[name] = {
            "triangles": len(mesh.faces),
            "minutes": round(timings[name] / 60.0, 4),
            "metrics": per_metric,
        }

    exact: dict[str, object] = {}
    for name in ("box_unit", "sliver_patch"):
        mesh = meshes[name]
        worst = 0.0
        detail: list[dict[str, object]] = []
        for index, face in enumerate(mesh.faces):
            measured = triangle_quality(
                mesh.vertices[face[0]], mesh.vertices[face[1]], mesh.vertices[face[2]], index=index
            )
            exact_values = quality_oracle.oracle_metrics_exact(
                mesh.vertices[face[0]], mesh.vertices[face[1]], mesh.vertices[face[2]]
            )
            element = {}
            for metric in quality_oracle.METRICS:
                if metric not in exact_values:
                    continue
                difference = abs(measured.value(metric) - exact_values[metric])
                element[metric] = difference
                worst = max(worst, difference)
                if difference > _tolerance(metric, EXACT_TOLERANCE):
                    failures.append(
                        f"exact/{name}/face{index}/{metric}: {difference:.3e} > "
                        f"{_tolerance(metric, EXACT_TOLERANCE):g}"
                    )
            detail.append({"index": index, "deviations": element})
        exact[name] = {"max_abs_deviation": worst, "elements": detail}

    # 4. upstream ------------------------------------------------------------------
    upstream = _probe_upstream()
    vtk_result = _vtk_compare(meshes)
    upstream["vtk"]["comparison"] = vtk_result

    # 5. dialect audit -------------------------------------------------------------
    dialect: dict[str, object] = {}
    for name, mesh in meshes.items():
        report = mesh_quality(mesh)
        validation = validation_module.measure_quality(mesh)
        facet = gate_module._facet_aspect(mesh)
        parametric = mesh_quality_report(mesh)
        conditioning = conditioning_measure(mesh)
        dialect[name] = {
            "library_max_aspect_ratio": report.summary["aspect_ratio"].maximum,
            "library_max_edge_ratio": report.summary["edge_ratio"].maximum,
            "library_min_radius_ratio_aspect": 1.0 / report.summary["radius_ratio"].minimum
            if report.summary["radius_ratio"].minimum > 0.0
            else math.inf,
            "validation_max_aspect_lmax2_over_2A": validation.max_aspect_ratio,
            "gate_facet_aspect_lmax2_over_A": facet,
            "gate_over_validation": facet / validation.max_aspect_ratio
            if validation.max_aspect_ratio
            else None,
            "parametric_max_aspect_lmax_over_lmin": parametric["max_aspect_ratio"],
            "parametric_min_scaled_jacobian": parametric["min_scaled_jacobian"],
            "validation_min_angle_deg": validation.min_angle_degrees,
            "library_min_angle_deg": report.summary["min_angle_deg"].minimum,
            "mesh_conditioning_measurements": {
                name: getattr(conditioning, name)
                for name in dir(conditioning)
                if not name.startswith("_") and isinstance(getattr(conditioning, name), (int, float))
            },
        }

    # 6. gate and its sensitivity --------------------------------------------------
    gate_results: dict[str, object] = {}
    for name, mesh in meshes.items():
        verdict = quality_gate(mesh)
        gate_results[name] = {
            "passed": verdict.passed,
            "checks": [check.to_dict() for check in verdict.checks],
        }
    sweep = {
        "min_angle_deg": {
            str(limit): {
                name: quality_gate(mesh, criteria={"min_angle_deg": limit}).passed
                for name, mesh in meshes.items()
            }
            for limit in (0.5, 1.0, 2.0, 5.0, 10.0)
        },
        "aspect_ratio": {
            str(limit): {
                name: quality_gate(mesh, criteria={"aspect_ratio": limit}).passed
                for name, mesh in meshes.items()
            }
            for limit in (10.0, 25.0, 50.0, 100.0)
        },
    }

    packet = {
        "generated_by": "tools/mesh_quality_crosscheck.py",
        "criteria": {
            "route_tolerance": ROUTE_TOLERANCE,
            "exact_tolerance": EXACT_TOLERANCE,
            "analytic_tolerance": ANALYTIC_TOLERANCE,
            "default_gate": DEFAULT_CRITERIA,
        },
        "ideal_values": ideal_check,
        "analytic_elements": analytic,
        "differential": differential,
        "exact_arithmetic": exact,
        "upstream": upstream,
        "dialect_audit": dialect,
        "mesh_quality_seconds": timings,
        "failures": failures,
        "passed": not failures,
    }
    (outdir / "crosscheck.json").write_text(
        json.dumps(packet, indent=2, sort_keys=True, default=str) + "\n"
    )
    (outdir / "gate.json").write_text(
        json.dumps({"gate": gate_results, "sensitivity": sweep}, indent=2, sort_keys=True, default=str)
        + "\n"
    )

    if not args.quiet:
        print(f"{'mesh':16s} {'tris':>7s} {'worst diff':>12s} {'min angle':>10s} "
              f"{'aspect':>9s} {'gamma':>10s} {'s':>7s}")
        for name, item in differential.items():
            worst = max(
                (value["difference"] for value in item["metrics"].values()),  # type: ignore[index]
                key=lambda value: (math.inf if not math.isfinite(value) else value),
                default=0.0,
            )
            report = mesh_quality(meshes[name])
            print(
                f"{name:16s} {item['triangles']:7d} {worst:12.3e} "
                f"{report.summary['min_angle_deg'].minimum:10.4f} "
                f"{report.summary['aspect_ratio'].maximum:9.4f} "
                f"{report.summary['gamma'].minimum:10.6f} {timings[name]:7.3f}"
            )
        print()
        print("upstream oracles:", {name: item["available"] for name, item in upstream.items() if name in ("vtk", "gmsh", "meshmagick")})
        print("gates:", {name: item["passed"] for name, item in gate_results.items()})
        print()
        if failures:
            print("FAILURES:")
            for line in failures:
                print("  -", line)
        else:
            print(
                f"PASS: every measure within {ROUTE_TOLERANCE:g} of the independent route, "
                f"{EXACT_TOLERANCE:g} of exact arithmetic on rational meshes, and "
                f"{ANALYTIC_TOLERANCE:g} of the closed forms."
            )
        print(f"artifacts: {outdir/'crosscheck.json'}, {outdir/'gate.json'}")
    return 0 if not failures else 1


if __name__ == "__main__":
    raise SystemExit(main())