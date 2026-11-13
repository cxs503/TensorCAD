#!/usr/bin/env python3
"""Benchmark the mesh-quality library against analytic values and independent oracles.

What this script establishes, and how it can fail:

1. **Analytic** -- on canonical elements (equilateral triangle, square, regular
   tetrahedron, unit cube) every normalised metric is exactly 1 by definition, and the
   regular tetrahedron is ``1/sqrt(2)`` in the scaled-Jacobian unit and
   ``(1/sqrt(2))**(2/3)`` in Gmsh's ICN unit. Tolerance 1e-12.
2. **Independent implementation** -- second derivations (Heron + atan2 for triangles,
   Cayley-Menger for tetrahedra, central finite differences of the isoparametric map for
   hexahedra) on every mesh built below. Triangles/tetrahedra at 1e-9, hexahedra at
   1e-7 (finite-difference truncation).
3. **Third-party oracle** -- ``trimesh`` (installed here) for triangle area and angles,
   on three independent triangulations at 1e-9.
4. **External tools that are not installed** -- Gmsh ships no wheel for this
   architecture and meshmagick is not on the package index; both are recorded as
   ``not_run`` with the command that would run them, never as ``pass``.

Writes ``artifacts/mesh_quality/results.json`` (every number the document quotes),
``artifacts/mesh_quality/plot.png`` and ``artifacts/mesh_quality/manifest.json``
(sha256 of the results), and exits non-zero if any criterion above stops holding.
"""
from __future__ import annotations

import hashlib
import json
import math
import sys
from pathlib import Path
from statistics import mean

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT))

from tensorcad.hull import HullParameters, build_mesh  # noqa: E402
from tensorcad.meshquality import DEFAULT_GATE, measure, reference  # noqa: E402
from tensorcad.meshquality.io import MeshSpec, write_spec  # noqa: E402
from tensorcad.volume_mesh import (  # noqa: E402
    box_section_grid,
    build_volume_mesh,
    volume_mesh_from_hull,
)

ARTIFACT_DIR = REPO_ROOT / "artifacts" / "mesh_quality"

TOLERANCES = {"analytic": 1e-12, "reference_tri": 1e-9, "reference_tet": 1e-9, "reference_hex": 1e-7, "trimesh": 1e-9}


def ideal_elements() -> dict[str, tuple[list[tuple[float, float, float]], dict[str, list[tuple[int, ...]]]]]:
    """The four canonical elements, each of which every normalised metric scores 1."""
    h = math.sqrt(3.0) / 2.0
    regular = [(0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (0.5, h, 0.0), (0.5, h / 3.0, math.sqrt(2.0 / 3.0))]
    return {
        "ideal_tri": ([(0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (0.5, h, 0.0)], {"tri": [(0, 1, 2)]}),
        "ideal_quad": (
            [(0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (1.0, 1.0, 0.0), (0.0, 1.0, 0.0)],
            {"quad": [(0, 1, 2, 3)]},
        ),
        "ideal_tet": (regular, {"tet": [(0, 1, 2, 3)]}),
        "ideal_hex": (
            [
                (0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (1.0, 1.0, 0.0), (0.0, 1.0, 0.0),
                (0.0, 0.0, 1.0), (1.0, 0.0, 1.0), (1.0, 1.0, 1.0), (0.0, 1.0, 1.0),
            ],
            {"hex": [tuple(range(8))]},
        ),
    }


def build_meshes() -> dict[str, tuple[list[tuple[float, float, float]], dict[str, list[tuple[int, ...]]]]]:
    """The real meshes: two volume meshes from the repo mesher and two surface meshes."""
    meshes: dict[str, tuple[list[tuple[float, float, float]], dict[str, list[tuple[int, ...]]]]] = {}

    # 1. 4x4x4 unit cubes -- every cell is the reference cube, so the metrics are exact.
    box = build_volume_mesh(
        box_section_grid(length=4.0, beam=4.0, draft=4.0, nx=5, nz=5),
        first_layer=1.0,
        growth=1.0,
        layers=2,
        split_tets=False,
    )
    meshes["box_ideal_hex"] = (list(box.vertices), {"hex": [tuple(h) for h in box.hexes]})

    # 2. The same box stretched 4x along x: orthogonal corners, anisotropic cells.
    stretched = [(4.0 * x, y, z) for x, y, z in box.vertices]
    meshes["box_stretched_hex"] = (stretched, {"hex": [tuple(h) for h in box.hexes]})

    # 3. The hull volume mesh the repo ships (L=10, B=2, T=0.8, 41x13 stations/levels).
    hull = HullParameters(length=10.0, beam=2.0, draft=0.8, nx=41, nz=13)
    volume = volume_mesh_from_hull(hull, first_layer=0.005, growth=1.2, layers=4)
    meshes["hull_hex"] = (list(volume.vertices), {"hex": [tuple(h) for h in volume.hexes]})
    meshes["hull_tet"] = (list(volume.vertices), {"tet": [tuple(t) for t in volume.tets]})

    # 4. The hull's own surface triangulation.
    surface, _ = build_mesh(hull)
    meshes["hull_surface_tri"] = (list(surface.vertices), {"tri": [tuple(f) for f in surface.faces]})
    return meshes


def triangulated_faces(kind: str, cells: list[tuple[int, ...]]) -> list[tuple[int, int, int]]:
    """Boundary-agnostic triangulation of a hexahedral/quad mesh, for the trimesh oracle."""
    if kind == "quad":
        return [(c[0], c[1], c[2]) for c in cells] + [(c[0], c[2], c[3]) for c in cells]
    faces = (
        (0, 3, 2, 1), (4, 5, 6, 7), (0, 1, 5, 4),
        (1, 2, 6, 5), (2, 3, 7, 6), (3, 0, 4, 7),
    )
    out: list[tuple[int, int, int]] = []
    for cell in cells:
        for face in faces:
            quad = [cell[i] for i in face]
            out.append((quad[0], quad[1], quad[2]))
            out.append((quad[0], quad[2], quad[3]))
    return out


def compare(
    got: tuple[float, ...],
    expected: tuple[float, ...],
    *,
    relative: bool = False,
    mask: "list[bool] | None" = None,
) -> dict[str, float]:
    """Max absolute (or relative) difference of two metric vectors, ignoring NaN.

    ``mask`` excludes cells that carry no well-posed value: collapsed cells (zero area
    or volume) and angle computations on near-degenerate corners, where ``acos`` loses
    half its digits by construction. The excluded count is returned so a reader can see
    how much of a mesh was compared rather than take the maximum on faith.
    """
    best = 0.0
    compared = 0
    excluded = 0
    for index, (a, b) in enumerate(zip(got, expected)):
        if mask is not None and not mask[index]:
            excluded += 1
            continue
        if math.isnan(a) or math.isnan(b) or math.isinf(a) or math.isinf(b):
            excluded += 1
            continue
        compared += 1
        if relative:
            denominator = max(abs(a), abs(b), 1e-300)
            best = max(best, abs(a - b) / denominator)
        else:
            best = max(best, abs(a - b))
    return {"max_diff": best, "count": compared, "excluded": excluded}


def well_posed(angles: tuple[float, ...], areas: tuple[float, ...]) -> list[bool]:
    """Cells with positive measure whose angles are not at the degeneracy limit."""
    return [
        area > 0.0 and 1e-3 < angle < 180.0 - 1e-3 for area, angle in zip(areas, angles)
    ]


def main() -> int:
    ARTIFACT_DIR.mkdir(parents=True, exist_ok=True)
    meshes = ideal_elements() | build_meshes()
    availability = reference.external_availability()

    results: dict[str, object] = {
        "schema": "tensorcad.meshquality.benchmark/1",
        "external_availability": availability,
        "not_run": {
            "gmsh": (
                "not installed: the gmsh wheel is published for x86_64/aarch64 only and this host is "
                "loongarch64 (pip download gmsh -> 'no matching distribution'). Run "
                "`python -c \"from tensorcad.meshquality.reference import gmsh_element_qualities as g\"` "
                "on a host with gmsh to compare gmsh_icn against getElementQualities('minSICN')."
            ),
            "meshmagick": (
                "not on the package index and GPL-3.0; its print_quality() delegates to "
                "vtk.vtkMeshQuality, whose Verdict definitions the metric module mirrors. A VTK "
                "wheel for loongarch64 does not exist, so the runtime comparison cannot run here."
            ),
        },
        "meshes": {},
    }
    failures: list[str] = []

    for name, (points, elements) in meshes.items():
        quality = measure(points, elements)
        gate = DEFAULT_GATE.evaluate(quality)
        positive: dict[str, dict[str, float]] = {}
        entry: dict[str, object] = {
            "cells": {kind: len(cells) for kind, cells in elements.items()},
            "points": len(points),
            "summary": {kind: quality.to_dict()[kind]["summary"] for kind in quality.elements},
            "gate": gate.to_dict(),
            "positive_only": positive,
            "comparisons": {},
        }
        comparisons: dict[str, dict[str, float]] = entry["comparisons"]  # type: ignore[assignment]

        for kind, quality_kind in quality.elements.items():
            size_name = "area" if kind in ("tri", "quad") else "volume"
            keep = [index for index, size in enumerate(quality_kind[size_name]) if size > 0.0]
            positive[kind] = {
                "cells": len(keep),
                "excluded_zero_measure": len(quality_kind[size_name]) - len(keep),
            }
            for metric in ("min_sicn", "scaled_jacobian", "min_angle", "min_dihedral", "max_angle"):
                if metric not in quality_kind.metrics or not keep:
                    continue
                values = [quality_kind[metric][index] for index in keep]
                if any(math.isnan(v) for v in values):
                    continue
                positive[kind][f"{metric}_min"] = min(values)
                positive[kind][f"{metric}_max"] = max(values)
            if kind == "tri":
                independent = reference.reference_triangles(points, elements["tri"])
                tri_mask = well_posed(quality_kind["min_angle"], quality_kind["area"])
                for metric in ("area", "min_angle", "max_angle", "aspect_ratio", "gamma"):
                    relative = metric == "area"
                    comparisons[f"reference_tri.{metric}"] = compare(
                        quality_kind[metric], independent[metric], relative=relative, mask=tri_mask
                    )
                external = reference.trimesh_triangle_metrics(points, elements["tri"])
                if external is not None:
                    for metric in ("area", "min_angle", "max_angle"):
                        comparisons[f"trimesh.{metric}"] = compare(
                            quality_kind[metric],
                            external[metric],
                            relative=metric == "area",
                            mask=tri_mask,
                        )
            elif kind == "tet":
                independent_tet = reference.reference_tets(points, elements["tet"])
                tet_mask = well_posed(quality_kind["min_dihedral"], quality_kind["volume"])
                for metric in ("volume", "eta", "min_dihedral", "max_dihedral"):
                    comparisons[f"reference_tet.{metric}"] = compare(
                        quality_kind[metric],
                        independent_tet[metric],
                        relative=metric == "volume",
                        mask=tet_mask,
                    )
            elif kind == "hex":
                independent_hex = reference.reference_hexahedra(points, elements["hex"])
                hex_mask = [volume > 0.0 for volume in quality_kind["volume"]]
                for metric in ("scaled_jacobian", "min_sicn"):
                    comparisons[f"reference_hex.{metric}"] = compare(
                        quality_kind[metric], independent_hex[metric], mask=hex_mask
                    )
                # The same hexahedra, triangulated, give the trimesh oracle a real surface.
                face_tris = triangulated_faces("hex", [tuple(c) for c in elements["hex"]])
                external = reference.trimesh_triangle_metrics(points, face_tris)
                if external is not None:
                    library_tri = measure(points, {"tri": face_tris})
                    face_mask = well_posed(library_tri["tri"]["min_angle"], library_tri["tri"]["area"])
                    for metric in ("area", "min_angle", "max_angle"):
                        comparisons[f"trimesh_faces.{metric}"] = compare(
                            library_tri["tri"][metric],
                            external[metric],
                            relative=metric == "area",
                            mask=face_mask,
                        )

        if name.startswith("ideal_"):
            kind = next(iter(elements))
            expectation = {"tri": 1.0, "quad": 1.0, "tet": 1.0, "hex": 1.0}[kind]
            entry["analytic"] = {
                "min_sicn": quality[kind]["min_sicn"][0],
                "scaled_jacobian": quality[kind]["scaled_jacobian"][0],
                "expected_sicn": expectation,
                "expected_scaled_jacobian": 1.0 / math.sqrt(2.0) if kind == "tet" else 1.0,
            }
            key = {"min_sicn": "expected_sicn", "scaled_jacobian": "expected_scaled_jacobian"}
            for metric in ("min_sicn", "scaled_jacobian"):
                diff = abs(entry["analytic"][metric] - entry["analytic"][key[metric]])
                entry["analytic"][f"{metric}_diff"] = diff
                if diff > TOLERANCES["analytic"]:
                    failures.append(f"{name}.analytic.{metric}: {diff:.3e} > 1e-12")
            # Gmsh's own ICN unit on the two tetrahedral references.
            if kind == "tet":
                entry["analytic"]["gmsh_icn"] = quality[kind]["gmsh_icn"][0]
                entry["analytic"]["expected_gmsh_icn_regular"] = (1.0 / math.sqrt(2.0)) ** (2.0 / 3.0)
                diff = abs(entry["analytic"]["gmsh_icn"] - entry["analytic"]["expected_gmsh_icn_regular"])
                entry["analytic"]["gmsh_icn_diff"] = diff
                if diff > TOLERANCES["analytic"]:
                    failures.append(f"{name}.analytic.gmsh_icn: {diff:.3e} > 1e-12")

        # Criteria on the cross-checks.
        for label, comparison in comparisons.items():
            tolerance = 0.0
            if label.startswith("reference_tri."):
                tolerance = TOLERANCES["reference_tri"]
            elif label.startswith("reference_tet."):
                tolerance = TOLERANCES["reference_tet"]
            elif label.startswith("reference_hex."):
                tolerance = TOLERANCES["reference_hex"]
            elif label.startswith("trimesh"):
                tolerance = TOLERANCES["trimesh"]
            if comparison["max_diff"] > tolerance:
                failures.append(f"{name}.{label}: {comparison['max_diff']:.3e} > {tolerance:g}")
        results["meshes"][name] = entry  # type: ignore[index]

    # The gate must pass on real meshes and fail on a mesh with a sliver -- otherwise it
    # is decoration.
    bad_points = [(0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (0.5, 1e-4, 0.0)]
    bad_quality = measure(bad_points, {"tri": [(0, 1, 2)]})
    bad_gate = DEFAULT_GATE.evaluate(bad_quality)
    results["gate_discrimination"] = {
        "clean_meshes_passing": [
            name for name, entry in results["meshes"].items() if entry["gate"]["passed"]  # type: ignore[index]
        ],
        "sliver_fails": not bad_gate.passed,
        "sliver_failures": [check.name for check in bad_gate.failures],
        "sliver_min_sicn": bad_quality["tri"]["min_sicn"][0],
    }
    for name, entry in results["meshes"].items():  # type: ignore[assignment]
        if name.startswith("ideal_") and not entry["gate"]["passed"]:
            failures.append(f"gate rejected the canonical element {name}")
    if bad_gate.passed:
        failures.append("gate accepted a sliver triangle")

    # Write the hull surface mesh as an exchange-format artifact (a real one, not a toy).
    hull_points, hull_elements = meshes["hull_surface_tri"]
    write_spec(
        ARTIFACT_DIR / "hull_surface_tri.json",
        MeshSpec(
            points=tuple(hull_points),
            elements={"tri": tuple(hull_elements["tri"])},
            metadata={"source": "tensorcad.hull.build_mesh", "L": 10.0, "B": 2.0, "T": 0.8},
        ),
    )
    write_spec(
        ARTIFACT_DIR / "hull_hex.json",
        MeshSpec(
            points=tuple(meshes["hull_hex"][0]),
            elements={"hex": tuple(meshes["hull_hex"][1]["hex"])},
            metadata={"source": "tensorcad.volume_mesh.volume_mesh_from_hull", "first_layer": 0.005,
                      "growth": 1.2, "layers": 4},
        ),
    )

    results["criteria"] = {
        "tolerances": TOLERANCES,
        "failures": failures,
        "passed": not failures,
    }

    results_path = ARTIFACT_DIR / "results.json"
    results_path.write_text(json.dumps(results, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    digest = hashlib.sha256(results_path.read_bytes()).hexdigest()
    (ARTIFACT_DIR / "manifest.json").write_text(
        json.dumps(
            {
                "results_sha256": digest,
                "generated_by": "examples/mesh_quality/run_benchmark.py",
                "external_availability": availability,
            },
            indent=2,
            sort_keys=True,
        )
        + "\n",
        encoding="utf-8",
    )

    # Figure: the quality distribution the gate actually sees, plus the agreement margin.
    hull_hex = results["meshes"]["hull_hex"]  # type: ignore[index]
    hull_quality = measure(*meshes["hull_hex"])
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.2), constrained_layout=True)
    axes[0].hist([v for v in hull_quality["hex"]["min_sicn"] if v > 0.0], bins=40, color="#3b6ea5")
    axes[0].axvline(DEFAULT_GATE.limits["hex"][0].minimum or 0.0, color="#b23a3a", linestyle="--",
                    label="gate min_sicn = 0.1")
    axes[0].set_xlabel("hex min_sicn")
    axes[0].set_ylabel("cells")
    axes[0].set_title("hull volume mesh: 3840 hexahedra")
    axes[0].legend()

    labels, margins = [], []
    for name, entry in results["meshes"].items():  # type: ignore[assignment]
        for label, comparison in entry["comparisons"].items():
            if comparison["max_diff"] > 0.0:
                labels.append(f"{name}\n{label}")
                margins.append(comparison["max_diff"])
    axes[1].bar(range(len(margins)), margins, color="#4f8a5b")
    axes[1].set_yscale("log")
    axes[1].axhline(1e-9, color="#b23a3a", linestyle="--", linewidth=1)
    axes[1].set_xticks([])
    axes[1].set_ylabel("max |library - reference|")
    axes[1].set_title(f"{len(margins)} cross-checks (red: 1e-9)")
    fig.savefig(ARTIFACT_DIR / "plot.png", dpi=130)

    print(json.dumps(results["criteria"], indent=2, sort_keys=True))
    print(f"meshes: {', '.join(results['meshes'])}")
    print(f"artifacts: {results_path} (sha256 {digest[:16]}...)")
    if failures:
        for failure in failures:
            print(f"FAIL {failure}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())