"""A small deterministic parametric hull-surface generator.

This first example generates a sampled, symmetric displacement-hull surface.
It is a geometry/API smoke test, not a validated naval-architecture hull form.
Only the Python standard library is required.
"""
from __future__ import annotations

import csv
import json
import math
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any


@dataclass(frozen=True)
class HullParameters:
    """Principal dimensions in metres and sampling resolution."""

    length: float = 10.0
    beam: float = 2.0
    draft: float = 0.8
    nx: int = 41
    nz: int = 13

    def validate(self) -> None:
        for name in ("length", "beam", "draft"):
            value = getattr(self, name)
            if not math.isfinite(value) or value <= 0:
                raise ValueError(f"{name} must be a finite positive number")
        if self.nx < 3 or self.nz < 3:
            raise ValueError("nx and nz must both be >= 3")


@dataclass(frozen=True)
class HullModel:
    parameters: HullParameters
    # Vertex coordinates in metres. Each side is sampled on an x-z grid.
    vertices: tuple[tuple[float, float, float], ...]
    # OBJ-style 1-based triangular face indices.
    faces: tuple[tuple[int, int, int], ...]
    offsets: tuple[tuple[float, float], ...]
    coordinate_frame: str = "x=longitudinal, y=transverse, z=vertical; origin at aft waterline centre"
    units: str = "m"

    def metadata(self) -> dict[str, Any]:
        return {
            "schema_version": "1.0",
            "model_type": "sampled_parametric_displacement_hull",
            "parameters": asdict(self.parameters),
            "units": self.units,
            "coordinate_frame": self.coordinate_frame,
            "vertex_count": len(self.vertices),
            "triangle_count": len(self.faces),
            "offset_count": len(self.offsets),
            "validation": validate_hull(self),
            "engineering_status": (
                "geometry demonstration only; hydrostatics, fairness, stability, "
                "resistance, and classification-rule checks are not implemented"
            ),
        }


def generate_hull(parameters: HullParameters = HullParameters()) -> HullModel:
    """Generate a symmetric sampled hull surface from principal dimensions.

    The half-breadth distribution is analytic and intentionally simple:
    breadth tapers to zero at both ends and grows from keel to waterline.
    The resulting surface is useful for exercising reproducible geometry
    workflows, but is not intended as a production hull design.
    """
    parameters.validate()
    p = parameters
    vertices: list[tuple[float, float, float]] = []
    offsets: list[tuple[float, float]] = []

    # Two independent surfaces (port and starboard) make the output easy to
    # consume without a CAD kernel. Coordinates are generated deterministically.
    for side in (-1.0, 1.0):
        for i in range(p.nx):
            x = p.length * i / (p.nx - 1)
            longitudinal = max(0.0, math.sin(math.pi * i / (p.nx - 1))) ** 0.7
            for k in range(p.nz):
                vertical_fraction = k / (p.nz - 1)
                z = -p.draft + p.draft * vertical_fraction
                half_breadth = 0.5 * p.beam * longitudinal * math.sqrt(vertical_fraction)
                y = side * half_breadth
                vertices.append((x, y, z))
                if side == 1.0:
                    offsets.append((x, 2.0 * half_breadth))

    faces: list[tuple[int, int, int]] = []
    grid_size = p.nx * p.nz
    for side_index in range(2):
        base = side_index * grid_size
        for i in range(p.nx - 1):
            for k in range(p.nz - 1):
                a = base + i * p.nz + k + 1
                b = base + (i + 1) * p.nz + k + 1
                c = base + (i + 1) * p.nz + (k + 1) + 1
                d = base + i * p.nz + (k + 1) + 1
                faces.extend(((a, b, c), (a, c, d)))

    model = HullModel(
        parameters=p,
        vertices=tuple(vertices),
        faces=tuple(faces),
        offsets=tuple(offsets),
    )
    errors = validate_hull(model)
    if errors:
        raise ValueError("generated hull failed validation: " + "; ".join(errors))
    return model


def validate_hull(model: HullModel) -> list[str]:
    """Return validation errors; an empty list means basic checks passed."""
    errors: list[str] = []
    p = model.parameters
    try:
        p.validate()
    except ValueError as exc:
        errors.append(str(exc))
    if not model.vertices:
        errors.append("no vertices generated")
    if any(len(v) != 3 or not all(math.isfinite(c) for c in v) for v in model.vertices):
        errors.append("vertex coordinates must be finite 3D points")
    if len(model.vertices) != 2 * p.nx * p.nz:
        errors.append("unexpected vertex count")
    if len(model.offsets) != p.nx:
        errors.append("unexpected longitudinal offset count")
    if any(width < -1e-12 or not math.isfinite(width) for _, width in model.offsets):
        errors.append("offset breadths must be finite and non-negative")
    if any(len(face) != 3 or min(face) < 1 or max(face) > len(model.vertices) for face in model.faces):
        errors.append("face indices must be valid 1-based triangle indices")
    expected_faces = 4 * (p.nx - 1) * (p.nz - 1)
    if len(model.faces) != expected_faces:
        errors.append("unexpected triangle count")
    return errors


def export_hull(model: HullModel, output_dir: str | Path) -> dict[str, str]:
    """Write OBJ mesh, longitudinal offsets CSV, and reproducibility metadata."""
    output = Path(output_dir)
    output.mkdir(parents=True, exist_ok=True)
    obj_path = output / "hull.obj"
    csv_path = output / "offsets.csv"
    metadata_path = output / "metadata.json"

    with obj_path.open("w", encoding="utf-8", newline="\n") as stream:
        stream.write("# TensorCAD sampled parametric hull; units: metres\n")
        for x, y, z in model.vertices:
            stream.write(f"v {x:.9g} {y:.9g} {z:.9g}\n")
        for a, b, c in model.faces:
            stream.write(f"f {a} {b} {c}\n")

    with csv_path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(("x_m", "overall_breadth_m"))
        writer.writerows((f"{x:.9g}", f"{width:.9g}") for x, width in model.offsets)

    with metadata_path.open("w", encoding="utf-8", newline="\n") as stream:
        json.dump(model.metadata(), stream, indent=2, sort_keys=True)
        stream.write("\n")

    return {"mesh": str(obj_path), "offsets": str(csv_path), "metadata": str(metadata_path)}
