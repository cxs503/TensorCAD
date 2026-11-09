"""Boolean field composition: closed-form volumes, convergence, SUBOFF equivalence.

Run from the repository root::

    python examples/csg_boolean/run_example.py

Writes ``docs/csg_boolean/results.json`` and ``docs/csg_boolean/plot.png``. Every number
in the results file is computed here from the library; nothing is transcribed.
"""
from __future__ import annotations

import json
import math
import pathlib

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

from tensorcad import csg  # noqa: E402
from tensorcad.geometry_field import (  # noqa: E402
    CylinderField,
    PlaneField,
    SphereField,
    to_mesh,
)
from tensorcad.mesh import total_signed_volume  # noqa: E402
from tensorcad.suboff_appendages import SuboffFullField  # noqa: E402
from tensorcad.voxel import to_voxel_grid  # noqa: E402

OUT = pathlib.Path(__file__).resolve().parents[2] / "docs" / "csg_boolean"
RESOLUTIONS = (48, 64, 96, 128)


def sphere_cap_volume(radius: float, height: float) -> float:
    return math.pi * height**2 * (3.0 * radius - height) / 3.0


def sphere_zone_volume(radius: float, half_height: float) -> float:
    return 2.0 * math.pi * (radius**2 * half_height - half_height**3 / 3.0)


def bicylinder_volume(radius: float) -> float:
    return 16.0 * radius**3 / 3.0


def sphere_lens_volume(radius: float, separation: float) -> float:
    return math.pi * (4.0 * radius + separation) * (2.0 * radius - separation) ** 2 / 12.0


def shell_volume(outer: float, inner: float) -> float:
    return 4.0 * math.pi / 3.0 * (outer**3 - inner**3)


def ring_volume(radius: float, hole: float) -> float:
    """Ball minus a coaxial through-cylinder: ``4 pi/3 (R^2 - r^2)^{3/2}``."""
    return 4.0 / 3.0 * math.pi * (radius**2 - hole**2) ** 1.5


def cases() -> dict[str, tuple[object, float, str]]:
    radius = 1.0
    sphere = SphereField(radius=radius, center=(0.0, 0.0, 0.0))
    top = PlaneField(normal=(0.0, 0.0, 1.0), offset=0.5)
    upper = csg.difference(sphere, top)
    zone = csg.intersection(
        sphere,
        PlaneField(normal=(0.0, 0.0, 1.0), offset=0.4),
        PlaneField(normal=(0.0, 0.0, -1.0), offset=0.4),
    )
    r = 0.6
    bicylinder = csg.intersection(
        CylinderField(radius=r, height=4.0, axis=(0.0, 0.0, 1.0)),
        CylinderField(radius=r, height=4.0, axis=(0.0, 1.0, 0.0)),
    )
    lens = csg.intersection(
        SphereField(radius=0.75, center=(0.0, 0.0, 0.0)),
        SphereField(radius=0.75, center=(1.0, 0.0, 0.0)),
    )
    shell = csg.difference(
        SphereField(radius=1.0, center=(0.0, 0.0, 0.0)),
        SphereField(radius=0.6, center=(0.0, 0.0, 0.0)),
    )
    ring = csg.difference(
        SphereField(radius=0.6, center=(0.0, 0.0, 0.0)),
        CylinderField(radius=0.3, height=2.0, axis=(0.0, 0.0, 1.0)),
    )
    return {
        "cap (sphere R=1 minus half-space z>=0.5)": (
            upper, sphere_cap_volume(1.0, 0.5), "pi h^2 (3R-h)/3, h = 0.5"),
        "zone (sphere R=1, |z| <= 0.4)": (
            zone, sphere_zone_volume(1.0, 0.4), "2 pi (R^2 t - t^3/3)"),
        "bicylinder (two cylinders r=0.6, perpendicular)": (
            bicylinder, bicylinder_volume(0.6), "16 r^3 / 3 (Steinmetz)"),
        "lens (two spheres R=0.75, centres 1.0 apart)": (
            lens, sphere_lens_volume(0.75, 1.0), "pi (4R+d)(2R-d)^2/12"),
        "shell (sphere R=1 minus concentric R=0.6)": (
            shell, shell_volume(1.0, 0.6), "4 pi/3 (Ro^3 - Ri^3)"),
        "ring (sphere R=0.6 minus coaxial cylinder r=0.3)": (
            ring, ring_volume(0.6, 0.3), "4 pi/3 (R^2 - r^2)^{3/2}"),
    }


def table() -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    longest = max(len(name) for name in cases())
    for name, (field, reference, formula) in cases().items():
        errors = []
        volumes = {}
        for resolution in RESOLUTIONS:
            volume = to_voxel_grid(field, resolution).volume()
            volumes[str(resolution)] = volume
            errors.append((volume - reference) / reference)
        rows.append(
            {
                "case": name,
                "formula": formula,
                "closed_form_m3": reference,
                "volumes_m3": volumes,
                "relative_error": dict(zip((str(r) for r in RESOLUTIONS), errors)),
                "error_at_128": errors[-1],
                # The partial-fill rule is second order but not strictly monotone
                # case by case (voxel phase noise), so the claim is the trend that
                # matters: the finest grid is closer than the coarsest one.
                "shrinks_from_coarsest": abs(errors[-1]) < abs(errors[0]),
            }
        )
        print(
            f"{name:<{longest}}  closed form {reference:9.6f}  "
            + "  ".join(f"{resolution:>3}:{100 * error:+8.4f} %"
                        for resolution, error in zip(RESOLUTIONS, errors))
        )
    return rows


def suboff_equivalence() -> dict[str, object]:
    full = SuboffFullField()
    composed = csg.union(full.hull(), full.appendages())
    reference = to_voxel_grid(full, 64)
    measured = to_voxel_grid(composed, 64)
    return {
        "volume_m3": reference.volume(),
        "volume_difference_m3": measured.volume() - reference.volume(),
        "sdf_grids_identical": measured.sdf_digest() == reference.sdf_digest(),
        "sdf_sha256": reference.sdf_digest(),
    }


def drilled_mesh() -> dict[str, object]:
    body = csg.difference(
        SphereField(radius=0.6, center=(0.0, 0.0, 0.0)),
        CylinderField(radius=0.3, height=2.0, axis=(0.0, 0.0, 1.0)),
    )
    volumes = {}
    for resolution in (48, 64, 96, 128):
        volumes[str(resolution)] = total_signed_volume(to_mesh(body, resolution).mesh)
    reference = ring_volume(0.6, 0.3)
    return {
        "closed_form_m3": reference,
        "mesh_volumes_m3": volumes,
        "relative_error": {key: (value - reference) / reference for key, value in volumes.items()},
    }


def figure(rows: list[dict[str, object]], drilled: dict[str, object]) -> None:
    fig, axes = plt.subplots(1, 3, figsize=(16.0, 4.6))
    grid = np.linspace(-6.5, 5.0, 460)
    fine = np.linspace(-0.85, 0.85, 340)

    # (a) SUBOFF plan view: union composed here, not SuboffFullField.
    full = SuboffFullField()
    composed = csg.union(full.hull(), full.appendages())
    xx, yy = np.meshgrid(grid, np.zeros_like(grid), indexing="ij")
    xx, zz = np.meshgrid(grid, np.linspace(-0.75, 0.75, 340), indexing="ij")
    plane = np.stack([xx.ravel(), np.zeros(xx.size), zz.ravel()], axis=1)
    values = np.asarray(composed.numpy_callable()(plane)).reshape(xx.shape)
    axes[0].contourf(xx, zz, np.minimum(values, 0.12), levels=24, cmap="viridis")
    axes[0].contour(xx, zz, values, levels=[0.0], colors="white", linewidths=1.1)
    axes[0].set_title("union(hull, appendages), centreplane z = 0\n"
                      "white = zero level set (the surface)")
    axes[0].set_xlabel("x [m] (SUBOFF 4.356 m)")
    axes[0].set_ylabel("z [m]")

    # (b) ball minus through-cylinder: the cut that changes topology.
    body = csg.difference(
        SphereField(radius=0.6, center=(0.0, 0.0, 0.0)),
        CylinderField(radius=0.3, height=2.0, axis=(0.0, 0.0, 1.0)),
    )
    xx, zz = np.meshgrid(fine, fine, indexing="ij")
    plane = np.stack([xx.ravel(), np.zeros(xx.size), zz.ravel()], axis=1)
    values = np.asarray(body.numpy_callable()(plane)).reshape(xx.shape)
    axes[1].contourf(xx, zz, np.clip(values, -0.15, 0.15), levels=25, cmap="coolwarm")
    axes[1].contour(xx, zz, values, levels=[0.0], colors="black", linewidths=1.1)
    axes[1].set_title("difference(ball R=0.6, cylinder r=0.3), y = 0\n"
                      "topology: genus 1 (V-E+F = 0), reported as such")
    axes[1].set_xlabel("x [m]")
    axes[1].set_ylabel("z [m]")
    axes[1].set_aspect("equal")

    # (c) convergence of the measured volume to the closed form.
    resolutions = RESOLUTIONS
    for row in rows:
        errors = [100.0 * row["relative_error"][str(r)] for r in resolutions]
        axes[2].plot(resolutions, errors, marker="o", label=str(row["case"]).split(" (")[0])
    mesh_errors = [100.0 * drilled["relative_error"][str(r)] for r in (48, 64, 96, 128)]
    axes[2].plot((48, 64, 96, 128), mesh_errors, marker="s", ls="--",
                 label="ring, extracted mesh")
    axes[2].axhline(0.0, color="black", linewidth=0.8)
    axes[2].set_xlabel("cells per axis (voxel) / resolution (mesh)")
    axes[2].set_ylabel("volume error vs closed form [%]")
    axes[2].set_title("boolean volumes converge to the closed forms")
    axes[2].legend(fontsize=7, loc="best")
    axes[2].grid(alpha=0.3)

    fig.tight_layout()
    OUT.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT / "plot.png", dpi=140)
    print(f"wrote {OUT / 'plot.png'}")


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    rows = table()
    suboff = suboff_equivalence()
    drilled = drilled_mesh()
    figure(rows, drilled)

    worst = max(abs(row["error_at_128"]) for row in rows)
    claims = {
        "every boolean volume within 1 % of its closed form at 128 cells": worst < 0.01,
        "every case is closer to its closed form at 128 cells than at 48": all(
            row["shrinks_from_coarsest"] for row in rows
        ),
        "composed SUBOFF union reproduces SuboffFullField bit for bit": (
            suboff["sdf_grids_identical"] and suboff["volume_difference_m3"] == 0.0
        ),
        "the drilled solid's extracted mesh is within 2 % of the closed form": (
            abs(drilled["relative_error"]["128"]) < 0.02
        ),
    }
    payload = {
        "schema": "tensorcad.csg_boolean.v1",
        "resolutions": list(RESOLUTIONS),
        "cases": rows,
        "worst_relative_error_at_128": worst,
        "suboff_equivalence": suboff,
        "drilled_solid": drilled,
        "claims": claims,
        "all_claims_hold": all(claims.values()),
    }
    (OUT / "results.json").write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    for claim, passed in claims.items():
        print(("PASS  " if passed else "FAIL  ") + claim)
    print(f"worst relative volume error at 128 cells: {100 * worst:.4f} %")
    print(f"wrote {OUT / 'results.json'}")
    return 0 if payload["all_claims_hold"] else 1


if __name__ == "__main__":
    raise SystemExit(main())