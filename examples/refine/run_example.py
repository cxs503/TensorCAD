"""Measure what ``tensorcad.refine`` claims, and fail the run if it stops holding.

Three questions, in the order a solver engineer asks them:

1. **Is splitting geometry-neutral?** Refine the same mesh to three targets without
   projection and compare volume, area and centroid against the input. The expected
   answer is round-off, not "small". If refinement moves the displacement, the operator
   is not the one this module documents.
2. **Is resolution accuracy?** Take the coarse sphere, refine, project onto the exact
   sphere, and watch the volume error fall. Splitting alone leaves the faceting error
   where it was; the projected column is the one that improves.
3. **Does it work on a hull, against an external number?** Sample the DARPA Suboff
   hull field, refine + project, and compare the displacement volume with the published
   0.699211 m^3. Third-party anchor, not a self-consistency check.

Run: ``python examples/refine/run_example.py --output results/refine``
"""

from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from tensorcad import mesh as meshlib
from tensorcad.geometry_field import SphereField, to_mesh
from tensorcad.refine import refine
from tensorcad.suboff import SuboffHullField

SPHERE_VOLUME = 4.0 / 3.0 * math.pi
SUBOFF_PUBLISHED_VOLUME = 0.699211  # Groves et al. (1989), AFF-1 bare hull


def relative(value: float, reference: float) -> float:
    return value / reference - 1.0


def sphere_sweep(rows: list[dict[str, object]]) -> bool:
    sampled = to_mesh(SphereField(radius=1.0), (18, 18, 18))
    mesh = sampled.mesh
    base_error = relative(meshlib.total_signed_volume(mesh), SPHERE_VOLUME)
    print(f"sphere: {len(mesh.faces)} faces, volume error {100 * base_error:+.4f} %")
    ok = True
    for target in (0.15, 0.10, 0.07):
        _plain, plain_report = refine(mesh, target)
        projected, projection_report = refine(
            mesh, target, project=SphereField(radius=1.0), project_iterations=3
        )
        error = relative(meshlib.total_signed_volume(projected), SPHERE_VOLUME)
        row = {
            "case": "sphere",
            "target": target,
            "faces": len(projected.faces),
            "max_edge": projection_report.after.max_edge,
            "volume_change_plain": plain_report.volume_relative_change,
            "area_change_plain": plain_report.area_relative_change,
            "centroid_shift_plain": plain_report.centroid_shift,
            "volume_error_before": base_error,
            "volume_error_after": error,
            "max_residual_after": (projection_report.projection or {}).get("max_residual_after"),
            "max_displacement": (projection_report.projection or {}).get("max_displacement"),
            "watertight": meshlib.is_watertight(projected),
        }
        rows.append(row)
        print(
            f"  target {target:.2f} m: {row['faces']:6d} faces, max edge {row['max_edge']:.4f} m, "
            f"plain dV {100 * plain_report.volume_relative_change:+.2e} %, "
            f"projected volume error {100 * error:+.4f} % "
            f"(was {100 * base_error:+.4f} %), residual "
            f"{(projection_report.projection or {}).get('max_residual_after', float('nan')):.2e} m"
        )
        ok &= abs(plain_report.volume_relative_change) < 1e-12
        ok &= abs(plain_report.area_relative_change) < 1e-12
        ok &= plain_report.centroid_shift < 1e-12
        ok &= row["watertight"] is True
        ok &= abs(error) <= abs(base_error)
    return bool(ok)


def suboff_case(rows: list[dict[str, object]]) -> bool:
    field = SuboffHullField()
    sampled = to_mesh(field, (28, 6, 8))
    mesh = sampled.mesh
    base_error = relative(meshlib.total_signed_volume(mesh), SUBOFF_PUBLISHED_VOLUME)
    print(
        f"suboff: {len(mesh.faces)} faces, displacement {meshlib.total_signed_volume(mesh):.6f} m^3 "
        f"vs published {SUBOFF_PUBLISHED_VOLUME:.6f} m^3 ({100 * base_error:+.4f} %)"
    )
    ok = True
    for target in (0.12, 0.08):
        projected, report = refine(mesh, target, project=field, project_iterations=3)
        volume = meshlib.total_signed_volume(projected)
        error = relative(volume, SUBOFF_PUBLISHED_VOLUME)
        rows.append(
            {
                "case": "suboff_aff1",
                "target": target,
                "faces": len(projected.faces),
                "max_edge": report.after.max_edge,
                "volume": volume,
                "volume_error_before": base_error,
                "volume_error_after": error,
                "max_residual_after": (report.projection or {}).get("max_residual_after"),
                "max_displacement": (report.projection or {}).get("max_displacement"),
                "min_angle_deg": report.after.min_angle_deg,
                "watertight": meshlib.is_watertight(projected),
            }
        )
        print(
            f"  target {target:.2f} m: {len(projected.faces):6d} faces, displacement "
            f"{volume:.6f} m^3 ({100 * error:+.4f} % vs published)"
        )
        ok &= report.converged and rows[-1]["watertight"] is True
        ok &= abs(error) <= abs(base_error)
    return bool(ok)


def adaptive_case(rows: list[dict[str, object]]) -> bool:
    """A size field that is fine at the bow and coarse aft; compare against uniform."""

    def size_field(point: tuple[float, float, float]) -> float:
        return 0.08 if point[0] < 0.5 * 4.356 else 0.50

    field = SuboffHullField()
    mesh = to_mesh(field, (28, 6, 8)).mesh
    adapted, adapted_report = refine(mesh, size_field)
    uniform, uniform_report = refine(mesh, 0.08)
    rows.append(
        {
            "case": "suboff_adaptive",
            "faces_adaptive": len(adapted.faces),
            "faces_uniform": len(uniform.faces),
            "saving": 1.0 - len(adapted.faces) / len(uniform.faces),
            "adaptive_converged": adapted_report.converged,
            "uniform_converged": uniform_report.converged,
        }
    )
    print(
        f"adaptive: {len(adapted.faces)} faces (h=0.08 forward, 0.50 aft) vs "
        f"{len(uniform.faces)} faces uniform at h=0.08 "
        f"= {100 * (1 - len(adapted.faces) / len(uniform.faces)):.1f} % fewer faces"
    )
    return bool(adapted_report.converged and len(adapted.faces) < len(uniform.faces))


def render_markdown(rows: list[dict[str, object]]) -> str:
    lines = [
        "| case | target [m] | faces | max edge [m] | plain dV | projected error | residual [m] |",
        "|---|---:|---:|---:|---:|---:|---:|",
    ]
    for row in rows:
        if row.get("case") == "suboff_adaptive":
            continue
        lines.append(
            "| {case} | {target:.2f} | {faces} | {max_edge:.4f} | ".format(**row)
            + (
                "{:+.2e} |".format(float(row["volume_change_plain"]))
                if "volume_change_plain" in row
                else "n/a |"
            )
            + " {:+.4f} % | {:+.2e} |".format(
                100 * float(row["volume_error_after"]), float(row["max_residual_after"] or 0.0)
            )
        )
    return "\n".join(lines)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=Path("results/refine"))
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)

    rows: list[dict[str, object]] = []
    ok = sphere_sweep(rows)
    ok &= suboff_case(rows)
    ok &= adaptive_case(rows)
    markdown = render_markdown(rows)
    print("\n" + markdown)
    (args.output / "results.json").write_text(
        json.dumps({"claims_hold": bool(ok), "rows": rows}, indent=2, sort_keys=True) + "\n"
    )
    (args.output / "table.md").write_text(markdown + "\n")
    print(f"\nclaims_hold = {ok} -> {args.output / 'results.json'}")
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
