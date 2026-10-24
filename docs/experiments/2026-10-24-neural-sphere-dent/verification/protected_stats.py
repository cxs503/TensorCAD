#!/usr/bin/env python3
"""Independent protected-region statistics for the committed dent candidate.

Samples 4096 directions on a fixed seed, locates the outermost zero crossing on
a dense radial grid, and reports mean/p95/p99/max displacement of the candidate
against the replayed parent (see ``verify_profile.py``). Also reports the
baseline field's own radial form error, so displacement caused by the edit is
not confused with the model's intrinsic error.

Usage::

    python verification/verify_profile.py --parent-out /tmp/parent.safetensors
    python verification/protected_stats.py --parent /tmp/parent.safetensors
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import torch
from safetensors.torch import load_file

EXPERIMENT_DIR = Path(__file__).resolve().parents[1]
REPO_ROOT = Path(__file__).resolve().parents[4]
sys.path.insert(0, str(REPO_ROOT / "examples" / "neural_geometry"))
from train_sphere import ConditionalTSDF  # noqa: E402

RADIUS = 0.5
SIGMA = 0.22


def roots(model: ConditionalTSDF, directions: torch.Tensor, grid: torch.Tensor) -> torch.Tensor:
    out = torch.full((len(directions),), float("nan"))
    points = (directions[:, None, :] * grid[None, :, None]).reshape(-1, 3)
    values = torch.empty(len(points))
    with torch.no_grad():
        for start in range(0, len(points), 200_000):
            chunk = points[start:start + 200_000]
            values[start:start + 200_000] = model(chunk, torch.full((len(chunk),), RADIUS))
    values = values.reshape(len(directions), len(grid))
    for i in range(len(directions)):
        row = values[i]
        hits = torch.nonzero(row[:-1] * row[1:] <= 0, as_tuple=False).flatten()
        if len(hits) == 0:
            continue
        index = int(hits[0])
        y0, y1 = float(row[index]), float(row[index + 1])
        x0, x1 = float(grid[index]), float(grid[index + 1])
        frac = 0.5 if abs(y1 - y0) < 1e-12 else -y0 / (y1 - y0)
        out[i] = x0 + frac * (x1 - x0)
    return out


def stats(mask: torch.Tensor, values: torch.Tensor) -> dict[str, float | int]:
    selected = values[mask]
    if len(selected) == 0:
        return {"count": 0}
    absolute = selected.abs()
    quantiles = torch.quantile(absolute, torch.tensor([0.5, 0.95, 0.99]), interpolation="linear")
    return {
        "count": len(selected),
        "mean_signed_m": float(selected.mean()),
        "mean_abs_m": float(absolute.mean()),
        "median_abs_m": float(quantiles[0]),
        "p95_abs_m": float(quantiles[1]),
        "p99_abs_m": float(quantiles[2]),
        "max_abs_m": float(absolute.max()),
        "outward_fraction": float((selected > 0).float().mean()),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--experiment-dir", type=Path, default=EXPERIMENT_DIR)
    parser.add_argument("--parent", type=Path, required=True,
                        help="parent weights written by verify_profile.py --parent-out")
    parser.add_argument("--candidate", type=Path, default=None)
    parser.add_argument("--output", type=Path, default=None)
    args = parser.parse_args()
    output = args.output or (args.experiment_dir / "verification" / "protected_region_stats.json")
    candidate_path = args.candidate or (
        args.experiment_dir / "bundles/sphere-dent-seed-503.ngo/weights.safetensors"
    )

    torch.set_num_threads(max(1, min(4, torch.get_num_threads())))
    parent = ConditionalTSDF()
    parent.load_state_dict(load_file(str(args.parent)))
    candidate = ConditionalTSDF()
    candidate.load_state_dict(load_file(str(candidate_path)))
    parent.eval()
    candidate.eval()

    grid = torch.linspace(0.20, 0.80, 1201)
    generator = torch.Generator().manual_seed(12345)
    directions = torch.randn((4096, 3), generator=generator)
    directions = directions / torch.linalg.vector_norm(directions, dim=1, keepdim=True)

    parent_roots = roots(parent, directions, grid)
    candidate_roots = roots(candidate, directions, grid)
    valid = torch.isfinite(parent_roots) & torch.isfinite(candidate_roots)
    angles = torch.acos(directions[:, 0].clamp(-1.0, 1.0))
    displacement = candidate_roots - parent_roots
    baseline_form_error = parent_roots - RADIUS
    analytic_target = RADIUS - 0.06 * torch.exp(-0.5 * (angles / SIGMA) ** 2)

    dent = valid & (angles < SIGMA)
    protected = valid & (angles > 2.5 * SIGMA)
    transition = valid & (angles >= SIGMA) & (angles <= 2.5 * SIGMA)

    worst_index = displacement[protected].abs().argmax()
    report = {
        "method": "4096 fixed-seed directions, first zero crossing on a 1201-point radial grid",
        "directions": len(directions),
        "coverage": float(valid.float().mean()),
        "baseline_form_error_m": {
            "mean_abs": float(baseline_form_error[valid].abs().mean()),
            "p95_abs": float(torch.quantile(baseline_form_error[valid].abs(), 0.95)),
            "max_abs": float(baseline_form_error[valid].abs().max()),
            "signed_mean": float(baseline_form_error[valid].mean()),
        },
        "dent_region": stats(dent, displacement),
        "transition_region": stats(transition, displacement),
        "protected_region": stats(protected, displacement),
        "protected_vs_analytic_surface": stats(protected, candidate_roots - analytic_target),
        "worst_protected": {
            "angle_deg": float(torch.rad2deg(angles[protected][worst_index])),
            "displacement_m": float(displacement[protected][worst_index]),
        },
    }
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(report, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
