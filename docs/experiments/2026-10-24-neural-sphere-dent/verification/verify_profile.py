#!/usr/bin/env python3
"""Independent angular profile of the baseline and edited neural sphere.

This does *not* reuse the experiment's own ``ray_surface_radii`` helper: it
scans a dense radial grid and bisects the zero crossing, then reports the
surface radius at fixed polar angles measured from the +x dent axis.

It also replays the base training of ``edit_sphere_dent.py`` from its seed
(1200 steps) so the parent field can be compared with the committed candidate
without trusting any value printed by the experiment.

Usage::

    python verification/verify_profile.py --parent-out /tmp/parent.safetensors
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import sys
from pathlib import Path

import torch
from safetensors.torch import load_file, save_file

EXPERIMENT_DIR = Path(__file__).resolve().parents[1]
REPO_ROOT = Path(__file__).resolve().parents[4]
sys.path.insert(0, str(REPO_ROOT / "examples" / "neural_geometry"))
from train_sphere import ConditionalTSDF, reference_tsdf  # noqa: E402

SEED = 503
BASE_STEPS = 1200
RADIUS = 0.5


def rebuild_parent(steps: int, seed: int) -> ConditionalTSDF:
    torch.manual_seed(seed)
    torch.set_num_threads(max(1, min(4, torch.get_num_threads())))
    model = ConditionalTSDF()
    optimizer = torch.optim.Adam(model.parameters(), lr=1e-3)
    for _ in range(steps):
        xyz = torch.rand((2048, 3)) * 1.6 - 0.8
        radii = torch.rand(2048) * 0.30 + 0.35
        loss = torch.mean((model(xyz, radii) - reference_tsdf(xyz, radii, 0.15)) ** 2)
        optimizer.zero_grad(set_to_none=True)
        loss.backward()
        optimizer.step()
    model.eval()
    return model


def surface_radius(model: ConditionalTSDF, angle_deg: float, azimuth_deg: float = 0.0) -> float:
    theta = math.radians(angle_deg)
    phi = math.radians(azimuth_deg)
    direction = torch.tensor([[
        math.cos(theta),
        math.sin(theta) * math.cos(phi),
        math.sin(theta) * math.sin(phi),
    ]])
    grid = torch.linspace(0.20, 0.80, 2401)
    points = direction * grid[:, None]
    state = torch.full((len(points),), RADIUS)
    with torch.no_grad():
        values = model(points, state)
    crossings = torch.nonzero(values[:-1] * values[1:] <= 0, as_tuple=False).flatten()
    if len(crossings) == 0:
        return float("nan")
    index = int(crossings[0])
    lo, hi = float(grid[index]), float(grid[index + 1])
    with torch.no_grad():
        left_value = float(model(direction * lo, torch.full((1,), RADIUS)))
    for _ in range(60):
        mid = 0.5 * (lo + hi)
        with torch.no_grad():
            value = float(model(direction * mid, torch.full((1,), RADIUS)))
        if left_value * value <= 0:
            hi = mid
        else:
            lo = mid
            left_value = value
    return 0.5 * (lo + hi)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--experiment-dir", type=Path, default=EXPERIMENT_DIR)
    parser.add_argument("--parent-out", type=Path, default=Path("/tmp/tc-parent.safetensors"))
    parser.add_argument("--output", type=Path, default=None)
    args = parser.parse_args()
    output = args.output or (args.experiment_dir / "verification" / "angular_profile.json")

    parent = rebuild_parent(BASE_STEPS, SEED)
    tensors = {k: v.detach().cpu().contiguous() for k, v in parent.state_dict().items()}
    save_file(tensors, str(args.parent_out))
    parent_hash = hashlib.sha256(args.parent_out.read_bytes()).hexdigest()

    paths = {
        "baseline_rev1": args.experiment_dir / "bundles/neural-sphere.ngo/weights.safetensors",
        "dent_rev2": args.experiment_dir / "bundles/sphere-dent-seed-503.ngo/weights.safetensors",
    }
    models = {}
    for name, path in paths.items():
        model = ConditionalTSDF()
        model.load_state_dict(load_file(str(path)))
        model.eval()
        models[name] = model

    rows = []
    for angle in (0.0, 5.0, 10.0, 15.0, 20.0, 30.0, 45.0, 90.0, 180.0):
        baseline = surface_radius(models["baseline_rev1"], angle)
        baseline_azimuth = surface_radius(models["baseline_rev1"], angle, 90.0)
        candidate = surface_radius(models["dent_rev2"], angle)
        target = RADIUS - 0.06 * math.exp(-0.5 * (math.radians(angle) / 0.22) ** 2)
        rows.append({
            "angle_deg": angle,
            "baseline_radius_m": baseline,
            "baseline_radius_azimuth90_m": baseline_azimuth,
            "candidate_radius_m": candidate,
            "analytic_dent_radius_m": target,
            "displacement_m": candidate - baseline,
            "candidate_vs_analytic_m": candidate - target,
        })

    report = {
        "method": "dense radial scan (0.20-0.80 m, 2401 samples) plus 60-step bisection; "
                  "parent replayed from seed 503 with 1200 base steps",
        "parent_rebuilt_sha256": parent_hash,
        "candidate_sha256": hashlib.sha256(paths["dent_rev2"].read_bytes()).hexdigest(),
        "baseline_sha256": hashlib.sha256(paths["baseline_rev1"].read_bytes()).hexdigest(),
        "profile": rows,
    }
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(report, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    print(json.dumps(report, indent=2, allow_nan=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
