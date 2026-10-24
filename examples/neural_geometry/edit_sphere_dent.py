"""Revisioned online-fit experiment: add a local dent to a learned sphere.

The baseline and edited geometry are the same ConditionalTSDF architecture.
The baseline weights are retained as an immutable checkpoint. A candidate is
fine-tuned with a local dent target plus replay/preservation constraints, then
accepted only when measured dent and protected-region checks pass.

Run: python examples/neural_geometry/edit_sphere_dent.py --output results/sphere-dent.ngo
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import sys
import time
from pathlib import Path

import torch
from safetensors.torch import save_file
from torch import Tensor

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from train_sphere import ConditionalTSDF, reference_tsdf
from tensorcad.neural_geometry import load_bundle


def dent_tsdf(points: Tensor, radius: Tensor, depth: float, sigma: float, truncation: float) -> Tensor:
    """Reference target for a smooth inward dent centered on the +x pole."""
    norm = torch.linalg.vector_norm(points, dim=1).clamp_min(1e-8)
    cosine = (points[:, 0] / norm).clamp(-1.0, 1.0)
    angle = torch.acos(cosine)
    radial_reduction = depth * torch.exp(-0.5 * (angle / sigma) ** 2)
    return torch.clamp(norm - (radius.reshape(-1) - radial_reduction), -truncation, truncation)


def ray_surface_radii(
    model: ConditionalTSDF, directions: Tensor, radius_value: float, r_min: float, r_max: float
) -> tuple[Tensor, Tensor]:
    radial = torch.linspace(r_min, r_max, 161)
    points = (directions[:, None, :] * radial[None, :, None]).reshape(-1, 3)
    state = torch.full((len(points),), radius_value)
    with torch.no_grad():
        values = model(points, state).reshape(len(directions), len(radial))
    found = torch.zeros(len(directions), dtype=torch.bool)
    roots = torch.full((len(directions),), float("nan"))
    for i, row in enumerate(values):
        indices = torch.nonzero(row[:-1] * row[1:] <= 0, as_tuple=False).flatten()
        if len(indices) == 0:
            continue
        idx = int(indices[0])
        y0, y1 = float(row[idx]), float(row[idx + 1])
        x0, x1 = float(radial[idx]), float(radial[idx + 1])
        frac = 0.5 if abs(y1 - y0) < 1e-12 else -y0 / (y1 - y0)
        roots[i] = x0 + frac * (x1 - x0)
        found[i] = True
    return roots, found


def json_safe_metrics(metrics: dict[str, float]) -> dict[str, float | int | None]:
    """Replace non-finite measurements with null so rollback reports remain valid JSON."""
    safe: dict[str, float | int | None] = {}
    for name, value in metrics.items():
        if isinstance(value, float) and not math.isfinite(value):
            safe[name] = None
        else:
            safe[name] = value
    return safe


def build_bundle(output: Path, model: ConditionalTSDF, radius: float, depth: float, sigma: float,
                 steps: int, seed: int, elapsed: float, metrics: dict[str, float], passed: bool) -> str:
    output.mkdir(parents=True, exist_ok=True)
    weights_path = output / "weights.safetensors"
    save_file({k: v.detach().cpu().contiguous() for k, v in model.state_dict().items()}, str(weights_path))
    digest = hashlib.sha256(weights_path.read_bytes()).hexdigest()
    (output / "state.json").write_text(json.dumps({
        "schema_version": "0.1.0",
        # `field_state` is what NeuralGeometryField conditions on; without it the bundle
        # verifies but cannot be evaluated (`field_from_object` raises a format error).
        "field_state": {"radius": radius},
        "controls": {
            "radius": {"value": radius, "unit": "m", "edit_mode": "conditioned"},
            "dent_depth": {"value": depth, "unit": "m", "edit_mode": "online_fit"},
            "dent_sigma": {"value": sigma, "unit": "rad", "edit_mode": "online_fit"},
        },
    }, indent=2) + "\n", encoding="utf-8")
    (output / "training.json").write_text(json.dumps({
        "schema_version": "0.1.0", "method": "online_fit_with_replay",
        "optimizer": "Adam", "seed": seed, "steps": steps, "learning_rate": 2e-4,
        "parent_revision": 1, "compute": {"device": "cpu", "elapsed_seconds": elapsed},
        "edit": {"type": "smooth_local_inward_dent", "center": "+x", "depth": depth, "sigma_radians": sigma},
    }, indent=2) + "\n", encoding="utf-8")
    (output / "validation.json").write_text(json.dumps({
        "schema_version": "0.1.0", "revision": 2, "weights_sha256": digest,
        "status": "passed" if passed else "failed", "metrics": metrics,
        "tolerances": {"dent_depth_error_m": 0.025, "protected_mean_displacement_m": 0.012,
                       "ray_coverage": 0.98},
    }, indent=2) + "\n", encoding="utf-8")
    (output / "manifest.json").write_text(json.dumps({
        "format": "tensorcad.neural-geometry", "format_version": "0.1.0",
        "object_id": f"urn:tensorcad:geometry:sphere-dent-seed-{seed}", "revision": 2,
        "model": {"family": "conditional_implicit_field", "architecture_id": "tensorcad.mlp_fourier.v1",
                  "architecture_config": {"coordinate_dim": 3, "state_dim": 1, "hidden_width": 128,
                                          "hidden_layers": 4, "coordinate_encoding": "fourier", "output_dim": 1},
                  "weights": {"path": "weights.safetensors", "sha256": digest}},
        "field": {"quantity": "truncated_signed_distance", "inside_sign": "negative",
                  "surface_isovalue": 0.0, "length_unit": "m", "truncation_distance": 0.15},
        "coordinates": {"dimension": 3, "frame": "object_local",
                        "bounds": {"min": [-radius, -radius, -radius],
                                   "max": [radius, radius, radius]}},
        "state_file": "state.json", "training_file": "training.json", "validation_file": "validation.json",
        "state": {"field_state": {"radius": radius}},
        "parent_revision": 1, "producer": {"name": "TensorCAD", "experiment": "online-dent-v0.1"},
    }, indent=2) + "\n", encoding="utf-8")
    return digest


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=Path("results/sphere-dent.ngo"))
    parser.add_argument("--base-steps", type=int, default=1200)
    parser.add_argument("--edit-steps", type=int, default=1000)
    parser.add_argument("--seed", type=int, default=503)
    parser.add_argument("--radius", type=float, default=0.5)
    parser.add_argument("--depth", type=float, default=0.06)
    parser.add_argument("--sigma", type=float, default=0.22)
    args = parser.parse_args()
    if args.base_steps < 1 or args.edit_steps < 1 or args.radius <= 0 or args.depth <= 0 or args.sigma <= 0:
        parser.error("step counts, radius, depth and sigma must be positive")

    torch.manual_seed(args.seed)
    torch.set_num_threads(max(1, min(4, torch.get_num_threads())))
    model = ConditionalTSDF()
    optimizer = torch.optim.Adam(model.parameters(), lr=1e-3)
    truncation = 0.15
    for _ in range(args.base_steps):
        xyz = torch.rand((2048, 3)) * 1.6 - 0.8
        radii = torch.rand(2048) * 0.30 + 0.35
        loss = torch.mean((model(xyz, radii) - reference_tsdf(xyz, radii, truncation)) ** 2)
        optimizer.zero_grad(set_to_none=True)
        loss.backward()
        optimizer.step()

    # Immutable parent snapshot: this is the rollback point and replay teacher.
    parent = ConditionalTSDF()
    parent.load_state_dict({k: v.detach().clone() for k, v in model.state_dict().items()})
    parent.eval()
    candidate = model
    optimizer = torch.optim.Adam(candidate.parameters(), lr=2e-4)
    start = time.perf_counter()
    for _ in range(args.edit_steps):
        # Bias edit samples toward the surface, where the zero level set lives.
        # Half the samples cover the whole sphere; half resolve the target dent
        # cone so a small local edit is not drowned out by the volume loss.
        count = 2048
        general = torch.randn((count // 2, 3))
        general = general / torch.linalg.vector_norm(general, dim=1, keepdim=True).clamp_min(1e-8)
        angle = torch.rand(count - count // 2) * (2.5 * args.sigma)
        azimuth = torch.rand_like(angle) * (2.0 * math.pi)
        dent_directions = torch.stack((
            torch.cos(angle),
            torch.sin(angle) * torch.cos(azimuth),
            torch.sin(angle) * torch.sin(azimuth),
        ), dim=1)
        directions_edit = torch.cat((general, dent_directions), dim=0)
        radial_offset = (torch.rand(count) - 0.5) * 0.30
        xyz = directions_edit * (args.radius + radial_offset)[:, None]
        radii = torch.full((count,), args.radius)
        target = dent_tsdf(xyz, radii, args.depth, args.sigma, truncation)
        candidate_prediction = candidate(xyz, radii)
        target_loss = torch.mean((candidate_prediction - target) ** 2)

        # Replay points from the parent field; explicitly protect the far side
        # of the sphere, while allowing the target edit to change the +x pole.
        directions = torch.randn((2048, 3))
        directions = directions / torch.linalg.vector_norm(directions, dim=1, keepdim=True).clamp_min(1e-8)
        replay_radius = torch.rand(2048) * 0.20 + args.radius - 0.10
        replay_points = directions * replay_radius[:, None]
        replay_state = torch.full((len(replay_points),), args.radius)
        with torch.no_grad():
            teacher = parent(replay_points, replay_state)
        replay_prediction = candidate(replay_points, replay_state)
        cosine = directions[:, 0]
        protected = cosine < math.cos(2.5 * args.sigma)
        if torch.any(protected):
            preservation_loss = torch.mean((replay_prediction[protected] - teacher[protected]) ** 2)
        else:
            preservation_loss = torch.mean((replay_prediction - teacher) ** 2)

        loss = target_loss + 3.0 * preservation_loss
        optimizer.zero_grad(set_to_none=True)
        loss.backward()
        optimizer.step()
    elapsed = time.perf_counter() - start
    candidate.eval()

    generator = torch.Generator().manual_seed(args.seed + 99)
    directions = torch.randn((4096, 3), generator=generator)
    directions = directions / torch.linalg.vector_norm(directions, dim=1, keepdim=True)
    parent_roots, parent_found = ray_surface_radii(parent, directions, args.radius, args.radius - 0.14, args.radius + 0.08)
    candidate_roots, candidate_found = ray_surface_radii(candidate, directions, args.radius, args.radius - 0.14, args.radius + 0.08)
    cosine = directions[:, 0].clamp(-1, 1)
    angles = torch.acos(cosine)
    target_roots = args.radius - args.depth * torch.exp(-0.5 * (angles / args.sigma) ** 2)
    dent_mask = angles < args.sigma
    protected_mask = angles > 2.5 * args.sigma
    common = parent_found & candidate_found
    dent_valid = common & dent_mask
    protected_valid = common & protected_mask
    dent_errors = torch.abs(candidate_roots[dent_valid] - target_roots[dent_valid])
    protected_displacements = torch.abs(candidate_roots[protected_valid] - parent_roots[protected_valid])
    achieved_depths = parent_roots[dent_valid] - candidate_roots[dent_valid]
    metrics = {
        "ray_coverage": float(common.float().mean()),
        "dent_depth_target_m": args.depth,
        "dent_depth_measured_mean_m": float(achieved_depths.mean()) if len(achieved_depths) else float("nan"),
        "dent_surface_error_mean_m": float(dent_errors.mean()) if len(dent_errors) else float("inf"),
        "protected_mean_displacement_m": float(protected_displacements.mean()) if len(protected_displacements) else float("inf"),
        "protected_p95_displacement_m": float(torch.quantile(protected_displacements, 0.95)) if len(protected_displacements) else float("inf"),
        "dent_ray_count": int(dent_valid.sum()),
        "protected_ray_count": int(protected_valid.sum()),
    }
    passed = (
        metrics["ray_coverage"] >= 0.98
        and metrics["dent_surface_error_mean_m"] <= 0.025
        and metrics["protected_mean_displacement_m"] <= 0.012
        and all(math.isfinite(v) for v in metrics.values() if isinstance(v, float))
    )
    if passed:
        build_bundle(args.output, candidate, args.radius, args.depth, args.sigma, args.edit_steps,
                     args.seed, elapsed, metrics, True)
        loaded = load_bundle(args.output)
        assert loaded.revision == 2
        result = {"transaction": "committed", "bundle": str(args.output), "revision": loaded.revision,
                  "weights_sha256": loaded.weights.sha256, **json_safe_metrics(metrics)}
    else:
        # Roll back by retaining the parent model; never publish failed candidate as active.
        result = {"transaction": "rolled_back", "reason": "validation thresholds not met",
                  "active_revision": 1, "candidate_not_published": True, **json_safe_metrics(metrics)}
    print(json.dumps(result, indent=2, allow_nan=False))
    if not passed:
        raise SystemExit("online dent candidate failed validation; parent revision retained")


if __name__ == "__main__":
    main()
