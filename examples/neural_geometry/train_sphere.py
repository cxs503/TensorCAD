"""Train and validate one state-conditioned neural TSDF on spheres.

Install with `pip install -e ".[neural]"`, then run this script. The analytic
sphere distance is used only to generate labels and validate predictions; all
reported predictions come from the trained neural field.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import time
from pathlib import Path
import sys

import torch
from safetensors.torch import save_file
from torch import Tensor, nn

# Allow direct execution from a source checkout.
ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from tensorcad.neural_geometry import load_bundle  # noqa: E402


class ConditionalTSDF(nn.Module):
    """MLP mapping Fourier-encoded xyz and radius state to a scalar TSDF."""

    def __init__(self, width: int = 128, layers: int = 4, bands: int = 6) -> None:
        super().__init__()
        self.bands = bands
        input_dim = 3 * bands * 2 + 1
        blocks: list[nn.Module] = []
        current = input_dim
        for _ in range(layers):
            blocks.extend((nn.Linear(current, width), nn.SiLU()))
            current = width
        blocks.append(nn.Linear(current, 1))
        self.network = nn.Sequential(*blocks)

    def forward(self, points: Tensor, radius: Tensor) -> Tensor:
        frequencies = (2.0 ** torch.arange(self.bands, device=points.device, dtype=points.dtype)) * math.pi
        phase = points.unsqueeze(-1) * frequencies
        encoded = torch.cat((torch.sin(phase), torch.cos(phase)), dim=-1).flatten(start_dim=1)
        return self.network(torch.cat((encoded, radius.reshape(-1, 1)), dim=1)).squeeze(1)


def reference_tsdf(points: Tensor, radius: Tensor, truncation: float) -> Tensor:
    """Reference oracle used only for supervision and validation."""
    return torch.clamp(torch.linalg.vector_norm(points, dim=1) - radius.reshape(-1), -truncation, truncation)


def evaluate(model: ConditionalTSDF, radius_value: float, seed: int, truncation: float) -> dict[str, float]:
    generator = torch.Generator(device="cpu").manual_seed(seed)
    points = torch.rand((50000, 3), generator=generator) * 1.6 - 0.8
    radius = torch.full((len(points),), radius_value)
    with torch.no_grad():
        predicted = model(points, radius)
        expected = reference_tsdf(points, radius, truncation)
        mae = torch.mean(torch.abs(predicted - expected)).item()
        rmse = torch.sqrt(torch.mean((predicted - expected) ** 2)).item()

        # Find each ray's zero crossing using a fixed radial grid and linear
        # interpolation. This measures the learned surface, not the oracle.
        directions = torch.randn((2048, 3), generator=generator)
        directions = directions / torch.linalg.vector_norm(directions, dim=1, keepdim=True)
        radial_samples = torch.linspace(max(0.01, radius_value - 0.12), radius_value + 0.12, 97)
        samples = (directions[:, None, :] * radial_samples[None, :, None]).reshape(-1, 3)
        radii = torch.full((len(samples),), radius_value)
        values = model(samples, radii).reshape(len(directions), len(radial_samples))
        crossings: list[float] = []
        for row in values:
            indices = torch.nonzero(row[:-1] * row[1:] <= 0, as_tuple=False).flatten()
            if len(indices) == 0:
                crossings.append(float("nan"))
                continue
            index = int(indices[0])
            left_value = float(row[index])
            right_value = float(row[index + 1])
            left_radius = float(radial_samples[index])
            right_radius = float(radial_samples[index + 1])
            denom = right_value - left_value
            fraction = 0.5 if abs(denom) < 1e-12 else -left_value / denom
            crossings.append(left_radius + fraction * (right_radius - left_radius))
        valid = [abs(value - radius_value) for value in crossings if math.isfinite(value)]
        coverage = len(valid) / len(crossings)
        boundary_mae = sum(valid) / len(valid) if valid else float("inf")
        boundary_p95 = sorted(valid)[int(0.95 * (len(valid) - 1))] if valid else float("inf")
        boundary_max = max(valid) if valid else float("inf")
    return {
        "field_mae": mae,
        "field_rmse": rmse,
        "boundary_radial_error_mean": boundary_mae,
        "boundary_radial_error_p95": boundary_p95,
        "boundary_radial_error_max": boundary_max,
        "ray_zero_crossing_coverage": coverage,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=Path("results/neural-sphere.ngo"))
    parser.add_argument("--steps", type=int, default=1800)
    parser.add_argument("--seed", type=int, default=503)
    parser.add_argument("--radius", type=float, default=0.537, help="unseen radius used for validation")
    args = parser.parse_args()
    if args.steps < 1 or not 0.35 <= args.radius <= 0.65:
        parser.error("--steps must be positive and --radius must be within [0.35, 0.65]")

    torch.manual_seed(args.seed)
    torch.set_num_threads(max(1, min(4, torch.get_num_threads())))
    model = ConditionalTSDF()
    optimizer = torch.optim.Adam(model.parameters(), lr=1e-3)
    truncation = 0.15
    start = time.perf_counter()
    final_loss = float("nan")
    model.train()
    for step in range(args.steps):
        points = torch.rand((2048, 3)) * 1.6 - 0.8
        radii = torch.rand(2048) * 0.30 + 0.35
        labels = reference_tsdf(points, radii, truncation)
        prediction = model(points, radii)
        loss = torch.mean((prediction - labels) ** 2)
        optimizer.zero_grad(set_to_none=True)
        loss.backward()
        optimizer.step()
        final_loss = float(loss.detach())
    elapsed = time.perf_counter() - start
    model.eval()
    metrics = evaluate(model, args.radius, args.seed + 1, truncation)

    output = args.output
    output.mkdir(parents=True, exist_ok=True)
    weights_path = output / "weights.safetensors"
    save_file({key: value.detach().cpu().contiguous() for key, value in model.state_dict().items()}, str(weights_path))
    digest = hashlib.sha256(weights_path.read_bytes()).hexdigest()
    (output / "state.json").write_text(json.dumps({
        "schema_version": "0.1.0",
        "controls": {"radius": {"value": args.radius, "unit": "m", "role": "global_parameter", "edit_mode": "conditioned"}},
    }, indent=2) + "\n", encoding="utf-8")
    (output / "training.json").write_text(json.dumps({
        "schema_version": "0.1.0",
        "method": "supervised_tsdf",
        "optimizer": "Adam",
        "seed": args.seed,
        "steps": args.steps,
        "learning_rate": 1e-3,
        "final_batch_mse": final_loss,
        "compute": {"device": "cpu", "elapsed_seconds": elapsed},
    }, indent=2) + "\n", encoding="utf-8")
    tolerances = {"field_mae": 0.008, "boundary_radial_error_mean": 0.01, "ray_zero_crossing_coverage": 0.95}
    passed = (
        metrics["field_mae"] <= tolerances["field_mae"]
        and metrics["boundary_radial_error_mean"] <= tolerances["boundary_radial_error_mean"]
        and metrics["ray_zero_crossing_coverage"] >= tolerances["ray_zero_crossing_coverage"]
    )
    (output / "validation.json").write_text(json.dumps({
        "schema_version": "0.1.0",
        "revision": 1,
        "weights_sha256": digest,
        "status": "passed" if passed else "failed",
        "metrics": metrics,
        "tolerances": tolerances,
        "checks": [
            {"name": "finite_field_values", "passed": all(math.isfinite(value) for value in metrics.values())},
            {"name": "unseen_radius_field_accuracy", "passed": metrics["field_mae"] <= tolerances["field_mae"]},
            {"name": "surface_zero_crossing_coverage", "passed": metrics["ray_zero_crossing_coverage"] >= tolerances["ray_zero_crossing_coverage"]},
        ],
    }, indent=2) + "\n", encoding="utf-8")
    (output / "manifest.json").write_text(json.dumps({
        "format": "tensorcad.neural-geometry",
        "format_version": "0.1.0",
        "object_id": f"urn:tensorcad:geometry:sphere-seed-{args.seed}",
        "revision": 1,
        "model": {
            "family": "conditional_implicit_field",
            "architecture_id": "tensorcad.mlp_fourier.v1",
            "architecture_config": {
                "coordinate_dim": 3, "state_dim": 1, "hidden_width": 128,
                "hidden_layers": 4, "coordinate_encoding": "fourier", "output_dim": 1,
            },
            "weights": {"path": "weights.safetensors", "sha256": digest},
        },
        "field": {
            "quantity": "truncated_signed_distance", "inside_sign": "negative",
            "surface_isovalue": 0.0, "length_unit": "m", "truncation_distance": truncation,
        },
        "coordinates": {"dimension": 3, "frame": "object_local"},
        "state_file": "state.json", "training_file": "training.json", "validation_file": "validation.json",
        "producer": {"name": "TensorCAD", "experiment": "neural-sphere-v0.1"},
    }, indent=2) + "\n", encoding="utf-8")

    # Round-trip the metadata and hash check through the production bundle loader.
    loaded = load_bundle(output)
    report = {
        "bundle": str(output),
        "revision": loaded.revision,
        "architecture_id": loaded.architecture_id,
        "steps": args.steps,
        "elapsed_seconds": elapsed,
        "validation_status": "passed" if passed else "failed",
        **metrics,
    }
    print(json.dumps(report, indent=2))
    if not passed:
        raise SystemExit("neural sphere validation thresholds were not met; see validation.json")


if __name__ == "__main__":
    main()
