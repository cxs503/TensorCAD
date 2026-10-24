"""The learned field: one state-conditioned implicit architecture for every object.

The field maps a point ``x`` and an editable state ``c`` to a truncated signed
distance: negative inside, positive outside, ``|value| <= truncation`` near the
surface. Nothing here is shape-specific — a sphere, an ellipsoid and a dented
sphere are the same network with different weights and states, which is the
whole point: the geometry *is* the parameters plus the field, and a mesh is only
what you get when you ask for one.

Reference signed-distance functions in this module are *oracles*: they generate
training labels and expected values for tests. They are never used at query time
to answer what the geometry is.
"""
from __future__ import annotations

import math
import time
from collections.abc import Callable, Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import Any, cast

import numpy as np
import numpy.typing as npt
import torch
from safetensors.torch import load_file as load_safetensors
from safetensors.torch import save_file as save_safetensors
from torch import Tensor, nn

from ..neural_geometry import sha256_file

#: Same parameterisation as the extraction module, for the same reason.
Array = npt.NDArray[Any]

ARCHITECTURE_ID = "tensorcad.mlp_fourier.v1"
DEFAULT_TRUNCATION = 0.15


def _to_tensor(array: Array) -> Tensor:
    """Convert numpy -> torch without assuming the optional NumPy bridge works.

    ``torch`` wheels are compiled against a specific NumPy ABI; a mismatched pair
    (for example torch 2.2 with NumPy 2.x) still runs tensor maths but raises
    "Numpy is not available" from ``torch.from_numpy``. Falling back to a list
    conversion keeps the object usable instead of making an environment detail
    look like a geometry failure.
    """
    try:
        return torch.from_numpy(np.ascontiguousarray(array, dtype=np.float64)).float()
    except (RuntimeError, TypeError):  # pragma: no cover - depends on the installed wheel
        return torch.tensor(array.tolist(), dtype=torch.float32)


def _to_numpy(tensor: Tensor) -> Array:
    try:
        return np.asarray(tensor.detach().double().numpy(), dtype=np.float64)
    except (RuntimeError, TypeError):  # pragma: no cover - depends on the installed wheel
        return np.array(tensor.detach().double().tolist(), dtype=np.float64)


@dataclass(frozen=True)
class FieldConfig:
    """Architecture hyper-parameters, persisted verbatim in the manifest."""

    width: int = 128
    layers: int = 4
    bands: int = 6
    state_dim: int = 1

    def to_dict(self) -> dict[str, Any]:
        return {
            "coordinate_dim": 3,
            "state_dim": self.state_dim,
            "hidden_width": self.width,
            "hidden_layers": self.layers,
            "coordinate_encoding": "fourier",
            "fourier_bands": self.bands,
            "output_dim": 1,
        }

    @classmethod
    def from_dict(cls, data: dict[str, Any]) -> FieldConfig:
        return cls(
            width=int(data.get("hidden_width", 128)),
            layers=int(data.get("hidden_layers", 4)),
            bands=int(data.get("fourier_bands", 6)),
            state_dim=int(data.get("state_dim", 1)),
        )


class ConditionalTSDF(nn.Module):  # type: ignore[misc]  # nn is untyped by design
    """MLP over Fourier-encoded coordinates plus the state vector, one scalar out."""

    def __init__(self, config: FieldConfig | None = None) -> None:
        super().__init__()
        self.config = config or FieldConfig()
        input_dim = 3 * self.config.bands * 2 + self.config.state_dim
        blocks: list[nn.Module] = []
        current = input_dim
        for _ in range(self.config.layers):
            blocks.extend((nn.Linear(current, self.config.width), nn.SiLU()))
            current = self.config.width
        blocks.append(nn.Linear(current, 1))
        self.network = nn.Sequential(*blocks)

    def forward(self, points: Tensor, state: Tensor) -> Tensor:
        frequencies = (
            2.0 ** torch.arange(self.config.bands, device=points.device, dtype=points.dtype)
        ) * math.pi
        phase = points.unsqueeze(-1) * frequencies
        encoded = torch.cat((torch.sin(phase), torch.cos(phase)), dim=-1).flatten(start_dim=1)
        conditioned = torch.cat((encoded, state.reshape(len(points), -1)), dim=1)
        return cast(Tensor, self.network(conditioned).squeeze(1))


# --------------------------------------------------------------------------- #
# reference oracles (labels and expected values only)
# --------------------------------------------------------------------------- #


def sphere_tsdf(points: Tensor, radius: Tensor, truncation: float) -> Tensor:
    return torch.clamp(
        torch.linalg.vector_norm(points, dim=1) - radius.reshape(-1), -truncation, truncation
    )


def dented_sphere_tsdf(
    points: Tensor, radius: Tensor, depth: float, sigma: float, truncation: float
) -> Tensor:
    """Sphere with a smooth Gaussian inward dent centred on the +x pole."""
    norm = torch.linalg.vector_norm(points, dim=1).clamp_min(1e-8)
    cosine = (points[:, 0] / norm).clamp(-1.0, 1.0)
    angle = torch.acos(cosine)
    reduction = depth * torch.exp(-0.5 * (angle / sigma) ** 2)
    return torch.clamp(norm - (radius.reshape(-1) - reduction), -truncation, truncation)


def dent_target_radius(
    radius: float, angle: Array, depth: float, sigma: float
) -> Array:
    """Analytic radius of the dented sphere along ``angle`` from the dent axis."""
    return radius - depth * np.exp(-0.5 * (angle / sigma) ** 2)


# --------------------------------------------------------------------------- #
# training
# --------------------------------------------------------------------------- #


@dataclass(frozen=True)
class TrainOutcome:
    model: ConditionalTSDF
    elapsed_seconds: float
    final_loss: float
    steps: int
    seed: int


def set_threads(limit: int = 4) -> None:
    torch.set_num_threads(max(1, min(limit, torch.get_num_threads())))


def train_base_field(
    *,
    steps: int,
    seed: int,
    radius_range: tuple[float, float] = (0.35, 0.65),
    batch: int = 2048,
    learning_rate: float = 1e-3,
    truncation: float = DEFAULT_TRUNCATION,
    config: FieldConfig | None = None,
    log_every: int = 0,
    log: Callable[[str], None] | None = None,
) -> TrainOutcome:
    """Fit the conditional field to the analytic sphere family over a radius range."""
    if steps < 1:
        raise ValueError("steps must be positive")
    torch.manual_seed(seed)
    set_threads()
    model = ConditionalTSDF(config)
    optimizer = torch.optim.Adam(model.parameters(), lr=learning_rate)
    low, high = radius_range
    start = time.perf_counter()
    final_loss = float("nan")
    model.train()
    for step in range(steps):
        points = torch.rand((batch, 3)) * 1.6 - 0.8
        radii = torch.rand(batch) * (high - low) + low
        labels = sphere_tsdf(points, radii, truncation)
        loss = torch.mean((model(points, radii) - labels) ** 2)
        optimizer.zero_grad(set_to_none=True)
        loss.backward()
        optimizer.step()
        final_loss = float(loss.detach())
        if log is not None and log_every and step % log_every == 0:
            log(f"step {step:5d}  loss {final_loss:.6f}")
    model.eval()
    return TrainOutcome(
        model=model,
        elapsed_seconds=time.perf_counter() - start,
        final_loss=final_loss,
        steps=steps,
        seed=seed,
    )


def fine_tune_dent(
    parent: ConditionalTSDF,
    *,
    radius: float,
    depth: float,
    sigma: float,
    steps: int,
    seed: int,
    learning_rate: float = 2e-4,
    truncation: float = DEFAULT_TRUNCATION,
    batch: int = 2048,
    replay_weight: float = 3.0,
) -> TrainOutcome:
    """Fit a local dent while replaying the parent field outside the edit cone.

    Half of each batch lands inside the dent cone so a small local edit is not
    drowned out by the volume term; the replay term keeps the parent's answer on
    the rest of the surface, which is what makes the edit *local*.
    """
    if steps < 1:
        raise ValueError("steps must be positive")
    torch.manual_seed(seed)
    set_threads()
    candidate = ConditionalTSDF(parent.config)
    candidate.load_state_dict({key: value.detach().clone() for key, value in parent.state_dict().items()})
    candidate.train()
    optimizer = torch.optim.Adam(candidate.parameters(), lr=learning_rate)
    start = time.perf_counter()
    final_loss = float("nan")
    for _ in range(steps):
        half = batch // 2
        general = torch.randn((half, 3))
        general = general / torch.linalg.vector_norm(general, dim=1, keepdim=True).clamp_min(1e-8)
        angle = torch.rand(batch - half) * (2.5 * sigma)
        azimuth = torch.rand_like(angle) * (2.0 * math.pi)
        dent_directions = torch.stack(
            (
                torch.cos(angle),
                torch.sin(angle) * torch.cos(azimuth),
                torch.sin(angle) * torch.sin(azimuth),
            ),
            dim=1,
        )
        directions = torch.cat((general, dent_directions), dim=0)
        radial = radius + (torch.rand(batch) - 0.5) * 0.30
        points = directions * radial[:, None]
        state = torch.full((batch,), radius)
        target = dented_sphere_tsdf(points, state, depth, sigma, truncation)
        target_loss = torch.mean((candidate(points, state) - target) ** 2)

        replay_directions = torch.randn((batch, 3))
        replay_directions = replay_directions / torch.linalg.vector_norm(
            replay_directions, dim=1, keepdim=True
        ).clamp_min(1e-8)
        replay_radius = torch.rand(batch) * 0.20 + radius - 0.10
        replay_points = replay_directions * replay_radius[:, None]
        with torch.no_grad():
            teacher = parent(replay_points, state)
        prediction = candidate(replay_points, state)
        protected = replay_directions[:, 0] < math.cos(2.5 * sigma)
        if torch.any(protected):
            preservation = torch.mean((prediction[protected] - teacher[protected]) ** 2)
        else:
            preservation = torch.mean((prediction - teacher) ** 2)

        loss = target_loss + replay_weight * preservation
        optimizer.zero_grad(set_to_none=True)
        loss.backward()
        optimizer.step()
        final_loss = float(loss.detach())
    candidate.eval()
    return TrainOutcome(
        model=candidate,
        elapsed_seconds=time.perf_counter() - start,
        final_loss=final_loss,
        steps=steps,
        seed=seed,
    )


# --------------------------------------------------------------------------- #
# weights and querying
# --------------------------------------------------------------------------- #


def save_weights(model: ConditionalTSDF, path: str | Path) -> tuple[str, int]:
    """Write safetensors weights; returns ``(sha256, size_bytes)``."""
    target = Path(path)
    target.parent.mkdir(parents=True, exist_ok=True)
    tensors = {key: value.detach().cpu().contiguous() for key, value in model.state_dict().items()}
    save_safetensors(tensors, str(target))
    return sha256_file(target), target.stat().st_size


def load_weights(path: str | Path, config: FieldConfig | None = None) -> ConditionalTSDF:
    model = ConditionalTSDF(config)
    model.load_state_dict(load_safetensors(str(path)))
    model.eval()
    return model


def field_callable(
    model: ConditionalTSDF, state_value: float
) -> Callable[[Array], Array]:
    """Wrap the torch model as a numpy ``points -> values`` callable."""

    def evaluate(points: Array) -> Array:
        tensor = _to_tensor(points)
        state = torch.full((len(tensor),), float(state_value), dtype=torch.float32)
        with torch.no_grad():
            values = model(tensor, state)
        return _to_numpy(cast(Tensor, values))

    return evaluate


def gradient_callable(
    model: ConditionalTSDF, state_value: float
) -> Callable[[Array], Array]:
    """Wrap the model's exact gradient (torch autograd) as a numpy callable.

    The field is a float32 network, so a central difference of its output is only
    accurate for steps around ``1e-4`` of the object's size; below that the
    difference of two nearly equal float32 values is noise. Autograd has no such
    problem, which is why the protocol's numerical fallback is not the default
    here.
    """

    def evaluate(points: Array) -> Array:
        tensor = _to_tensor(points).requires_grad_(True)
        state = torch.full((len(tensor),), float(state_value), dtype=torch.float32)
        with torch.enable_grad():
            values = model(tensor, state)
            (gradient,) = torch.autograd.grad(cast(Tensor, values).sum(), tensor)
        return _to_numpy(cast(Tensor, gradient))

    return evaluate


def _sample_directions(count: int, generator: torch.Generator) -> Tensor:
    directions = torch.randn((count, 3), generator=generator)
    return cast(
        Tensor,
        directions / torch.linalg.vector_norm(directions, dim=1, keepdim=True).clamp_min(1e-12),
    )


def surface_radii(
    model: ConditionalTSDF,
    directions: Tensor,
    state_value: float,
    radial_range: tuple[float, float],
    *,
    samples: int = 161,
) -> tuple[Tensor, Tensor]:
    """Radius of the first zero crossing along each direction (NaN where none)."""
    radial = torch.linspace(radial_range[0], radial_range[1], samples, dtype=torch.float32)
    points = (directions[:, None, :] * radial[None, :, None]).reshape(-1, 3)
    state = torch.full((len(points),), float(state_value), dtype=torch.float32)
    values = torch.empty(len(points))
    with torch.no_grad():
        for start in range(0, len(points), 200_000):
            block = points[start : start + 200_000]
            values[start : start + 200_000] = model(block, state[start : start + 200_000])
    values = values.reshape(len(directions), samples)

    sign_change = values[:, :-1] * values[:, 1:] <= 0
    found = sign_change.any(dim=1)
    index = torch.where(found, sign_change.float().argmax(dim=1), torch.zeros_like(found, dtype=torch.long))
    left_value = values.gather(1, index[:, None]).squeeze(1)
    right_value = values.gather(1, (index + 1)[:, None]).squeeze(1)
    denominator = right_value - left_value
    fraction = torch.where(
        denominator.abs() < 1e-12,
        torch.full_like(denominator, 0.5),
        -left_value / torch.where(denominator.abs() < 1e-12, torch.ones_like(denominator), denominator),
    )
    left_radius = radial[index]
    right_radius = radial[index + 1]
    roots = left_radius + fraction * (right_radius - left_radius)
    return torch.where(found, roots, torch.full_like(roots, float("nan"))), found


@dataclass(frozen=True)
class FieldMetrics:
    """Measured behaviour of a field, used for both reporting and acceptance."""

    field_mae: float
    field_rmse: float
    boundary_error_mean: float
    boundary_error_p95: float
    boundary_error_max: float
    ray_coverage: float
    rays: int

    def to_dict(self) -> dict[str, float | int]:
        return {
            "field_mae_m": self.field_mae,
            "field_rmse_m": self.field_rmse,
            "boundary_error_mean_m": self.boundary_error_mean,
            "boundary_error_p95_m": self.boundary_error_p95,
            "boundary_error_max_m": self.boundary_error_max,
            "ray_coverage": self.ray_coverage,
            "ray_count": self.rays,
        }


def measure_sphere(
    model: ConditionalTSDF,
    radius: float,
    *,
    seed: int,
    truncation: float = DEFAULT_TRUNCATION,
    points: int = 50_000,
    rays: int = 2048,
) -> FieldMetrics:
    """Compare the field against the analytic sphere it was fitted to."""
    generator = torch.Generator().manual_seed(seed)
    sample_points = torch.rand((points, 3), generator=generator) * 1.6 - 0.8
    state = torch.full((points,), radius)
    with torch.no_grad():
        predicted = model(sample_points, state)
        expected = sphere_tsdf(sample_points, state, truncation)
        mae = torch.mean(torch.abs(predicted - expected)).item()
        rmse = torch.sqrt(torch.mean((predicted - expected) ** 2)).item()

    directions = _sample_directions(rays, generator)
    roots, found = surface_radii(model, directions, radius, (max(0.01, radius - 0.20), radius + 0.20))
    valid = found & torch.isfinite(roots)
    errors = (roots[valid] - radius).abs()
    coverage = float(valid.float().mean())
    if len(errors) == 0:
        return FieldMetrics(mae, rmse, float("inf"), float("inf"), float("inf"), coverage, int(rays))
    return FieldMetrics(
        field_mae=float(mae),
        field_rmse=float(rmse),
        boundary_error_mean=float(errors.mean()),
        boundary_error_p95=float(torch.quantile(errors, 0.95)),
        boundary_error_max=float(errors.max()),
        ray_coverage=coverage,
        rays=int(rays),
    )


def _dent_cone_directions(sigma: float, *, per_ring: int = 16, rings: int = 4) -> Tensor:
    """Fixed polar grid inside the dent cone, expressed about the +x axis."""
    fractions = (np.arange(rings, dtype=np.float64) + 0.5) / rings
    angles = fractions * min(sigma, math.pi) * 0.98
    azimuths = np.arange(per_ring, dtype=np.float64) * (2.0 * math.pi / per_ring)
    theta_grid, phi_grid = np.meshgrid(angles, azimuths, indexing="ij")
    theta = theta_grid.reshape(-1)
    phi = phi_grid.reshape(-1)
    directions = np.stack(
        (np.cos(theta), np.sin(theta) * np.cos(phi), np.sin(theta) * np.sin(phi)), axis=1
    )
    return torch.tensor(directions, dtype=torch.float32)


def measure_dent_edit(
    parent: ConditionalTSDF,
    candidate: ConditionalTSDF,
    *,
    radius: float,
    depth: float,
    sigma: float,
    seed: int,
    rays: int = 4096,
    margin: float = 0.14,
    band: float = 0.08,
) -> dict[str, float]:
    """Measure a dent edit in the terms the acceptance gates actually use.

    Reported separately: the peak depth on the dent axis, the error against the
    analytic dented surface inside the dent cone, and the displacement of the
    protected region measured *relative to the parent*, so the field's own
    approximation error is not counted as edit drift.
    """
    generator = torch.Generator().manual_seed(seed)
    # The dent cone is narrow, so random rays alone may contain none of it and
    # the dent metrics would silently become infinite when the ray count is
    # small. Always add a deterministic stratified cone so the measured region
    # is the region that was edited, independent of the sampling seed.
    cone = _dent_cone_directions(sigma, per_ring=16, rings=4)
    directions = torch.cat((cone, _sample_directions(rays, generator)), dim=0)
    span = (max(0.01, radius - margin), radius + band)
    parent_roots, parent_found = surface_radii(parent, directions, radius, span, samples=241)
    candidate_roots, candidate_found = surface_radii(candidate, directions, radius, span, samples=241)
    common = parent_found & candidate_found

    angles = torch.acos(directions[:, 0].clamp(-1.0, 1.0))
    target = torch.tensor(
        dent_target_radius(radius, _to_numpy(angles), depth, sigma), dtype=torch.float32
    )
    dent_mask = common & (angles < sigma)
    protected_mask = common & (angles > 2.5 * sigma)

    def safe_mean(values: Tensor) -> float:
        return float(values.mean()) if len(values) else float("inf")

    def safe_quantile(values: Tensor, q: float) -> float:
        return float(torch.quantile(values, q)) if len(values) else float("inf")

    dent_errors = (candidate_roots[dent_mask] - target[dent_mask]).abs()
    displacement = candidate_roots[protected_mask] - parent_roots[protected_mask]
    axis = torch.tensor([[1.0, 0.0, 0.0]])
    axis_parent, axis_parent_found = surface_radii(parent, axis, radius, span, samples=241)
    axis_candidate, axis_candidate_found = surface_radii(candidate, axis, radius, span, samples=241)
    peak = (
        float(axis_parent[0] - axis_candidate[0])
        if bool(axis_parent_found[0] and axis_candidate_found[0])
        else float("nan")
    )

    return {
        "ray_coverage": float(common.float().mean()),
        "ray_count": int(rays),
        "dent_peak_target_m": depth,
        "dent_peak_measured_m": peak,
        "dent_peak_error_m": abs(peak - depth) if math.isfinite(peak) else float("inf"),
        "dent_region_error_mean_m": safe_mean(dent_errors),
        "dent_region_ray_count": int(dent_mask.sum()),
        "protected_mean_displacement_m": float(displacement.mean()) if len(displacement) else float("inf"),
        "protected_mean_abs_displacement_m": safe_mean(displacement.abs()),
        "protected_p95_abs_displacement_m": safe_quantile(displacement.abs(), 0.95),
        "protected_p99_abs_displacement_m": safe_quantile(displacement.abs(), 0.99),
        "protected_max_abs_displacement_m": (
            float(displacement.abs().max()) if len(displacement) else float("inf")
        ),
        "protected_ray_count": int(protected_mask.sum()),
    }


def environment() -> dict[str, Any]:
    """Recorded with every revision so a run can be attributed to a stack."""
    return {
        "python": ".".join(str(part) for part in __import__("sys").version_info[:3]),
        "torch": torch.__version__,
        "numpy": np.__version__,
        "device": "cuda" if torch.cuda.is_available() else "cpu",
        "threads": torch.get_num_threads(),
    }


def normalise_direction(direction: Sequence[float]) -> tuple[float, float, float]:
    vector = np.asarray(direction, dtype=np.float64).reshape(3)
    norm = float(np.linalg.norm(vector))
    if norm <= 0.0 or not math.isfinite(norm):
        raise ValueError("edit direction must be a finite non-zero vector")
    unit = vector / norm
    return (float(unit[0]), float(unit[1]), float(unit[2]))
