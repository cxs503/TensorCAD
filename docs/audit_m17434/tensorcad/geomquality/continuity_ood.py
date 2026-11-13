"""Out-of-distribution and resolution audit of the calibrated G0/G1/G2 classifier.

:mod:`.continuity` calibrates its two limits on ``label_seam_dataset`` -- cubic plate
patches joined at ``u = 1``, planar, straight seam, one defect kind per case -- and
scores them on a held-out split of the *same* generator.  That establishes interpolation
within one family; it says nothing about the two things a shipyard model contains and the
generator does not:

**Regime.**  The curvature criterion is a *relative* jump, but the ratio is taken against
``max(|kappa_a|, |kappa_b|, 1/scale)``, so below a curvature of about ``1/scale`` the
criterion is effectively *absolute*: the ``1/scale`` floor dominates and a large
relative mismatch is reported clean.  This module maps that transition instead of
assuming it away -- :func:`curvature_scale_sweep` walks the absolute curvature level from
``1e-4`` to ``1e1`` per unit length at fixed relative mismatch and locates where the
documented "curvatures off by 30 %" reading stops holding.

**Resolution.**  :func:`~.continuity.seam_metrics` samples each seam isocurve at
``samples = 17`` points and takes the maximum.  A defect narrower than one sample
interval is invisible to that maximum, no matter how large it is, so a local ripple
along the seam can pass a "G2" verdict.  :func:`localized_curvature_defect` builds such a
ripple analytically (a curvature discontinuity whose magnitude is a bump function of the
position *along* the seam), :func:`sampling_alias_sweep` measures it at 17 samples and at
a dense reference, and :func:`~.continuity.seam_metrics_adaptive` is the resolution fix
the sweep is used to justify.

The classifier itself is not touched: every number here is measured from the same
analytic surfaces at two resolutions, so the dense measurement is the ground truth and
the coarse one is the production default.
"""
from __future__ import annotations

import math
from typing import Any

from ..nurbs import NURBSSurface
from .brep import Face
from .continuity import (
    DEFAULT_ANGLE_LIMIT_DEG,
    DEFAULT_CURVATURE_LIMIT,
    classify,
    seam_metrics,
    seam_metrics_adaptive,
)

Vec3 = tuple[float, float, float]

#: Resolution used as the "truth" in the aliasing study: 2049 samples across a seam whose
#: defect is at most 1/8 of the seam long, i.e. hundreds of samples across the defect.
DENSE_SAMPLES = 2049

#: Production default, quoted so the sweeps cannot drift away from the code they audit.
PRODUCTION_SAMPLES = 17


# --------------------------------------------------------------------------- #
# seam constructors (independent of ``continuity.make_seam_case``)
# --------------------------------------------------------------------------- #
def _clamped_uniform_knots(spans: int, degree: int) -> tuple[float, ...]:
    """Clamped knot vector with ``spans`` equal spans on ``[0, 1]``."""
    if spans < 1:
        raise ValueError("spans must be >= 1")
    if degree < 1:
        raise ValueError("degree must be >= 1")
    interior = [index / spans for index in range(1, spans)]
    return tuple([0.0] * (degree + 1) + interior + [1.0] * (degree + 1))


def _clustered_knots(
    spans: int,
    degree: int,
    centre: float,
    half_width: float,
    cluster: int,
) -> tuple[float, ...]:
    """Clamped uniform knots plus ``cluster`` extra spans inside ``centre +- half_width``.

    The extra spans put control points where the defect is, so the *geometry* resolves a
    narrow bump and the sampler resolution is the only thing left that can miss it.
    """
    low = max(0.0, centre - half_width)
    high = min(1.0, centre + half_width)
    if not low < high:
        raise ValueError("cluster window is empty")
    interior = [index / spans for index in range(1, spans)]
    interior.extend(low + (high - low) * index / (cluster + 1) for index in range(1, cluster + 1))
    interior = sorted(set(interior))
    return tuple([0.0] * (degree + 1) + interior + [1.0] * (degree + 1))

def _greville(knots: tuple[float, ...], degree: int) -> list[float]:
    count = len(knots) - degree - 1
    return [sum(knots[i + 1 : i + degree + 1]) / degree for i in range(count)]


def _face(surface: NURBSSurface, index: int) -> Face:
    return Face(index=index, surface=surface, boundary=())


def curved_plate(
    curvature: float,
    length: float,
    breadth: float,
    *,
    incoming: bool = True,
) -> NURBSSurface:
    """A quadratic plate whose across-seam curvature at the seam is ``curvature``.

    ``incoming`` places the plate on ``x in [-length, 0]`` and the shape is
    ``z = kappa x^2 / 2``; the outgoing plate is the same profile on ``x in [0, length]``.
    Both have a horizontal tangent at ``x = 0``, so a pair with different ``curvature``
    is a G1 seam by construction whose curvature jump is exactly the difference of the
    two arguments.
    """
    half = 0.5 * curvature * length * length
    if incoming:
        xs = (-length, -0.5 * length, 0.0)
        zs = (half, 0.0, 0.0)
    else:
        xs = (0.0, 0.5 * length, length)
        zs = (0.0, 0.0, half)
    net = [[(x, 0.0, z), (x, breadth, z)] for x, z in zip(xs, zs)]
    return NURBSSurface(
        degree=(2, 1),
        control_net=net,
        weights=[[1.0, 1.0] for _ in net],
        knots=((0.0, 0.0, 0.0, 1.0, 1.0, 1.0), (0.0, 0.0, 1.0, 1.0)),
    )


def straight_plate(length: float, breadth: float, *, incoming: bool = True) -> NURBSSurface:
    """Bilinear plate at ``z = 0``: the ``curvature = 0`` member of the family above."""
    xs = (-length, 0.0) if incoming else (0.0, length)
    net = [[(x, 0.0, 0.0), (x, breadth, 0.0)] for x in xs]
    return NURBSSurface(
        degree=(1, 1),
        control_net=net,
        weights=[[1.0, 1.0] for _ in net],
        knots=((0.0, 0.0, 1.0, 1.0), (0.0, 0.0, 1.0, 1.0)),
    )


def localized_curvature_defect(
    curvature: float,
    width: float,
    *,
    centre: float = 0.53125,
    across: float = 0.02,
    breadth: float = 1.0,
    spans: int = 16,
    cluster: int = 8,
    degree_v: int = 3,
) -> tuple[Face, Face, dict[str, Any]]:
    """A G1 seam whose curvature jump is a bump function of position along the seam.

    ``x(u) = across * u`` and ``z(u, v) = kappa across^2 u^2 G(v) / 2``, so the normal
    curvature across the seam is ``kappa * G(v)`` exactly while the deflection scales with
    ``across^2`` and stays small: a localized curvature discontinuity without the surface
    swinging out of plane.  ``G`` is a cubic B-spline through samples of a Gaussian of
    scale ``width`` (fraction of the seam length) centred at ``centre`` -- by default half
    way between two samples of the 17-point production grid, i.e. the worst case.

    The facing side is a straight plate, so the jump is ``kappa`` wherever ``G`` is
    nonzero and the seam is G1 by construction; how much of that the classifier sees is
    what the caller measures.

    The knot vector is uniform over ``spans`` spans *plus* a cluster of ``cluster`` spans
    around ``centre``. Without the cluster a narrow bump would be limited by the control
    net rather than by the sampler, and the audit would be measuring the wrong thing:
    the defect has to exist to be missed.
    """
    knots_v = _clustered_knots(spans, degree_v, centre, 3.0 * width, cluster)
    parameters = _greville(knots_v, degree_v)
    envelope = [math.exp(-((t - centre) / width) ** 2) for t in parameters]
    amplitude = 0.5 * curvature * across * across
    net: list[list[Vec3]] = []
    for column, x in enumerate((0.0, 0.5 * across, across)):
        row: list[Vec3] = []
        for index, t in enumerate(parameters):
            z = amplitude * envelope[index] if column == 2 else 0.0
            row.append((x, t * breadth, z))
        net.append(row)
    outgoing = NURBSSurface(
        degree=(2, degree_v),
        control_net=net,
        weights=[[1.0] * len(parameters) for _ in net],
        knots=((0.0, 0.0, 0.0, 1.0, 1.0, 1.0), knots_v),
    )
    info = {
        "curvature": curvature,
        "width": width,
        "centre": centre,
        "across": across,
        "spans_v": spans,
        "cluster": cluster,
        "envelope_peak": max(envelope),
        "label": "G1",
    }
    return _face(straight_plate(across, breadth, incoming=True), 0), _face(outgoing, 1), info


# --------------------------------------------------------------------------- #
# measurement wrapper
# --------------------------------------------------------------------------- #
def measure(
    face_a: Face,
    face_b: Face,
    *,
    samples: int,
    scale: float = 1.0,
    adaptive: bool = False,
    max_samples: int = 513,
    angle_limit_deg: float = DEFAULT_ANGLE_LIMIT_DEG,
    curvature_limit: float = DEFAULT_CURVATURE_LIMIT,
) -> dict[str, Any]:
    """``seam_metrics`` (optionally adaptive) plus the verdict it implies."""
    if adaptive:
        metrics = seam_metrics_adaptive(
            face_a, face_b, samples=samples, max_samples=max_samples, scale=scale
        )
    else:
        metrics = seam_metrics(face_a, face_b, samples=samples, scale=scale)
    verdict = classify(
        metrics, angle_limit_deg=angle_limit_deg, curvature_limit=curvature_limit
    )
    return {
        "samples": int(metrics.get("samples_used", samples)),
        "converged": bool(metrics.get("converged", True)),
        "max_normal_angle_deg": metrics["max_normal_angle_deg"],
        "max_curvature_jump": metrics["max_curvature_jump"],
        "max_curvature_jump_absolute": metrics.get(
            "max_curvature_jump_absolute", float("nan")
        ),
        "classification": verdict,
        # A seam is "flagged" when the classifier refuses to call it smooth.
        "flagged": verdict != "G2",
    }


# --------------------------------------------------------------------------- #
# study 1: where the "relative" criterion is really absolute
# --------------------------------------------------------------------------- #
def curvature_scale_sweep(
    *,
    levels: tuple[float, ...] = (1e-4, 1e-3, 3e-3, 1e-2, 3e-2, 0.1, 0.3, 1.0, 3.0, 10.0),
    length: float = 1.0,
    breadth: float = 1.0,
    mismatch: float = 1.0,
    samples: int = 65,
    scale: float = 1.0,
    angle_limit_deg: float = DEFAULT_ANGLE_LIMIT_DEG,
    curvature_limit: float = DEFAULT_CURVATURE_LIMIT,
) -> list[dict[str, Any]]:
    """Classification of a G1 seam at a prescribed *relative* curvature mismatch.

    ``mismatch`` is ``(k_b - k_a) / k_b``: ``1.0`` puts a curved plate against a flat one,
    ``0.3`` is the "curvatures off by 30 %" case the module docstring quotes.  ``levels``
    are the absolute curvature ``k_b`` in 1/length; ``scale`` is the seam scale the
    criterion normalises by (production passes the model length).
    """
    if not 0.0 < mismatch <= 1.0:
        raise ValueError("mismatch must be in (0, 1]")
    rows: list[dict[str, Any]] = []
    for level in levels:
        curved = curved_plate(level, length, breadth, incoming=False)
        if mismatch == 1.0:
            flat = straight_plate(length, breadth, incoming=True)
        else:
            flat = curved_plate(level * (1.0 - mismatch), length, breadth, incoming=True)
        measured = measure(
            _face(flat, 0),
            _face(curved, 1),
            samples=samples,
            scale=scale,
            angle_limit_deg=angle_limit_deg,
            curvature_limit=curvature_limit,
        )
        rows.append(
            {
                "level": level,
                "level_times_scale": level * scale,
                "mismatch": mismatch,
                "expected_curvature_jump": level * mismatch,
                **measured,
            }
        )
    return rows


def detection_threshold(
    rows: list[dict[str, Any]], *, key: str = "max_curvature_jump"
) -> float | None:
    """Smallest ``level_times_scale`` in ``rows`` whose seam is flagged, if any."""
    for row in sorted(rows, key=lambda item: item["level_times_scale"]):
        if row[key] > DEFAULT_CURVATURE_LIMIT:
            return row["level_times_scale"]
    return None


# --------------------------------------------------------------------------- #
# study 2: sampling aliasing along the seam
# --------------------------------------------------------------------------- #
def sampling_alias_sweep(
    *,
    widths: tuple[float, ...] = (1 / 8, 1 / 16, 1 / 32, 1 / 64, 1 / 128, 1 / 256),
    curvature: float = 2.0,
    centres: tuple[float, ...] = (0.53125, 0.5),
    samples: int = PRODUCTION_SAMPLES,
    dense: int = DENSE_SAMPLES,
    max_samples: int = 513,
    scale: float = 1.0,
    angle_limit_deg: float = DEFAULT_ANGLE_LIMIT_DEG,
    curvature_limit: float = DEFAULT_CURVATURE_LIMIT,
) -> list[dict[str, Any]]:
    """Measure a localized defect at the dense truth, the default grid and adaptively."""
    rows: list[dict[str, Any]] = []
    for centre in centres:
        for width in widths:
            face_a, face_b, info = localized_curvature_defect(
                curvature, width, centre=centre
            )
            truth = measure(
                face_a,
                face_b,
                samples=dense,
                scale=scale,
                angle_limit_deg=angle_limit_deg,
                curvature_limit=curvature_limit,
            )
            coarse = measure(
                face_a,
                face_b,
                samples=samples,
                scale=scale,
                angle_limit_deg=angle_limit_deg,
                curvature_limit=curvature_limit,
            )
            adaptive = measure(
                face_a,
                face_b,
                samples=samples,
                max_samples=max_samples,
                adaptive=True,
                scale=scale,
                angle_limit_deg=angle_limit_deg,
                curvature_limit=curvature_limit,
            )
            rows.append(
                {
                    "width": width,
                    "centre": centre,
                    "curvature": info["curvature"],
                    "dense": truth,
                    "coarse": coarse,
                    "adaptive": adaptive,
                    "coarse_missed": truth["flagged"] and not coarse["flagged"],
                    "adaptive_recovered": truth["flagged"] and adaptive["flagged"],
                    "coarse_capture": (
                        coarse["max_curvature_jump"] / truth["max_curvature_jump"]
                        if truth["max_curvature_jump"] > 0
                        else float("nan")
                    ),
                }
            )
    return rows


# --------------------------------------------------------------------------- #
# study 3: does adaptivity disturb the calibrated family?
# --------------------------------------------------------------------------- #
def adaptive_regression(
    *,
    cases: int = 100,
    samples: int = PRODUCTION_SAMPLES,
    max_samples: int = 513,
    scale: float = 1.0,
    seed: int = 20261107,
) -> dict[str, Any]:
    """Re-classify the production seam family with adaptive sampling.

    Adaptivity is only safe if it does not move verdicts on the family the thresholds
    were calibrated on; anything else would silently require recalibration.
    """
    import random

    from .continuity import make_seam_case

    generator = random.Random(seed)
    kinds = ("G0", "G1", "G2")
    changed: list[dict[str, Any]] = []
    samples_used: list[int] = []
    total = 0
    for index in range(cases):
        kind = kinds[index % 3]
        surface_a, surface_b, _ = make_seam_case(kind, generator, scale=scale)
        face_a, face_b = _face(surface_a, 0), _face(surface_b, 1)
        fixed = measure(face_a, face_b, samples=samples, scale=scale)
        adaptive = measure(
            face_a, face_b, samples=samples, max_samples=max_samples, adaptive=True, scale=scale
        )
        samples_used.append(adaptive["samples"])
        total += 1
        if fixed["classification"] != adaptive["classification"]:
            changed.append(
                {
                    "index": index,
                    "kind": kind,
                    "fixed": fixed["classification"],
                    "adaptive": adaptive["classification"],
                    "fixed_jump": fixed["max_curvature_jump"],
                    "adaptive_jump": adaptive["max_curvature_jump"],
                }
            )
    samples_used.sort()
    return {
        "cases": total,
        "changed": len(changed),
        "changed_rows": changed,
        "samples_used_median": samples_used[len(samples_used) // 2],
        "samples_used_max": samples_used[-1],
        "samples_used_mean": sum(samples_used) / len(samples_used),
    }


def ood_audit(*, dense: int = DENSE_SAMPLES, max_samples: int = 513) -> dict[str, Any]:
    """Every study in one dictionary, as the example writes it to ``results.json``."""
    from .continuity import continuity_evidence

    baseline = continuity_evidence(cases=60)
    scale_rows = {
        "full_mismatch": curvature_scale_sweep(mismatch=1.0),
        "thirty_percent": curvature_scale_sweep(mismatch=0.3),
    }
    alias_rows = sampling_alias_sweep(dense=dense, max_samples=max_samples)
    regression = adaptive_regression(max_samples=max_samples)
    return {
        "baseline": baseline,
        "curvature_scale": scale_rows,
        "curvature_scale_thresholds": {
            name: detection_threshold(rows) for name, rows in scale_rows.items()
        },
        "sampling_alias": alias_rows,
        "sampling_alias_summary": {
            "coarse_missed": sum(1 for row in alias_rows if row["coarse_missed"]),
            "adaptive_recovered": sum(1 for row in alias_rows if row["adaptive_recovered"]),
            "cases": len(alias_rows),
            "narrowest_width_coarse_missed": min(
                (row["width"] for row in alias_rows if row["coarse_missed"]), default=None
            ),
        },
        "adaptive_regression": regression,
        "settings": {
            "dense_samples": dense,
            "production_samples": PRODUCTION_SAMPLES,
            "max_samples": max_samples,
            "angle_limit_deg": DEFAULT_ANGLE_LIMIT_DEG,
            "curvature_limit": DEFAULT_CURVATURE_LIMIT,
        },
    }