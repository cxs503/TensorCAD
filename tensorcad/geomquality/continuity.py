"""Curvature continuity (G0/G1/G2) of CAD patch seams, with calibrated thresholds.

A fair hull has second-derivative continuity across every patch seam that is not a
deliberate knuckle.  A *nearly* G1 seam -- tangents aligned to five decimal places but
curvatures off by 30 % -- produces exactly the sort of ripple that a towing tank sees
and a thumbnail hides, and it survives every topology check.  So the seam is measured
with three numbers and classified:

``position gap``
    ``|S_a(seam) - S_b(seam)|``; a nonzero gap is a crack (handled by the topology
    audit) and here only serves as a sanity check.
``normal angle``
    ``angle(n_a, n_b)`` along the seam, in degrees.  Above the G1 limit the seam is
    **G0**: a tangent discontinuity.
``curvature jump``
    ``|kappa_a - kappa_b|`` of the normal curvature in the direction *across* the
    seam, relative to the larger of the two and to ``1 / scale``.  Above the G2 limit
    the seam is **G1**: tangent-continuous but curvature-discontinuous.

Gaussian and mean curvature distributions (:func:`curvature_distribution`) are
reported alongside from the analytic second derivatives, so the absolute curvature
level of the patch is visible and a jump can be read in context rather than as a bare
ratio.

The two limits are **calibrated** (:func:`calibrate`) on a labelled seam set and then
scored on a *held-out* set (:func:`score`), because a threshold chosen on the same
cases it is scored on is a tautology.
"""
from __future__ import annotations

import math
import random
from collections.abc import Sequence
from dataclasses import dataclass
from typing import Any

from ..nurbs import NURBSSurface
from ._eval import evaluate_fast, jacobian_fast
from .brep import BRep, Face

Vec3 = tuple[float, float, float]

#: Default classification limits; :func:`calibrate` replaces them with measured ones.
DEFAULT_ANGLE_LIMIT_DEG = 1.0
DEFAULT_CURVATURE_LIMIT = 0.10

#: A seam whose normals differ by more than this is a *designed* crease (deck edge, keel
#: knuckle, transom corner), not a continuity defect: it is reported as ``KNUCKLE`` and
#: no smoothness criterion is applied to it.
CREASE_LIMIT_DEG = 30.0


def _sub(a: Sequence[float], b: Sequence[float]) -> Vec3:
    return (a[0] - b[0], a[1] - b[1], a[2] - b[2])


def _dot(a: Sequence[float], b: Sequence[float]) -> float:
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]


def _norm(a: Sequence[float]) -> float:
    return math.sqrt(_dot(a, a))


def _angle_deg(a: Sequence[float], b: Sequence[float]) -> float:
    denominator = _norm(a) * _norm(b)
    if denominator == 0.0:
        return 0.0
    cosine = max(-1.0, min(1.0, _dot(a, b) / denominator))
    return math.degrees(math.acos(cosine))


def normal_curvature(surface: NURBSSurface, u: float, v: float, along: str = "u") -> float:
    """Curvature of the surface normal section in the ``u`` or ``v`` direction."""
    derivatives = surface.derivatives(u, v, 2)
    s_u, s_v = derivatives[(1, 0)], derivatives[(0, 1)]
    normal = _sub_cross(s_u, s_v)
    length = _norm(normal)
    if length == 0.0:
        return float("nan")
    normal = (normal[0] / length, normal[1] / length, normal[2] / length)
    if along == "u":
        second, first = derivatives[(2, 0)], s_u
    else:
        second, first = derivatives[(0, 2)], s_v
    speed = _dot(first, first)
    if speed == 0.0:
        return float("nan")
    return _dot(second, normal) / speed


def _sub_cross(a: Sequence[float], b: Sequence[float]) -> Vec3:
    return (
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    )


def curvature_distribution(
    surface: NURBSSurface, *, samples: int = 9
) -> dict[str, Any]:
    """Gaussian and mean curvature statistics over the patch (Gauss-point sampling)."""
    (u0, u1), (v0, v1) = surface.domain
    gaussian: list[float] = []
    mean: list[float] = []
    for i in range(samples):
        u = u0 + (u1 - u0) * (i + 0.5) / samples
        for j in range(samples):
            v = v0 + (v1 - v0) * (j + 0.5) / samples
            d = surface.derivatives(u, v, 2)
            s_u, s_v = d[(1, 0)], d[(0, 1)]
            s_uu, s_uv, s_vv = d[(2, 0)], d[(1, 1)], d[(0, 2)]
            normal = _sub_cross(s_u, s_v)
            length = _norm(normal)
            if length == 0.0:
                continue
            unit = (normal[0] / length, normal[1] / length, normal[2] / length)
            big_e, big_f, big_g = _dot(s_u, s_u), _dot(s_u, s_v), _dot(s_v, s_v)
            little_e = _dot(s_uu, unit)
            little_f = _dot(s_uv, unit)
            little_g = _dot(s_vv, unit)
            denominator = big_e * big_g - big_f * big_f
            if denominator <= 0.0:
                continue
            gaussian.append((little_e * little_g - little_f * little_f) / denominator)
            mean.append(
                (big_e * little_g - 2.0 * big_f * little_f + big_g * little_e)
                / (2.0 * denominator)
            )
    if not gaussian:
        return {"samples": 0}
    return {
        "samples": len(gaussian),
        "gaussian": _stats(gaussian),
        "mean": _stats(mean),
    }


def _stats(values: Sequence[float]) -> dict[str, float]:
    ordered = sorted(values)
    count = len(ordered)
    mean = sum(ordered) / count
    variance = sum((value - mean) ** 2 for value in ordered) / count
    return {
        "min": ordered[0],
        "max": ordered[-1],
        "mean": mean,
        "rms": math.sqrt(variance),
        "p05": ordered[max(0, int(0.05 * count) - 1)],
        "p95": ordered[min(count - 1, int(0.95 * count))],
    }


@dataclass
class SeamReport:
    """One measured patch seam."""

    face_a: int
    face_b: int
    samples: int
    max_gap: float
    max_normal_angle_deg: float
    max_curvature_jump: float
    mean_curvature_jump: float
    classification: str = "G2"

    def to_dict(self) -> dict[str, Any]:
        return {
            "face_a": self.face_a,
            "face_b": self.face_b,
            "max_gap": self.max_gap,
            "max_normal_angle_deg": self.max_normal_angle_deg,
            "max_curvature_jump": self.max_curvature_jump,
            "mean_curvature_jump": self.mean_curvature_jump,
            "classification": self.classification,
        }


def _boundary(face: Face, side: int, samples: int) -> tuple[list[Vec3], list[tuple[float, float]]]:
    """Sample one boundary isocurve of a face: ``(points, parameters)``."""
    (u0, u1), (v0, v1) = face.domain
    points: list[Vec3] = []
    parameters: list[tuple[float, float]] = []
    for k in range(samples):
        t = k / (samples - 1)
        if side == 0:
            parameter = (u0 + (u1 - u0) * t, v0)
        elif side == 1:
            parameter = (u1, v0 + (v1 - v0) * t)
        elif side == 2:
            parameter = (u1 - (u1 - u0) * t, v1)
        else:
            parameter = (u0, v1 - (v1 - v0) * t)
        points.append(evaluate_fast(face.surface, *parameter))
        parameters.append(parameter)
    return points, parameters


def _along(side: int) -> str:
    """Direction *across* the seam: opposite to the isocurve's own variation."""
    return "v" if side in (0, 2) else "u"


def _unit_normal(surface: NURBSSurface, parameter: tuple[float, float]) -> Vec3:
    """Unit normal on the fast path; raises when the parameterisation degenerates."""
    _, _, normal = jacobian_fast(surface, *parameter)
    length = math.sqrt(normal[0] ** 2 + normal[1] ** 2 + normal[2] ** 2)
    if length == 0.0:
        raise ValueError("degenerate parameterisation: no normal")
    return (normal[0] / length, normal[1] / length, normal[2] / length)


def seam_metrics(
    face_a: Face,
    face_b: Face,
    *,
    samples: int = 17,
    scale: float = 1.0,
    tolerance: float | None = None,
) -> dict[str, float]:
    """Measure the best-matching pair of boundary isocurves between two faces.

    The four isocurves of each face are sampled once; the pair whose endpoints coincide
    (within ``tolerance``, default ``1e-3 * scale``) and whose samples agree best is the
    seam.  At that seam the normal angle and the curvature jump in the across direction
    are evaluated, skipping samples where the parameterisation degenerates (a collapsed
    bow boundary has zero Jacobian and no normal to compare).
    """
    limit = tolerance if tolerance is not None else 1e-3 * max(scale, 1e-9)
    curves_a = {side: _boundary(face_a, side, samples) for side in range(4)}
    curves_b = {side: _boundary(face_b, side, samples) for side in range(4)}
    best: dict[str, Any] | None = None
    for side_a, (points_a, _) in curves_a.items():
        for side_b, (points_b, _) in curves_b.items():
            forward = max(math.dist(p, q) for p, q in zip(points_a, points_b))
            backward = max(math.dist(p, q) for p, q in zip(points_a, reversed(points_b)))
            for direction, gap in (("forward", forward), ("backward", backward)):
                if gap > limit:
                    continue
                if best is None or gap < best["max_gap"]:
                    best = {
                        "side_a": side_a,
                        "side_b": side_b,
                        "direction": direction,
                        "max_gap": gap,
                    }
    if best is None:
        return {
            "max_gap": math.inf,
            "max_normal_angle_deg": math.inf,
            "max_curvature_jump": math.inf,
            "mean_curvature_jump": math.inf,
            "samples": float(samples),
            "side_a": -1,
            "side_b": -1,
        }
    points_a, parameters_a = curves_a[best["side_a"]]
    points_b, parameters_b = curves_b[best["side_b"]]
    if best["direction"] == "backward":
        points_b = list(reversed(points_b))
        parameters_b = list(reversed(parameters_b))
    gaps: list[float] = []
    angles: list[float] = []
    jumps: list[float] = []
    for point_a, parameter_a, point_b, parameter_b in zip(
        points_a, parameters_a, points_b, parameters_b
    ):
        try:
            normal_a = _unit_normal(face_a.surface, parameter_a)
            normal_b = _unit_normal(face_b.surface, parameter_b)
        except Exception:
            continue  # degenerate isocurve point: no tangent plane exists to compare
        # A face with same_sense = False has the reverse of S_u x S_v as its outward
        # normal, and the sign of the normal curvature follows the same flip. Comparing
        # the raw surface normals would turn a perfectly smooth mirror seam into a 180 deg
        # "flip" and a spurious curvature jump of exactly 2.
        sign_a = 1.0 if face_a.same_sense else -1.0
        sign_b = 1.0 if face_b.same_sense else -1.0
        normal_a = tuple(sign_a * value for value in normal_a)
        normal_b = tuple(sign_b * value for value in normal_b)
        gaps.append(math.dist(point_a, point_b))
        angles.append(_angle_deg(normal_a, normal_b))
        kappa_a = sign_a * normal_curvature(face_a.surface, *parameter_a, along=_along(best["side_a"]))
        kappa_b = sign_b * normal_curvature(face_b.surface, *parameter_b, along=_along(best["side_b"]))
        if math.isfinite(kappa_a) and math.isfinite(kappa_b):
            reference = max(abs(kappa_a), abs(kappa_b), 1.0 / max(scale, 1e-9))
            jumps.append(abs(kappa_a - kappa_b) / reference)
    return {
        "max_gap": max(gaps) if gaps else math.inf,
        "max_normal_angle_deg": max(angles) if angles else math.inf,
        "max_curvature_jump": max(jumps) if jumps else math.inf,
        "mean_curvature_jump": sum(jumps) / len(jumps) if jumps else math.inf,
        "samples": float(len(gaps)),
        "side_a": float(best["side_a"]),
        "side_b": float(best["side_b"]),
    }


def classify(
    metrics: dict[str, float],
    *,
    angle_limit_deg: float = DEFAULT_ANGLE_LIMIT_DEG,
    curvature_limit: float = DEFAULT_CURVATURE_LIMIT,
) -> str:
    """``KNUCKLE`` (designed crease), ``G0``, ``G1`` or ``G2`` (continuous)."""
    if not math.isfinite(metrics.get("max_gap", math.inf)):
        return "NO_SEAM"
    if metrics["max_normal_angle_deg"] > CREASE_LIMIT_DEG:
        return "KNUCKLE"
    if metrics["max_normal_angle_deg"] > angle_limit_deg:
        return "G0"
    if metrics["max_curvature_jump"] > curvature_limit:
        return "G1"
    return "G2"


def audit_seams(
    brep: BRep,
    *,
    angle_limit_deg: float = DEFAULT_ANGLE_LIMIT_DEG,
    curvature_limit: float = DEFAULT_CURVATURE_LIMIT,
    samples: int = 17,
    tolerance: float = 1e-3,
) -> list[SeamReport]:
    """Classify every pair of faces that share an edge in the adjacency graph."""
    by_index = {face.index: face for face in brep.faces}
    reports: list[SeamReport] = []
    seen: set[tuple[int, int]] = set()
    for edge_index, uses in brep.edge_uses().items():
        if len(uses) != 2:
            continue
        (face_a, _), (face_b, _) = uses
        key = (min(face_a, face_b), max(face_a, face_b))
        if key in seen:
            continue
        seen.add(key)
        scale = max(brep.length, 1e-9)
        metrics = seam_metrics(
            by_index[face_a],
            by_index[face_b],
            samples=samples,
            scale=scale,
            tolerance=tolerance * scale,
        )
        if not math.isfinite(metrics["max_gap"]):
            continue
        reports.append(
            SeamReport(
                face_a=face_a,
                face_b=face_b,
                samples=int(metrics["samples"]),
                max_gap=metrics["max_gap"],
                max_normal_angle_deg=metrics["max_normal_angle_deg"],
                max_curvature_jump=metrics["max_curvature_jump"],
                mean_curvature_jump=metrics["mean_curvature_jump"],
                classification=classify(
                    metrics,
                    angle_limit_deg=angle_limit_deg,
                    curvature_limit=curvature_limit,
                ),
            )
        )
    return reports


# --------------------------------------------------------------------------- #
# labelled seam set: construction + calibration + scoring
# --------------------------------------------------------------------------- #
def _bezier_plate(control: Sequence[tuple[float, float]], width: float) -> NURBSSurface:
    """Plate patch ``S(u,v) = (X(u), v*width, Z(u))`` with a planar cubic profile."""
    control_net = [[(x, 0.0, z), (x, width, z)] for x, z in control]
    return NURBSSurface(
        degree=(3, 1),
        control_net=control_net,
        weights=[[1.0, 1.0] for _ in control_net],
        knots=[[0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0], [0.0, 0.0, 1.0, 1.0]],
    )


def make_seam_case(
    kind: str, rng: random.Random, *, scale: float = 1.0
) -> tuple[NURBSSurface, NURBSSurface, dict[str, Any]]:
    """Two plate patches joined at ``u = 1`` with a *known* continuity class.

    The second patch is built from the joint derivatives of the first:

    ``G2``  same tangent *and* same second derivative as a vector (equal speed, so the
            curvature matches exactly);
    ``G1``  same tangent, second control point displaced off the tangent line, which
            changes the curvature by a controlled amount;
    ``G0``  tangent rotated by a random angle, so the normals are discontinuous.

    Only the class label is stored; the measured metrics come from the surfaces alone.
    """
    if kind not in {"G0", "G1", "G2"}:
        raise ValueError(f"unknown seam class {kind!r}")
    length = scale * rng.uniform(0.5, 1.5)
    angle = rng.uniform(-1.0, 1.0)
    bend = rng.uniform(0.2, 2.0) * rng.choice((-1.0, 1.0)) * scale
    direction = (math.cos(angle), math.sin(angle))
    normal_dir = (-math.sin(angle), math.cos(angle))
    step = length / 3.0
    # First patch, ending at the origin: A'(1) = 3*step*direction, and the second
    # derivative set by `bend`.
    p3 = (0.0, 0.0)
    p2 = (p3[0] - step * direction[0], p3[1] - step * direction[1])
    second_derivative = (6.0 * bend * normal_dir[0], 6.0 * bend * normal_dir[1])
    p1 = (
        2.0 * p2[0] - p3[0] + second_derivative[0] / 6.0,
        2.0 * p2[1] - p3[1] + second_derivative[1] / 6.0,
    )
    p0 = (p1[0] - 2.0 * step * direction[0], p1[1] - 2.0 * step * direction[1])
    first = _bezier_plate([p0, p1, p2, p3], width=1.0)

    extend = rng.uniform(0.6, 1.6) * scale
    unit = direction
    if kind == "G0":
        twist = math.radians(rng.uniform(2.0, 25.0)) * rng.choice((-1.0, 1.0))
        unit = (
            direction[0] * math.cos(twist) - direction[1] * math.sin(twist),
            direction[0] * math.sin(twist) + direction[1] * math.cos(twist),
        )
    q0 = p3
    q1 = (q0[0] + unit[0] * step, q0[1] + unit[1] * step)
    if kind == "G2":
        # A G2 joint still carries modelling noise: perturb at the 1e-5 scale, i.e.
        # below the curvature limit, so the classifier is tested at its noise floor
        # instead of only on machine-exact seams.
        noise = rng.uniform(0.0, 5e-5) * scale * rng.choice((-1.0, 1.0))
        q2 = (
            2.0 * q1[0] - q0[0] + second_derivative[0] / 6.0 + noise * normal_dir[0],
            2.0 * q1[1] - q0[1] + second_derivative[1] / 6.0 + noise * normal_dir[1],
        )
    else:
        offset = rng.uniform(0.002, 1.2) * scale * rng.choice((-1.0, 1.0))
        q2 = (2.0 * q1[0] - q0[0] + offset * normal_dir[0], 2.0 * q1[1] - q0[1] + offset * normal_dir[1])
    q3 = (q2[0] + unit[0] * extend, q2[1] + unit[1] * extend)
    second = _bezier_plate([q0, q1, q2, q3], width=1.0)
    info = {
        "kind": kind,
        "bend": bend,
        "extend": extend,
        "twist_deg": 0.0 if kind != "G0" else math.degrees(
            math.atan2(unit[1], unit[0]) - math.atan2(direction[1], direction[0])
        ),
    }
    return first, second, info


def label_seam_dataset(
    cases: int = 100, seed: int = 20261107, scale: float = 1.0
) -> list[dict[str, Any]]:
    """Measure a random seam set with known G0/G1/G2 labels."""
    rng = random.Random(seed)
    rows: list[dict[str, Any]] = []
    kinds = ("G2", "G1", "G0")
    for index in range(cases):
        kind = kinds[index % 3]
        surface_a, surface_b, info = make_seam_case(kind, rng, scale=scale)
        metrics = seam_metrics(
            Face(index=0, surface=surface_a, boundary=()),
            Face(index=1, surface=surface_b, boundary=()),
            samples=17,
            scale=scale,
        )
        row = dict(info)
        row.update(metrics)
        rows.append(row)
    return rows


def _f1(truth: Sequence[bool], predicted: Sequence[bool]) -> dict[str, float]:
    true_positive = sum(1 for t, p in zip(truth, predicted) if t and p)
    false_positive = sum(1 for t, p in zip(truth, predicted) if not t and p)
    false_negative = sum(1 for t, p in zip(truth, predicted) if t and not p)
    precision = true_positive / (true_positive + false_positive) if true_positive + false_positive else 0.0
    recall = true_positive / (true_positive + false_negative) if true_positive + false_negative else 0.0
    f1 = 2 * precision * recall / (precision + recall) if precision + recall else 0.0
    return {
        "precision": precision,
        "recall": recall,
        "f1": f1,
        "tp": true_positive,
        "fp": false_positive,
        "fn": false_negative,
    }


def calibrate(
    calibration: Sequence[dict[str, Any]],
    *,
    angle_grid: Sequence[float] = (0.05, 0.1, 0.2, 0.5, 1.0, 2.0, 4.0),
    curvature_grid: Sequence[float] = (0.01, 0.02, 0.05, 0.1, 0.2, 0.35, 0.5),
) -> dict[str, Any]:
    """Grid search for the (angle, curvature) limits that maximise F1 on calibrators."""
    truth = [row["kind"] != "G2" for row in calibration]
    best: dict[str, Any] = {"f1": -1.0}
    for angle_limit in angle_grid:
        for curvature_limit in curvature_grid:
            predicted = [
                row["max_normal_angle_deg"] > angle_limit
                or row["max_curvature_jump"] > curvature_limit
                for row in calibration
            ]
            score = _f1(truth, predicted)
            if score["f1"] > best["f1"]:
                best = {
                    "angle_limit_deg": angle_limit,
                    "curvature_limit": curvature_limit,
                    "f1": score["f1"],
                    "precision": score["precision"],
                    "recall": score["recall"],
                }
    return best


def score(
    rows: Sequence[dict[str, Any]],
    *,
    angle_limit_deg: float,
    curvature_limit: float,
) -> dict[str, Any]:
    """Score the calibrated limits on held-out seams, overall and per class."""
    truth = [row["kind"] != "G2" for row in rows]
    predicted = [
        row["max_normal_angle_deg"] > angle_limit_deg
        or row["max_curvature_jump"] > curvature_limit
        for row in rows
    ]
    overall = _f1(truth, predicted)
    # Curvature-only view: does the curvature criterion alone separate G1 from G2?
    subset = [row for row in rows if row["kind"] in {"G1", "G2"}]
    curvature_only = _f1(
        [row["kind"] == "G1" for row in subset],
        [row["max_curvature_jump"] > curvature_limit for row in subset],
    )
    per_class: dict[str, dict[str, float]] = {}
    for kind in ("G0", "G1", "G2"):
        group = [row for row in rows if row["kind"] == kind]
        if not group:
            continue
        detected = sum(
            1
            for row in group
            if (row["max_normal_angle_deg"] > angle_limit_deg)
            or (row["max_curvature_jump"] > curvature_limit)
            if kind != "G2"
        )
        if kind == "G2":
            false_alarms = sum(
                1
                for row in group
                if (row["max_normal_angle_deg"] > angle_limit_deg)
                or (row["max_curvature_jump"] > curvature_limit)
            )
            per_class[kind] = {
                "count": float(len(group)),
                "false_alarm_rate": false_alarms / len(group),
            }
        else:
            per_class[kind] = {
                "count": float(len(group)),
                "detection_rate": detected / len(group),
            }
    return {
        "overall": overall,
        "curvature_only": curvature_only,
        "per_class": per_class,
        "angle_limit_deg": angle_limit_deg,
        "curvature_limit": curvature_limit,
    }


def continuity_evidence(
    *, cases: int = 100, calibration_fraction: float = 0.6, seed: int = 20261107
) -> dict[str, Any]:
    """End-to-end evidence: calibrate on one split, score on the other."""
    rows = label_seam_dataset(cases=cases, seed=seed)
    split = int(len(rows) * calibration_fraction)
    calibration, test = rows[:split], rows[split:]
    best = calibrate(calibration)
    scored = score(
        test,
        angle_limit_deg=best["angle_limit_deg"],
        curvature_limit=best["curvature_limit"],
    )
    jumps_g1 = [row["max_curvature_jump"] for row in test if row["kind"] == "G1"]
    jumps_g2 = [row["max_curvature_jump"] for row in test if row["kind"] == "G2"]
    angles_g0 = [row["max_normal_angle_deg"] for row in test if row["kind"] == "G0"]
    angles_smooth = [
        row["max_normal_angle_deg"] for row in test if row["kind"] in {"G1", "G2"}
    ]
    return {
        "cases": len(rows),
        "calibration_cases": len(calibration),
        "test_cases": len(test),
        "separation": {
            "curvature_jump_min_G1": min(jumps_g1) if jumps_g1 else None,
            "curvature_jump_max_G2": max(jumps_g2) if jumps_g2 else None,
            "normal_angle_min_G0_deg": min(angles_g0) if angles_g0 else None,
            "normal_angle_max_smooth_deg": max(angles_smooth) if angles_smooth else None,
        },
        "default_threshold_score": score(
            test,
            angle_limit_deg=DEFAULT_ANGLE_LIMIT_DEG,
            curvature_limit=DEFAULT_CURVATURE_LIMIT,
        )["overall"],
        "thresholds": {
            "angle_limit_deg": best["angle_limit_deg"],
            "curvature_limit": best["curvature_limit"],
        },
        "calibration_f1": best["f1"],
        "test": scored,
        "seed": seed,
    }
