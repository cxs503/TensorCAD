#!/usr/bin/env python3
"""Audit the three aspect-ratio dialects on real meshes and write the evidence artefact.

What this establishes, and how it can fail:

1. **Identity** -- on every mesh, the ``legacy`` value equals what
   ``tensorcad.validation.measure_quality`` reports and the ``edge_form``/**verdict**
   values equal what ``tensorcad.geomquality.quality`` / ``tensorcad.meshquality``
   report. Tolerance 1e-12 relative. A drift here means a fourth dialect appeared.
2. **Bounds** -- ``alpha = edge_form/legacy`` is inside ``[1/sqrt(3), sqrt(3)/2]`` and
   ``beta = verdict/edge_form >= sqrt(3)/2`` on every element of every mesh.
3. **Consequence** -- at the copied limit 25 the foreign-dialect gates pass elements the
   legacy gate rejects, and at :func:`conservative_limit` they pass none. The second
   number is 0 by proof; it is recomputed here as a self-check on the proof.

Writes ``artifacts/mesh_quality_dialects/audit.json``, ``.../dialects.png`` and
``.../manifest.json`` (sha256 of the JSON), and exits non-zero if any criterion fails.
"""
from __future__ import annotations

import hashlib
import json
import math
import random
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from tensorcad.geomquality.quality import triangle_quality  # noqa: E402
from tensorcad.hull import HullParameters, build_mesh  # noqa: E402
from tensorcad.meshquality import dialects as D  # noqa: E402
from tensorcad.meshquality import tri_metrics  # noqa: E402
from tensorcad.validation import measure_quality  # noqa: E402

ARTIFACT_DIR = REPO_ROOT / "artifacts" / "mesh_quality_dialects"
LIMIT = 25.0
TOL = 1e-12


def box_surface() -> tuple[list[tuple[float, float, float]], list[tuple[int, int, int]]]:
    points = [
        (0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (1.0, 1.0, 0.0), (0.0, 1.0, 0.0),
        (0.0, 0.0, 1.0), (1.0, 0.0, 1.0), (1.0, 1.0, 1.0), (0.0, 1.0, 1.0),
    ]
    faces: list[tuple[int, int, int]] = []
    for a, b, c, d in ((0, 3, 2, 1), (4, 5, 6, 7), (0, 1, 5, 4), (1, 2, 6, 5), (2, 3, 7, 6), (3, 0, 4, 7)):
        faces.append((a, b, c))
        faces.append((a, c, d))
    return points, faces


def hull_surface(nx: int, nz: int) -> tuple[list[tuple[float, float, float]], list[tuple[int, int, int]]]:
    mesh, _ = build_mesh(HullParameters(length=10.0, beam=2.0, draft=0.8, nx=nx, nz=nz))
    return list(mesh.vertices), [tuple(face) for face in mesh.faces]


def sliver_patch() -> tuple[list[tuple[float, float, float]], list[tuple[int, int, int]]]:
    """A fan whose outer ring is 400x wider than tall: the shape the gate exists for."""
    points = [(0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (0.0, 400.0, 0.0), (1.0, 400.0, 0.0)]
    return points, [(0, 1, 3), (0, 3, 2)]


def random_soup(count: int = 400) -> tuple[list[tuple[float, float, float]], list[tuple[int, int, int]]]:
    rng = random.Random(20261113)
    points: list[tuple[float, float, float]] = []
    faces: list[tuple[int, int, int]] = []
    for _ in range(count):
        base = (rng.uniform(-2.0, 2.0), rng.uniform(-2.0, 2.0), rng.uniform(-2.0, 2.0))
        scale = 10.0 ** rng.uniform(-3.0, 3.0)
        start = len(points)
        points.append(base)
        points.append(tuple(base[i] + scale * rng.uniform(-1.0, 1.0) for i in range(3)))  # type: ignore[arg-type]
        points.append(tuple(base[i] + scale * rng.uniform(-1.0, 1.0) for i in range(3)))  # type: ignore[arg-type]
        faces.append((start, start + 1, start + 2))
    return points, faces


def build_meshes() -> dict[str, tuple[list[tuple[float, float, float]], list[tuple[int, int, int]]]]:
    return {
        "box_unit": box_surface(),
        "hull_surface_coarse": hull_surface(21, 9),
        "hull_surface_fine": hull_surface(41, 13),
        "sliver_patch": sliver_patch(),
        "random_soup_400": random_soup(),
    }


def max_rel_diff(got: list[float], expected: list[float]) -> dict[str, object]:
    worst, compared, excluded = 0.0, 0, 0
    for a, b in zip(got, expected):
        if not (math.isfinite(a) and math.isfinite(b)):
            excluded += 1
            continue
        compared += 1
        worst = max(worst, abs(a - b) / max(abs(a), abs(b), 1e-300))
    return {"max_relative_difference": worst, "compared": compared, "excluded_non_finite": excluded}


def main() -> int:
    ARTIFACT_DIR.mkdir(parents=True, exist_ok=True)
    meshes = build_meshes()
    failures: list[str] = []
    per_mesh: dict[str, object] = {}
    curves: dict[str, dict[str, list[float]]] = {}

    for name, (points, faces) in meshes.items():
        values = D.aspect_dialects(points, faces)
        # 1. identity against the three production modules
        # the legacy dialect must reproduce what the production report already says
        from tensorcad.geometry import Mesh

        reported = measure_quality(Mesh(list(points), list(faces))).max_aspect_ratio
        gor = [triangle_quality(points[a], points[b], points[c]).aspect_ratio for a, b, c in faces]
        new = list(tri_metrics(points, faces)["aspect_ratio"])
        identity = {
            "legacy_vs_validation_measure_quality": {
                "reported_max": reported,
                "dialect_max": max(values["legacy"]),
                "max_relative_difference": abs(reported - max(values["legacy"])) / max(reported, max(values["legacy"]), 1e-300),
            },
            "edge_form_vs_geomquality": max_rel_diff(list(values["edge_form"]), gor),
            "verdict_vs_meshquality": max_rel_diff(list(values["verdict"]), new),
        }
        if any(item["max_relative_difference"] > TOL for item in identity.values()):
            failures.append(f"{name}: identity beyond {TOL}")

        # 2. bounds
        alpha = [a for a in values["alpha"] if math.isfinite(a)]
        beta = [b for b in values["beta"] if math.isfinite(b)]
        bound_ok = (
            (not alpha or (min(alpha) >= D.ALPHA_MIN - 1e-12 and max(alpha) <= D.ALPHA_MAX + 1e-12))
            and (not beta or min(beta) >= D.BETA_MIN - 1e-12)
        )
        if not bound_ok:
            failures.append(f"{name}: bound violated")

        # 3. consequence
        flips = {}
        for dialect in ("edge_form", "verdict"):
            safe = D.conservative_limit(LIMIT, target=dialect)
            naive = sum(1 for g, r in zip(values[dialect], values["legacy"]) if g <= LIMIT and r > LIMIT)
            guarded = sum(1 for g, r in zip(values[dialect], values["legacy"]) if g <= safe and r > LIMIT)
            flips[dialect] = {"at_copied_limit_25": naive, "at_conservative_limit": guarded, "conservative_limit": safe}
            if guarded != 0:
                failures.append(f"{name}: conservative limit {dialect} passed {guarded} elements")
        per_mesh[name] = {
            "triangles": len(faces),
            "identity": identity,
            "alpha_min": min(alpha) if alpha else None,
            "alpha_max": max(alpha) if alpha else None,
            "beta_min": min(beta) if beta else None,
            "max_legacy": max(values["legacy"]),
            "max_edge_form": max(values["edge_form"]),
            "max_verdict": max(values["verdict"]),
            "flips": flips,
        }
        curves[name] = {
            "legacy": [v for v in values["legacy"] if math.isfinite(v)],
            "edge_form": [v for v in values["edge_form"] if math.isfinite(v)],
            "verdict": [v for v in values["verdict"] if math.isfinite(v)],
            "alpha": alpha,
        }

    # the two canonical modules, on the element where they disagree most visibly
    right = [(0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (1.0, 1.0, 0.0)]
    counterexample = {
        "element": right,
        "dialects": {k: v[0] for k, v in D.aspect_dialects(right, [(0, 1, 2)]).items()},
        "geomquality_aspect_ratio": triangle_quality(*right).aspect_ratio,
        "meshquality_aspect_ratio": tri_metrics(right, [(0, 1, 2)])["aspect_ratio"][0],
        "auditor_counterexample": {
            "element": [(0.0, 0.0, 0.0), (2.0, 0.0, 0.0), (1.0, 0.0741, 0.0)],
            **{k: v[0] for k, v in D.aspect_dialects(
                [(0.0, 0.0, 0.0), (2.0, 0.0, 0.0), (1.0, 0.0741, 0.0)], [(0, 1, 2)]
            ).items()},
        },
    }

    audit = {
        "schema": "tensorcad.meshquality.dialects.audit/1",
        "limit": LIMIT,
        "bounds": D.transfer_ratio_bounds(),
        "conservative_limits_from_legacy_25": {
            "edge_form": D.conservative_limit(LIMIT, target="edge_form"),
            "verdict": D.conservative_limit(LIMIT, target="verdict"),
        },
        "dialects": {k: dict(v) for k, v in D.ASPECT_DIALECTS.items()},
        "meshes": per_mesh,
        "disagreement_between_canonical_modules": counterexample,
        "failures": failures,
        "verdict": "PASS" if not failures else "FAIL",
    }
    (ARTIFACT_DIR / "audit.json").write_text(json.dumps(audit, indent=2, sort_keys=True))

    figure = plot(curves, counterexample)
    (ARTIFACT_DIR / "manifest.json").write_text(
        json.dumps(
            {
                "audit_sha256": hashlib.sha256((ARTIFACT_DIR / "audit.json").read_bytes()).hexdigest(),
                "figure": str(figure.name),
                "figure_bytes": figure.stat().st_size,
                "repro_cmd": "python tools/aspect_dialect_audit.py",
            },
            indent=2,
            sort_keys=True,
        )
    )
    print(json.dumps({"verdict": audit["verdict"], "failures": failures, "figure": str(figure)}, indent=2))
    return 0 if not failures else 1


def plot(curves: dict[str, dict[str, list[float]]], counterexample: dict[str, object]) -> Path:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(1, 3, figsize=(15.0, 4.6))
    colors = {"box_unit": "#1f77b4", "hull_surface_coarse": "#2ca02c", "hull_surface_fine": "#d62728",
              "sliver_patch": "#9467bd", "random_soup_400": "#7f7f7f"}

    ax = axes[0]
    for name, data in curves.items():
        n = min(len(data["legacy"]), len(data["edge_form"]))
        ax.scatter(data["legacy"][:n], data["edge_form"][:n], s=6, alpha=0.5, color=colors.get(name), label=name)
    xs = [10.0 ** (k / 20.0) for k in range(-20, 81)]
    ax.plot(xs, [x * D.ALPHA_MIN for x in xs], "k--", lw=1.0, label=r"$\alpha_{\min}=1/\sqrt{3}$")
    ax.plot(xs, [x * D.ALPHA_MAX for x in xs], "k:", lw=1.0, label=r"$\alpha_{\max}=\sqrt{3}/2$")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("legacy dialect  $l_{max}^2/2A$")
    ax.set_ylabel("edge_form dialect  $l_{max}p/(4\\sqrt{3}A)$")
    ax.set_title("(a) legacy vs edge_form: the band the proof predicts")
    ax.legend(fontsize=7, loc="upper left")
    ax.grid(alpha=0.3)

    ax = axes[1]
    for name, data in curves.items():
        n = min(len(data["legacy"]), len(data["verdict"]))
        ax.scatter(data["legacy"][:n], data["verdict"][:n], s=6, alpha=0.5, color=colors.get(name))
    ax.plot(xs, [0.5 * x for x in xs], "k--", lw=1.0, label=r"$q/legacy \geq 1/2$")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("legacy dialect  $l_{max}^2/2A$")
    ax.set_ylabel("verdict dialect  $R/(2r)$")
    ax.set_title("(b) legacy vs verdict (the two canonical modules' number)")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)

    ax = axes[2]
    alpha_all = [a for data in curves.values() for a in data["alpha"]]
    ax.hist(alpha_all, bins=60, color="#4c72b0", alpha=0.85)
    ax.axvline(D.ALPHA_MIN, color="k", ls="--", lw=1.0, label=r"$1/\sqrt{3}$")
    ax.axvline(D.ALPHA_MAX, color="k", ls=":", lw=1.0, label=r"$\sqrt{3}/2$")
    ex = counterexample["auditor_counterexample"]  # type: ignore[index]
    ax.axvline(ex["alpha"], color="#d62728", lw=1.0, label=f"F1 counterexample alpha={ex['alpha']:.4f}")
    ax.set_xlabel(r"$\alpha = edge_form/legacy$")
    ax.set_ylabel("triangles")
    ax.set_title("(c) transfer factor, %d triangles" % len(alpha_all))
    ax.legend(fontsize=7)
    ax.grid(alpha=0.3)

    fig.suptitle("TensorCAD mesh-quality aspect-ratio dialects: exact algebra and a provable threshold transfer")
    fig.tight_layout()
    out = ARTIFACT_DIR / "dialects.png"
    fig.savefig(out, dpi=140)
    plt.close(fig)
    return out


if __name__ == "__main__":
    sys.exit(main())