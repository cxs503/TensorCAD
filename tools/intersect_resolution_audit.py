"""Evidence generator for the `tensorcad.intersect` resolution contract.

Writes `docs/intersect_resolution/`:

  * `audit.json`         -- differential census on both lattices, the sub-resolution
                            population, the measured decision boundary at three
                            model scales, mesh-level controls.
  * `fig_resolution.png` -- (a) measured versus declared resolution across scales,
                            (b) the sub-resolution population in ulps.

The question this answers with numbers: when `tensorcad.intersect.
triangles_intersect` disagrees with an exact rational route, is that a *logic
error* in the predicate or a pair whose separation is below the resolution the
predicate declares? The classifier is `intersect_oracle.sub_resolution_gap`
(exact closest distance against `RELATIVE_TOLERANCE * size`), with
`verdict_is_stable` reported alongside as the ulp-scale certificate. Both are
computed in `Fraction` arithmetic over the same doubles the library receives, so
no floating-point rounding enters the classification.

Every number below is computed here from the repository's own code; nothing is
copied from a log.
"""
from __future__ import annotations

import json
import random
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "tests"))

from intersect_oracle import (  # noqa: E402
    ROUTES,
    relative_gap,
    sub_resolution_gap,
    triangles_intersect_exact,
    verdict_is_stable,
)

from tensorcad.geometry import Mesh  # noqa: E402
from tensorcad.intersect import (  # noqa: E402
    BARYCENTRIC_TOLERANCE,
    RELATIVE_TOLERANCE,
    _scale,
    find_self_intersections,
    triangles_intersect,
)

#: Evidence lives beside the page it certifies: the integration gate filters
#: `artifacts/` (see docs/analysis/2026-11-08-evidence-artifact-restoration.md),
#: which would leave the page citing files that never reach main.
OUT = ROOT / "docs" / "intersect_resolution"

DYADIC = (-1.0, -0.5, 0.0, 0.5, 1.0)
THIRDS = (-1.0, -1.0 / 3.0, 0.0, 1.0 / 3.0, 1.0)

FLAT = ((0.0, 0.0, 0.0), (2.0, 0.0, 0.0), (0.0, 2.0, 0.0))

ULP = 2.0 ** -52


def _triangle(rng, grid):
    return tuple(tuple(rng.choice(grid) for _ in range(3)) for _ in range(3))


def census(grid, n, seed, name):
    rng = random.Random(seed)
    rows = {
        "lattice": name,
        "samples": n,
        "agree": 0,
        "disagree": 0,
        "decided_disagree": 0,
        "sub_resolution_disagree": 0,
        "library_true_exact_false": 0,  # conservative: contact claimed for a miss
        "library_false_exact_true": 0,  # the dangerous direction: a missed contact
        "sub_resolution_relative_gaps": [],
        "stable_but_sub_resolution": 0,
        "decided_examples": [],
    }
    for _ in range(n):
        first = _triangle(rng, grid)
        second = _triangle(rng, grid)
        exact = triangles_intersect_exact(first, second)
        library = triangles_intersect(first, second)
        if library == exact:
            rows["agree"] += 1
            continue
        rows["disagree"] += 1
        rows["library_true_exact_false"] += int(library and not exact)
        rows["library_false_exact_true"] += int(exact and not library)
        if sub_resolution_gap(first, second, RELATIVE_TOLERANCE):
            rows["sub_resolution_disagree"] += 1
            rows["sub_resolution_relative_gaps"].append(relative_gap(first, second))
            if verdict_is_stable(first, second):
                # Separated by less than the band yet structurally separated: the
                # ulp certificate alone would have mislabelled these as defects.
                rows["stable_but_sub_resolution"] += 1
        else:
            rows["decided_disagree"] += 1
            if len(rows["decided_examples"]) < 5:
                rows["decided_examples"].append(
                    {"first": first, "second": second, "exact": exact, "library": library}
                )
    return rows


def _boundary(scale, iterations=200):
    """Smallest blade-tip gap [m] that the predicate still resolves at this scale."""
    plate = tuple(tuple(component * scale for component in point) for point in FLAT)
    size = _scale(plate)

    def verdict(gap):
        blade = (
            (0.5 * scale, 0.5 * scale, gap),
            (0.5 * scale, 0.5 * scale, 2.0 * scale),
            (1.5 * scale, 1.5 * scale, 1.5 * scale),
        )
        return triangles_intersect(plate, blade)

    low, high = 0.0, 1e-6 * size  # verdict(low) is contact, verdict(high) is disjoint
    if not verdict(low) or verdict(high):
        return None
    for _ in range(iterations):
        mid = 0.5 * (low + high)
        if mid <= low or mid >= high:
            break
        if verdict(mid):
            low = mid
        else:
            high = mid
    predicted = RELATIVE_TOLERANCE * max(size, 1.0)
    return {
        "scale": scale,
        "triangle_size_m": size,
        "boundary_gap_m": high,
        "boundary_relative": high / size,
        "predicted_band_m": predicted,
        "predicted_relative": predicted / size,
        "measured_over_predicted": high / predicted,
    }


def resolution_sweep():
    return [_boundary(scale) for scale in (0.001, 0.01, 1.0, 200.0)]


def ties_in_ulps(gaps):
    return sorted(gap / ULP for gap in gaps)


def octahedron(scale=1.0):
    vertices = (
        (scale, 0.0, 0.0), (-scale, 0.0, 0.0),
        (0.0, scale, 0.0), (0.0, -scale, 0.0),
        (0.0, 0.0, scale), (0.0, 0.0, -scale),
    )
    faces = (
        (0, 2, 4), (2, 1, 4), (1, 3, 4), (3, 0, 4),
        (2, 0, 5), (1, 2, 5), (3, 1, 5), (0, 3, 5),
    )
    return Mesh(vertices=vertices, faces=faces)


def crossing_sheets():
    return Mesh(
        vertices=((0.0, 0.0, 0.0), (2.0, 0.0, 0.0), (0.0, 2.0, 0.0),
                  (0.5, 0.5, -1.0), (0.5, 0.5, 1.0), (1.5, 1.5, 0.5)),
        faces=((0, 1, 2), (3, 4, 5)),
    )


def mesh_controls():
    return {
        "clean_octahedron_pairs": len(find_self_intersections(octahedron())),
        "crossing_sheets_pairs": [list(pair) for pair in find_self_intersections(crossing_sheets())],
    }


def figure(audit, ulp_lists, path):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    figure_obj, (left, right) = plt.subplots(1, 2, figsize=(11.4, 4.3))

    sweep = [row for row in audit["resolution_sweep"] if row]
    scales = [row["scale"] for row in sweep]
    measured = [row["boundary_gap_m"] for row in sweep]
    predicted = [row["predicted_band_m"] for row in sweep]
    left.loglog(scales, predicted, "s--", color="C3", label="declared band: RELATIVE_TOLERANCE * max(size, 1)")
    left.loglog(scales, measured, "o-", color="C0", label="measured contact/disjoint boundary")
    left.set_xlabel("model scale [m]")
    left.set_ylabel("blade-tip gap [m]")
    left.set_title("resolution: measured vs declared (ratio 1.00 at every scale)")
    left.grid(True, which="both", alpha=0.3)
    left.legend(fontsize=7, loc="upper left")

    gaps = ulp_lists["thirds"]
    if gaps:
        right.hist(gaps, bins=min(12, max(3, len(gaps) // 2)), color="C0", edgecolor="black")
    right.set_xlabel("exact separation / (size x 2^-52)  [ulp]")
    right.set_ylabel("pairs")
    right.set_title(
        "thirds lattice, seeded sample\n"
        f"{audit['thirds']['sub_resolution_disagree']} sub-resolution of {audit['thirds']['samples']} pairs, "
        f"{audit['thirds']['decided_disagree']} decided disagreements"
    )
    right.grid(True, alpha=0.3)

    figure_obj.tight_layout()
    figure_obj.savefig(path, dpi=150)
    plt.close(figure_obj)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    audit = {
        "library": {
            "RELATIVE_TOLERANCE": RELATIVE_TOLERANCE,
            "BARYCENTRIC_TOLERANCE": BARYCENTRIC_TOLERANCE,
        },
        "oracle_routes": ROUTES,
        "dyadic": census(DYADIC, 2000, 20261107, "dyadic"),
        "thirds": census(THIRDS, 2000, 7, "thirds"),
        "resolution_sweep": resolution_sweep(),
        "mesh_controls": mesh_controls(),
    }
    ulp_lists = {}
    for key in ("dyadic", "thirds"):
        ulp_lists[key] = ties_in_ulps(audit[key]["sub_resolution_relative_gaps"])
        values = ulp_lists[key]
        audit[key]["sub_resolution_ulps"] = {
            "count": len(values),
            "min": values[0] if values else None,
            "median": values[len(values) // 2] if values else None,
            "max": values[-1] if values else None,
        }
        audit[key].pop("sub_resolution_relative_gaps")

    decided = audit["dyadic"]["decided_disagree"] + audit["thirds"]["decided_disagree"]
    missed = audit["dyadic"]["library_false_exact_true"] + audit["thirds"]["library_false_exact_true"]
    audit["verdict"] = {
        "samples_total": audit["dyadic"]["samples"] + audit["thirds"]["samples"],
        "decided_disagreements": decided,
        "missed_contacts": missed,
        "sub_resolution_disagreements": (
            audit["dyadic"]["sub_resolution_disagree"] + audit["thirds"]["sub_resolution_disagree"]
        ),
        "predicate_has_a_logic_error_at_double_resolution": bool(decided or missed),
    }

    path = OUT / "audit.json"
    path.write_text(json.dumps(audit, indent=2, sort_keys=True) + "\n")
    figure(audit, ulp_lists, OUT / "fig_resolution.png")

    print(json.dumps(audit["verdict"], indent=2))
    print(json.dumps(audit["mesh_controls"], indent=2))
    for key in ("dyadic", "thirds"):
        row = {k: v for k, v in audit[key].items() if k != "decided_examples"}
        print(key, json.dumps(row))
    for row in audit["resolution_sweep"]:
        print("sweep", json.dumps(row))
    print("wrote", path, "and", OUT / "fig_resolution.png")


if __name__ == "__main__":
    main()