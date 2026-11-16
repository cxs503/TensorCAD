"""Out-of-sample audit of the ShipD rebuild pipeline (round 6).

Why this module exists
======================
Rounds 1--5 of the ``GAP-CAD-C2-SHIPD-REBUILD`` campaign report a pass rate over
500 ships and, from round 5, prove the gate instrument is *sharp* on mesh-level
defects (an injected inverted triangle makes ``G4`` fire).

Both statements are about **the generator**. ``tensorcad.shipd_population``
samples ``Cb`` and then derives ``Cm``/``Cwp`` from textbook relations with a
small scatter, and raises ``Cwp`` whenever the resulting ``Cp = Cb/Cm`` would
violate ``Cp <= Cwp``. The 500-ship table therefore contains no infeasible
coefficient triplet by construction, and the coefficient-calibrated family was
written to realise exactly the triplet it is handed. "100 % pass" is then a
statement about a closed loop -- generator, method and gate written by the same
hand, tuned until they agree -- and it carries no information about *other*
ships. A pass rate is evidence only when the population can fail.

This module supplies the two populations that can:

1. **Out of sample** (:func:`literature_corpus`): the six benchmark hulls whose
   principal particulars are already in the tree
   (``tensorcad/fleet/data/fleet_manifest.json``: KCS, KVLCC2, ONRT, JBC, DTC,
   Series 60). They are *not* products of the ShipD generator: they are the
   published workshop hulls a naval architect actually reasons about, and their
   coefficients were not chosen to make this pipeline look good. Small ``n``
   (6) -- which is why every record is reported individually rather than as a
   mean.

2. **Deliberately hostile** (:func:`hostile_coefficient_population`): draws
   ``Cb``, ``Cm`` and ``Cwp`` *independently* from plausible ranges, so a share
   of the draws violate ``Cp <= Cwp`` -- a waterplane cannot be finer than the
   prismatic coefficient, so those parameter vectors describe no ship at all.
   This turns the pass rate into a contingency table: does the pipeline accept
   the feasible draws and reject the impossible ones?

A third thing the campaign never measured is measured here: the **longitudinal
centre of buoyancy**. ``LCB`` is a principal particular of every real hull and
one of the 45 components of the MIT ShipD input vector, but it is absent from the
six scalars this pipeline consumes (``L, B, T, Cb, Cm, Cwp``). Every rebuilt
hull is section-symmetric fore and aft, so its ``LCB`` is 0.000 % ``L`` by
construction -- a property of the family, not of any hull. :func:`run_out_of_sample`
reports that number next to the published one instead of leaving it unmeasured.

Independence of the measurements
================================
The gate verdicts come from the code under audit (:func:`tensorcad.shipd_rebuild.rebuild_ship`).
The two quantities used to *judge* those verdicts are computed here, from the
triangle soup, by the divergence theorem -- :func:`mesh_volume_and_lcb` -- and not
taken from ``RebuildResult``, so a bug in the pipeline's own volume bookkeeping
cannot hide behind its own gate. The only thing imported from the audited side is
the surface itself.

Honest provenance and limits
============================
* No network and no local ShipD payload (``data/raw/shipd/`` is empty), so the
  MIT 30k distribution still cannot be sampled; :func:`literature_corpus` is the
  best available out-of-sample surrogate and it is six hulls, not 500.
* The corpus particulars are ``retrieval=offline_recollection, verified=false``
  in the manifest. They are recall, not a fetched data sheet.
* ``Cwp`` is not in the manifest for any hull; it is derived from
  :func:`schneekluth_cwp` and flagged ``cwp_source="schneekluth_derived"``. The
  ternary then decides which method the pipeline may use, so a hull whose real
  ``Cwp`` is finer than the regression implies could be rejected for a reason
  that is ours, not the hull's.
* Published ``LCB`` sign conventions differ (see
  :mod:`tensorcad.fleet.lcb_frames`). The manifest is read in the canonical frame
  documented there -- percent of ``Lpp``, positive forward, from the ``Lpp``
  mid-point -- and the comparison is reported as a **magnitude** ``|dLCB|`` so a
  frame disagreement cannot be mistaken for a small (signed) miss.
"""
from __future__ import annotations

import json
import math
from pathlib import Path
from typing import Any, Sequence

import numpy as np

from tensorcad.geometry import Mesh
from tensorcad.shipd_population import (
    holtrop_mennen_wetted_surface,
    nominal_performance_fields,
)
from tensorcad.shipd_rebuild import rebuild_ship

#: In-tree baseline fleet manifest (literature particulars + provenance).
CORPUS_PATH = Path(__file__).resolve().parent / "fleet" / "data" / "fleet_manifest.json"

#: Deviation threshold the campaign's own acceptance standard uses.
DEVIATION_LIMIT_PCT = 3.0

#: Ranges the hostile population draws from, per coefficient, independently.
HOSTILE_CB_RANGE = (0.40, 0.88)
HOSTILE_CM_RANGE = (0.85, 0.995)
HOSTILE_CWP_RANGE = (0.62, 0.92)


# --------------------------------------------------------------------------- #
# Corpus
# --------------------------------------------------------------------------- #
def schneekluth_cwp(cb: float) -> float:
    """Schneekluth & Bertram waterplane regression ``Cwp = (1 + 2*Cb)/3``."""
    return (1.0 + 2.0 * cb) / 3.0


def literature_corpus(path: str | Path | None = None) -> list[dict[str, Any]]:
    """Benchmark hulls from the in-tree fleet manifest, as rebuild input vectors.

    Each record keeps the provenance block it came with, adds the derived
    ``Cp = Cb/Cm`` and the feasibility flag the pipeline's method selection uses,
    and records where ``Cwp`` came from.
    """
    source = Path(path) if path is not None else CORPUS_PATH
    payload = json.loads(source.read_text())
    ships: list[dict[str, Any]] = []
    for record in payload["ships"]:
        cb = float(record["block_coefficient"])
        cm = float(record.get("midship_coefficient") or schneekluth_cwp(cb))
        cwp_published = record.get("waterplane_coefficient")
        if cwp_published is None:
            cwp = schneekluth_cwp(cb)
            cwp_source = "schneekluth_derived"
        else:
            cwp = float(cwp_published)
            cwp_source = "published"
        lcb = record.get("lcb_fraction")
        ships.append({
            "ship_id": record["name"],
            "ship_type": record.get("ship_class", record["name"]),
            "L": float(record["length"]),
            "B": float(record["beam"]),
            "T": float(record["draft"]),
            "Cb": cb,
            "Cm": cm,
            "Cwp": cwp,
            "cwp_source": cwp_source,
            "Cp": cb / cm,
            "cp_le_cwp": cb / cm <= cwp,
            "lcb_frac_L_published": None if lcb is None else float(lcb),
            "design_froude": record.get("design_froude"),
            "geometry_route": record.get("geometry_route"),
            "provenance": record.get("provenance", {}),
            "corpus": "literature_baseline_fleet",
        })
    return ships


def hostile_coefficient_population(n: int, *, seed: int = 20261116) -> list[dict[str, Any]]:
    """``n`` parameter vectors whose coefficients are drawn *independently*.

    The ShipD generator ties ``Cm`` and ``Cwp`` to ``Cb`` and repairs a tripped
    ``Cp <= Cwp`` violation before the pipeline ever sees it. Here the three
    coefficients are independent uniforms over plausible ranges, so the sample
    contains coefficient triplets no ship can have -- which is what makes a pass
    rate falsifiable.
    """
    if n <= 0:
        raise ValueError("n must be positive")
    rng = np.random.default_rng(seed)
    ships: list[dict[str, Any]] = []
    for index in range(n):
        cb = round(float(rng.uniform(*HOSTILE_CB_RANGE)), 5)
        cm = round(float(rng.uniform(*HOSTILE_CM_RANGE)), 5)
        cwp = round(float(rng.uniform(*HOSTILE_CWP_RANGE)), 5)
        ships.append({
            "ship_id": f"hostile_{index + 1:03d}",
            "ship_type": "hostile_draw",
            "L": float(rng.uniform(100.0, 350.0)),
            "B": float(rng.uniform(20.0, 55.0)),
            "T": float(rng.uniform(8.0, 20.0)),
            "Cb": cb,
            "Cm": cm,
            "Cwp": cwp,
            "cwp_source": "independent_draw",
            "Cp": cb / cm,
            "cp_le_cwp": cb / cm <= cwp,
            "lcb_frac_L_published": None,
            "design_froude": 0.20,
            "corpus": "hostile_independent_coefficients",
        })
    return ships


# --------------------------------------------------------------------------- #
# Independent mesh measurement
# --------------------------------------------------------------------------- #
def mesh_volume_and_lcb(mesh: Mesh) -> tuple[float, float]:
    """Signed volume [m^3] and LCB [fraction of L, +forward from mid-length].

    Divergence theorem on the closed triangle soup: each triangle ``(a, b, c)``
    contributes the signed tetrahedron ``(0, a, b, c)``, of volume
    ``a . (b x c) / 6`` and centroid ``(a + b + c) / 4``. ``LCB`` is the volume
    centroid's ``x`` referred to the mesh's own longitudinal mid-point, so the
    result does not depend on whether the generator put ``x = 0`` at the AP or at
    mid-length.
    """
    vertices = np.asarray(mesh.vertices, dtype=float)
    faces = np.asarray(mesh.faces, dtype=np.int64)
    if faces.size == 0:
        raise ValueError("mesh has no faces")
    a = vertices[faces[:, 0]]
    b = vertices[faces[:, 1]]
    c = vertices[faces[:, 2]]
    tet_volume = np.einsum("ij,ij->i", a, np.cross(b, c)) / 6.0
    volume = float(tet_volume.sum())
    if abs(volume) < 1e-12:
        raise ValueError("mesh encloses no volume")
    centroid_x = float((tet_volume * (a[:, 0] + b[:, 0] + c[:, 0]) / 4.0).sum()) / volume
    x_min = float(vertices[:, 0].min())
    x_max = float(vertices[:, 0].max())
    length = x_max - x_min
    if length <= 0.0:
        raise ValueError("mesh has zero longitudinal extent")
    return volume, (centroid_x - 0.5 * (x_min + x_max)) / length


# --------------------------------------------------------------------------- #
# One ship
# --------------------------------------------------------------------------- #
def audit_ship(ship: dict[str, Any], *, nx: int = 41, nz: int = 13,
               method_preference: str = "calibrated") -> dict[str, Any]:
    """Rebuild one ship and measure the result independently of its gates."""
    result = rebuild_ship(ship, nx=nx, nz=nz, method_preference=method_preference)
    failing = [gate.name for gate in result.gates if not gate.passed]

    volume_mesh = lcb_frac_L = float("nan")
    measured = False
    if result.mesh is not None:
        try:
            volume_mesh, lcb_frac_L = mesh_volume_and_lcb(result.mesh)
            measured = True
        except ValueError:
            measured = False

    target_volume = ship["Cb"] * ship["L"] * ship["B"] * ship["T"]
    volume_error_independent = (
        abs(volume_mesh - target_volume) / target_volume if measured and target_volume > 0 else float("nan")
    )

    published_lcb = ship.get("lcb_frac_L_published")
    delta_lcb = (
        abs(lcb_frac_L - published_lcb) if measured and published_lcb is not None else float("nan")
    )

    wetted_holtrop = holtrop_mennen_wetted_surface(ship)
    wetted_mesh = result.wetted_surface_area
    wetted_deviation_pct = (
        100.0 * (wetted_mesh - wetted_holtrop) / wetted_holtrop
        if wetted_mesh and wetted_holtrop > 0 else float("nan")
    )

    # `nominal_performance_fields` needs a design Froude number; the literature
    # hulls do not all publish one, and `dict.get(k, default)` does not replace an
    # explicit `None`, so the default is applied here and flagged in the record.
    # `Cp` and its feasibility flag are derived here rather than trusted from the
    # record: in-sample rows come straight from the generator and carry neither.
    cp = ship["Cb"] / ship["Cm"]
    cp_le_cwp = bool(ship.get("cp_le_cwp", cp <= ship["Cwp"]))

    nominal_input = dict(ship)
    froude_defaulted = nominal_input.get("design_froude") is None
    if froude_defaulted:
        nominal_input["design_froude"] = 0.20
    nominal = nominal_performance_fields(nominal_input)
    record = {
        "ship_id": ship["ship_id"],
        "ship_type": ship.get("ship_type"),
        "corpus": ship.get("corpus", "unspecified"),
        "L": ship["L"], "B": ship["B"], "T": ship["T"],
        "Cb": ship["Cb"], "Cm": ship["Cm"], "Cwp": ship["Cwp"],
        "Cp": cp,
        "cp_le_cwp": cp_le_cwp,
        "cwp_source": ship.get("cwp_source"),
        "method": result.method,
        "method_reason": result.method_reason,
        "verdict": result.verdict,
        "strict_verdict": result.strict_verdict,
        "failing_gates": failing,
        "error_message": result.error_message,
        "volume_error_reported": result.volume_error,
        "volume_error_independent": volume_error_independent,
        "volume_mesh_independent": volume_mesh,
        "volume_target": target_volume,
        "lcb_frac_L_mesh": lcb_frac_L,
        "lcb_frac_L_published": published_lcb,
        "delta_lcb_frac_L_abs": delta_lcb,
        "wetted_surface_mesh": wetted_mesh,
        "wetted_surface_holtrop": wetted_holtrop,
        "wetted_deviation_pct": wetted_deviation_pct,
        "design_froude_used": nominal_input["design_froude"],
        "design_froude_defaulted": froude_defaulted,
        "nominal_viscous_resistance_kN": nominal["viscous_resistance_kN"],
        "triangle_count": result.triangle_count,
        "lcb_frac_L_published_sign_note": (
            "magnitude only: manifest read in the canonical frame of "
            "tensorcad.fleet.lcb_frames (percent Lpp, +forward, from mid)"
        ),
    }
    record["category"] = classify_failure(record)
    return record


# --------------------------------------------------------------------------- #
# Failure-mode ledger
# --------------------------------------------------------------------------- #
def classify_failure(record: dict[str, Any]) -> str:
    """One label per ship; the ledger the acceptance criterion asks for.

    Order matters: a geometric rejection outranks a metric miss, because a hull
    that never became a mesh has no metric to be off target.
    """
    verdict = record["verdict"]
    if verdict == "REFUSED":
        if not record["cp_le_cwp"]:
            return "refused:infeasible_coefficients"
        return "refused:outside_method_windows"
    if verdict != "VERIFIED":
        if record["failing_gates"]:
            return "gate_failed:" + "+".join(record["failing_gates"])
        if record["error_message"]:
            # No mesh, no gate, but a stated reason: a parent-series or parameter
            # check refused the hull before geometry existed (route seen on
            # Cb > 0.80 hulls the Series 60 ladder does not cover).
            return "rejected:parameters"
        return "failed:no_mesh"
    if record["strict_verdict"] not in (None, "VERIFIED"):
        return "passed_gates_coefficient_miss"
    wetted = record["wetted_deviation_pct"]
    if not math.isnan(wetted) and abs(wetted) > DEVIATION_LIMIT_PCT:
        return "passed_gates_surface_off_target"
    return "clean"


def discrimination_table(records: Sequence[dict[str, Any]]) -> dict[str, Any]:
    """2x2 contingency: coefficient feasibility x pipeline acceptance.

    ``feasible`` means ``Cp <= Cwp`` -- a physical precondition, decided from the
    parameter vector alone, before any geometry exists. The two off-diagonal
    cells are the instrument's two error modes and both must be zero for the
    pass rate to mean anything.
    """
    table = {"feasible_pass": 0, "feasible_fail": 0, "infeasible_pass": 0, "infeasible_fail": 0}
    offenders: dict[str, list[str]] = {"false_accept": [], "false_reject": []}
    for record in records:
        passed = record["verdict"] == "VERIFIED"
        if record["cp_le_cwp"]:
            if passed:
                table["feasible_pass"] += 1
            else:
                table["feasible_fail"] += 1
                offenders["false_reject"].append(record["ship_id"])
        else:
            if passed:
                table["infeasible_pass"] += 1
                offenders["false_accept"].append(record["ship_id"])
            else:
                table["infeasible_fail"] += 1
    n = len(records)
    table["n"] = n
    table["pass_rate"] = (100.0 * (table["feasible_pass"] + table["infeasible_pass"]) / n) if n else float("nan")
    table["false_accept_ids"] = offenders["false_accept"]
    table["false_reject_ids"] = offenders["false_reject"]
    return table


# --------------------------------------------------------------------------- #
# Population runner
# --------------------------------------------------------------------------- #
def _finite(values: Sequence[float]) -> list[float]:
    return [v for v in values if isinstance(v, (int, float)) and not math.isnan(v)]


def summarise(records: Sequence[dict[str, Any]]) -> dict[str, Any]:
    """Aggregate ledger + the numbers a reader needs to judge the pipeline."""
    n = len(records)
    counts: dict[str, int] = {}
    for record in records:
        counts[record["category"]] = counts.get(record["category"], 0) + 1
    gate_failures: dict[str, int] = {}
    for record in records:
        for gate in record["failing_gates"]:
            gate_failures[gate] = gate_failures.get(gate, 0) + 1

    volume_errors = _finite([r["volume_error_independent"] for r in records])
    lcb_mesh = _finite([r["lcb_frac_L_mesh"] for r in records])
    delta_lcb = _finite([r["delta_lcb_frac_L_abs"] for r in records])
    wetted = _finite([r["wetted_deviation_pct"] for r in records])

    def stats(values: list[float]) -> dict[str, float]:
        if not values:
            return {"n": 0}
        return {
            "n": len(values),
            "min": min(values),
            "max": max(values),
            "mean": sum(values) / len(values),
            "mean_abs": sum(abs(v) for v in values) / len(values),
        }

    verified = sum(1 for r in records if r["verdict"] == "VERIFIED")
    return {
        "n_ships": n,
        "verified": verified,
        "pass_rate": (100.0 * verified / n) if n else float("nan"),
        "category_counts": counts,
        "gate_failure_counts": gate_failures,
        "volume_error_independent": stats(volume_errors),
        "lcb_frac_L_mesh": stats(lcb_mesh),
        "delta_lcb_frac_L_abs": stats(delta_lcb),
        "wetted_deviation_pct": stats(wetted),
        "lcb_over_limit_ids": [
            r["ship_id"] for r in records
            if not math.isnan(r["delta_lcb_frac_L_abs"]) and r["delta_lcb_frac_L_abs"] * 100.0 > DEVIATION_LIMIT_PCT
        ],
        "wetted_over_limit_ids": [
            r["ship_id"] for r in records
            if not math.isnan(r["wetted_deviation_pct"]) and abs(r["wetted_deviation_pct"]) > DEVIATION_LIMIT_PCT
        ],
    }


def run_out_of_sample(ships: Sequence[dict[str, Any]] | None = None, *, nx: int = 41,
                      nz: int = 13, method_preference: str = "calibrated") -> dict[str, Any]:
    """Audit the literature corpus (default) or any supplied parameter vectors."""
    population = list(ships) if ships is not None else literature_corpus()
    records = [
        audit_ship(ship, nx=nx, nz=nz, method_preference=method_preference)
        for ship in population
    ]
    return {
        "schema": "tensorcad.shipd_oos.v1",
        "anchor": "GAP-CAD-C2-SHIPD-REBUILD",
        "grid": {"nx": nx, "nz": nz},
        "method_preference": method_preference,
        "summary": summarise(records),
        "discrimination": discrimination_table(records),
        "records": records,
    }


def run_hostile_population(n: int = 120, *, seed: int = 20261116, nx: int = 41,
                           nz: int = 13, method_preference: str = "calibrated") -> dict[str, Any]:
    """Same audit over ``n`` independently drawn coefficient triplets."""
    population = hostile_coefficient_population(n, seed=seed)
    payload = run_out_of_sample(population, nx=nx, nz=nz, method_preference=method_preference)
    payload["population"] = {"kind": "hostile_independent_coefficients", "n": n, "seed": seed,
                             "cb_range": list(HOSTILE_CB_RANGE), "cm_range": list(HOSTILE_CM_RANGE),
                             "cwp_range": list(HOSTILE_CWP_RANGE)}
    return payload


def write_json(payload: dict[str, Any], path: str | Path) -> Path:
    """Write ``payload`` to ``path`` (parent directories created)."""
    destination = Path(path)
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_text(json.dumps(payload, indent=2, sort_keys=False))
    return destination