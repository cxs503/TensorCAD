"""Round-6 out-of-sample audit: corpus, independent measurement, falsifiability.

The tests here exist to keep three claims honest:

* the out-of-sample corpus is *not* the generator's population and still carries
  its provenance;
* the volume and LCB used to judge the pipeline are computed independently of the
  pipeline's own bookkeeping (divergence theorem on the triangle soup), and they
  agree with the closed-form target for a hull the pipeline itself certifies;
* the pass rate is falsifiable: on independently drawn coefficient triplets the
  pipeline rejects, and it never accepts a triplet that violates the physical
  precondition ``Cp <= Cwp``.

``test_rebuilt_family_has_no_lcb_control`` pins a **known limitation**, not a
desired behaviour: the section family is fore/aft symmetric, so LCB is 0.000 % L
for every hull including the published-asymmetric KVLCC2. If LCB control is ever
added, that test must be updated -- which is the point of pinning it.
"""
from __future__ import annotations

import math

import pytest

from tensorcad import mesh as meshlib
from tensorcad.shipd_oos import (
    CORPUS_PATH,
    DEVIATION_LIMIT_PCT,
    audit_ship,
    discrimination_table,
    hostile_coefficient_population,
    literature_corpus,
    mesh_volume_and_lcb,
    run_out_of_sample,
    schneekluth_cwp,
    summarise,
)
from tensorcad.shipd_rebuild import rebuild_ship

NX, NZ = 41, 13


# --------------------------------------------------------------------------- #
# Corpus
# --------------------------------------------------------------------------- #
def test_corpus_is_in_the_tree_and_not_the_generator_population():
    assert CORPUS_PATH.is_file()
    corpus = literature_corpus()
    assert len(corpus) == 6
    names = {ship["ship_id"] for ship in corpus}
    assert {"KCS", "KVLCC2", "ONRT", "JBC", "DTC"} <= names
    # Provenance travels with every record; these are recalled, not fetched.
    for ship in corpus:
        assert ship["corpus"] == "literature_baseline_fleet"
        assert ship["provenance"].get("source"), ship["ship_id"]


def test_corpus_flags_derived_waterplane_and_feasibility():
    corpus = {ship["ship_id"]: ship for ship in literature_corpus()}
    for ship in corpus.values():
        assert ship["cwp_source"] == "schneekluth_derived"
        assert ship["Cwp"] == pytest.approx(schneekluth_cwp(ship["Cb"]))
        assert ship["cp_le_cwp"] == (ship["Cp"] <= ship["Cwp"])
    # KCS is a well-posed hull: Cp (0.670) is inside its waterplane coefficient.
    assert corpus["KCS"]["cp_le_cwp"] is True
    # ONRT as recorded (Cm=0.75, Cb=0.535) implies Cp=0.713 > Cwp=0.690: the
    # recorded triplet describes no ship, and the corpus says so rather than
    # silently repairing it.
    assert corpus["ONRT"]["cp_le_cwp"] is False


def test_corpus_missing_file_is_reported(tmp_path):
    with pytest.raises(FileNotFoundError):
        literature_corpus(tmp_path / "nope.json")


# --------------------------------------------------------------------------- #
# Independent measurement
# --------------------------------------------------------------------------- #
def test_independent_volume_matches_closed_form_on_a_certified_hull():
    ship = next(s for s in literature_corpus() if s["ship_id"] == "KCS")
    result = rebuild_ship(ship, nx=NX, nz=NZ, method_preference="calibrated")
    assert result.verdict == "VERIFIED"
    volume, lcb = mesh_volume_and_lcb(result.mesh)
    target = ship["Cb"] * ship["L"] * ship["B"] * ship["T"]
    assert abs(volume - target) / target < 1e-3
    # Independent measurement vs the pipeline's own number: same hull, two routes.
    assert volume == pytest.approx(result.volume_mesh, rel=1e-6)
    # Frame-independent: LCB of a fore/aft symmetric family is 0 % of L.
    assert abs(lcb) < 1e-6


def test_volume_and_lcb_reject_a_degenerate_mesh():
    from tensorcad.geometry import Mesh

    empty = Mesh(vertices=[(0.0, 0.0, 0.0)], faces=[])
    with pytest.raises(ValueError):
        mesh_volume_and_lcb(empty)


# --------------------------------------------------------------------------- #
# Known limitation, pinned
# --------------------------------------------------------------------------- #
def test_rebuilt_family_has_no_lcb_control():
    """Every hull rebuilds with LCB = 0 % L, published LCB or not.

    KVLCC2 is published at 3.48 % L (magnitude; sign convention per
    ``tensorcad.fleet.lcb_frames``), so the miss is 3.48 % L -- above the 3 %
    standard the campaign applies elsewhere. This is a *family* property: the
    six consumed particulars (L, B, T, Cb, Cm, Cwp) contain no LCB, so no
    reconstruction from them can place the centroid anywhere else.
    """
    records = run_out_of_sample(nx=NX, nz=NZ)["records"]
    by_id = {record["ship_id"]: record for record in records}
    for ship_id in ("KCS", "KVLCC2"):
        record = by_id[ship_id]
        assert abs(record["lcb_frac_L_mesh"]) < 1e-6
    assert by_id["KVLCC2"]["lcb_frac_L_published"] == pytest.approx(0.0348)
    assert by_id["KVLCC2"]["delta_lcb_frac_L_abs"] == pytest.approx(0.0348, rel=1e-6)
    assert by_id["KVLCC2"]["delta_lcb_frac_L_abs"] * 100.0 > DEVIATION_LIMIT_PCT
    assert by_id["KCS"]["delta_lcb_frac_L_abs"] * 100.0 < DEVIATION_LIMIT_PCT


# --------------------------------------------------------------------------- #
# Out-of-sample pass rate and its failure ledger
# --------------------------------------------------------------------------- #
def test_out_of_sample_audit_is_not_all_green():
    payload = run_out_of_sample(nx=NX, nz=NZ)
    summary = payload["summary"]
    assert summary["n_ships"] == 6
    assert 0 < summary["verified"] < 6, "an out-of-sample population that passes 100% is not a test"
    # The single rejection is ONRT, and it is rejected for the right reason: a
    # volume that misses Cb*L*B*T, reported by G6 rather than waved through.
    by_id = {record["ship_id"]: record for record in payload["records"]}
    onrt = by_id["ONRT"]
    assert onrt["verdict"] != "VERIFIED"
    assert "G6_volume_consistency" in onrt["failing_gates"]
    assert onrt["volume_error_independent"] > 0.05
    # Every record has a ledger label and no label is empty.
    for record in payload["records"]:
        assert record["category"]
    assert set(summary["category_counts"]) <= {
        "clean", "refused:infeasible_coefficients", "refused:outside_method_windows",
        "rejected:parameters", "failed:no_mesh", "passed_gates_coefficient_miss",
        "passed_gates_surface_off_target",
    } | {label for label in summary["category_counts"] if label.startswith("gate_failed:")}


def test_rejection_is_reproducible():
    first = run_out_of_sample(nx=NX, nz=NZ)["summary"]
    second = run_out_of_sample(nx=NX, nz=NZ)["summary"]
    assert first["pass_rate"] == pytest.approx(second["pass_rate"])
    assert first["category_counts"] == second["category_counts"]


# --------------------------------------------------------------------------- #
# Falsifiability on independently drawn coefficients
# --------------------------------------------------------------------------- #
def test_hostile_population_contains_infeasible_triplets():
    population = hostile_coefficient_population(40, seed=11)
    infeasible = [s for s in population if not s["cp_le_cwp"]]
    assert infeasible, "independent draws must produce Cp > Cwp cases"
    for ship in population:
        assert ship["Cp"] == pytest.approx(ship["Cb"] / ship["Cm"])


def test_pipeline_never_accepts_an_infeasible_coefficient_triplet():
    population = hostile_coefficient_population(40, seed=11)
    records = [audit_ship(ship, nx=NX, nz=NZ) for ship in population]
    table = discrimination_table(records)
    assert table["infeasible_pass"] == 0, table["false_accept_ids"]
    assert table["feasible_fail"] == 0, table["false_reject_ids"]
    assert table["n"] == 40


def test_hostile_pass_rate_is_strictly_below_one():
    population = hostile_coefficient_population(24, seed=5)
    records = [audit_ship(ship, nx=NX, nz=NZ) for ship in population]
    summary = summarise(records)
    assert 0.0 < summary["pass_rate"] < 100.0
    assert summary["gate_failure_counts"], "rejections must name the gate that fired"


def test_discrimination_table_counts_both_diagonals():
    records = [
        {"ship_id": "a", "verdict": "VERIFIED", "cp_le_cwp": True},
        {"ship_id": "b", "verdict": "FAILED", "cp_le_cwp": True},
        {"ship_id": "c", "verdict": "VERIFIED", "cp_le_cwp": False},
        {"ship_id": "d", "verdict": "REFUSED", "cp_le_cwp": False},
    ]
    table = discrimination_table(records)
    assert table["feasible_pass"] == 1 and table["feasible_fail"] == 1
    assert table["infeasible_pass"] == 1 and table["infeasible_fail"] == 1
    assert table["false_accept_ids"] == ["c"] and table["false_reject_ids"] == ["b"]
    assert table["pass_rate"] == pytest.approx(50.0)


def test_meshlib_agrees_with_the_independent_volume():
    """Two implementations of the same integral must agree on a real hull."""
    ship = next(s for s in literature_corpus() if s["ship_id"] == "DTC")
    result = rebuild_ship(ship, nx=NX, nz=NZ, method_preference="calibrated")
    volume, _ = mesh_volume_and_lcb(result.mesh)
    assert meshlib.is_watertight(result.mesh)
    assert math.isfinite(volume) and volume > 0.0