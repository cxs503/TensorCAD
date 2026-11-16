"""Pin KVLCC2 ``lwl`` as a documented single-lineage gap (anchor REPRO-GEOM-FLEET-REFS).

The cross-source library's most impactful remaining cell is KVLCC2 ``lwl`` = 325.5 m:
three *documents* state it but all three are the same SIMMAN-SERIES *lineage*, so
criterion ① (≥2 independent lineages) fails while ①b (≥2 documents) passes.  These tests
pin that state so a later edit cannot quietly reclassify the three SIMMAN documents as two
independent sources, and so the two genuinely independent lineages already in the library
are shown *not* to contribute an ``lwl`` cell (the gap is real, not a parse miss).
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from tensorcad.fleet.crossref import lineage_of

LIBRARY = Path("benchmarks/verified/fleet_refs/crossref.json")
ACCEPTANCE = Path("benchmarks/verified/fleet_refs/acceptance.json")


@pytest.fixture(scope="module")
def library() -> dict:
    return json.loads(LIBRARY.read_text())


@pytest.fixture(scope="module")
def acceptance() -> dict:
    return json.loads(ACCEPTANCE.read_text())


def test_kvlcc2_lwl_has_three_documents_but_one_lineage(library):
    """①b (≥2 documents) passes; ① (≥2 independent lineages) fails — the honest gap."""
    row = library["report"]["KVLCC2"]["quantities"]["lwl"]
    documents = {entry["cell"]["source"] for entry in row["entries"]}
    assert len(documents) >= 2, "①b must already pass: ≥2 documents state Lwl"
    assert {"simman2020_summary", "simman2014_kvlcc2", "kvlcc2_osaka"} <= documents
    # ... but every document is the same lineage, so ① must fail.
    lineage = library["lineage"]["KVLCC2"]["lwl"]
    assert lineage["verdict"] == "SINGLE_SOURCE", lineage
    assert lineage["lineages"] == ["SIMMAN-SERIES"], lineage


def test_all_three_lwl_documents_collapse_to_one_lineage(library):
    """A second paper that reprints the same workshop table is not a second source."""
    for source in ("simman2020_summary", "simman2014_kvlcc2", "kvlcc2_osaka"):
        assert lineage_of(source) == "SIMMAN-SERIES", source


def test_independent_lineages_do_not_contribute_lwl(library):
    """The two genuinely independent families state Lpp/B/T/∇/Cb but not Lwl — real gap."""
    for source in ("kvlcc2_ugent", "kvlcc2_jstage"):
        assert lineage_of(source) != "SIMMAN-SERIES", source
    row = library["report"]["KVLCC2"]["quantities"]["lwl"]
    contributing = {entry["cell"]["source"] for entry in row["entries"]}
    assert "kvlcc2_ugent" not in contributing
    assert "kvlcc2_jstage" not in contributing


def test_criterion_one_remains_honest_fail_for_kvlcc2_lwl(acceptance):
    """Criterion ① must stay an honest FAIL for KVLCC2-``lwl``; do not reclassify."""
    one = acceptance["criteria"]["one_two_independent_lineages"]
    assert one["pass"] is False
    kvlcc2 = one["ships"]["KVLCC2"]
    assert "lwl" in kvlcc2["failing"]
