"""Verify the DTC lcb_percent cross-contamination defect (anchor REPRO-GEOM-FLEET-REFS).

The cross-source library must not let a KVLCC2-tagged source (``kvlcc2_ugent``, JSR-2017
multi-ship Table 1) place its DTC design-draught LCG into DTC's ``lcb_percent`` row as if it
were a DTC document.  This test pins the defect so the next round starts from evidence.
"""

from __future__ import annotations

import json
from pathlib import Path

LIBRARY = Path("benchmarks/verified/fleet_refs/crossref.json")


def test_dtc_lcb_percent_not_sourced_from_kvlcc2():
    lib = json.loads(LIBRARY.read_text())
    row = lib["report"]["DTC"]["quantities"]["lcb_percent"]
    sources = {e["cell"]["source"] for e in row["entries"]}
    # The DTC lcb_percent row must not be fed by a KVLCC2-tagged source.
    # (Defect: it currently is -- kvlcc2_ugent's DTC design-draught LCG is misattributed.)
    assert "kvlcc2_ugent" not in sources, (
        "DTC lcb_percent is cross-contaminated by kvlcc2_ugent "
        "(JSR-2017 multi-ship table DTC design-draught LCG)"
    )


def test_dtc_lwl_is_missing():
    lib = json.loads(LIBRARY.read_text())
    row = lib["report"]["DTC"]["quantities"]["lwl"]
    assert row["verdict"] == "MISSING", f"DTC lwl should be MISSING, got {row['verdict']}"
    assert not row["entries"], "DTC lwl has no documents"
