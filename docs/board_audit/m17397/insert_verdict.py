#!/usr/bin/env python3
"""Insert the board-audit maturity verdict for m17317 (subject) into
tensor_maturity_verdicts, then print the row back for the record.

Idempotent: one row per (repo, scope, verdict_by) — re-running updates it.
Run: python3 docs/board_audit/m17397/insert_verdict.py [--dry-run]
"""
from __future__ import annotations

import argparse
import json
import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
EVIDENCE = ROOT / "docs/board_audit/m17397/evidence.json"
REPO = "TensorCAD"
SCOPE = "m17317 GAP-CAD-T2-PARAM-SERIES"
MISSION = 17397          # this audit task
LEVEL = 2                # independent functional implementation; L3/L4 refused with reasons
BY = "board-audit"


def psql(sql: str) -> str:
    r = subprocess.run(["runuser", "-u", "postgres", "--", "psql", "-d", "tensorplatform",
                        "-tA", "-v", "ON_ERROR_STOP=1", "-c", sql],
                       capture_output=True, text=True)
    if r.returncode:
        raise SystemExit(f"psql failed: {r.stderr.strip()}")
    return r.stdout.strip()


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--dry-run", action="store_true")
    a = ap.parse_args()
    ev = json.loads(EVIDENCE.read_text())
    ev_text = json.dumps(ev, ensure_ascii=False).replace("'", "''")
    sql = (
        "INSERT INTO tensor_maturity_verdicts (repo, scope, mission_id, level, evidence, verdict_by) "
        f"VALUES ('{REPO}', '{SCOPE}', {MISSION}, {LEVEL}, '{ev_text}'::jsonb, '{BY}') "
        "ON CONFLICT (repo, scope, verdict_by) DO UPDATE SET "
        f"mission_id={MISSION}, level={LEVEL}, evidence='{ev_text}'::jsonb, gate_state='pending', "
        "gate_note=NULL, gate_level=NULL, created_at=now() RETURNING id, repo, scope, level, "
        "gate_state, gate_level;"
    )
    if a.dry_run:
        print(sql[:600], "...")
        return 0
    print(psql(sql))
    print(psql("select id, repo, scope, mission_id, level, verdict_by, gate_state, "
               "jsonb_object_keys(evidence) from tensor_maturity_verdicts "
               f"where repo='{REPO}' and scope='{SCOPE}' and verdict_by='{BY}';"))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
