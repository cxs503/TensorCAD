#!/usr/bin/env python3
"""Insert the board-audit maturity verdict for m17357 into tensor_maturity_verdicts.

Sandbox has no working psycopg driver, so the write goes through
`runuser -u postgres -- psql` (same path maturity_gate.py itself uses).
Idempotent on (repo, scope, verdict_by).

Before writing, the row is replayed through the real gate judge
(/root/tp-agent/maturity_gate.py:judge) so the claimed level and the machine
level are shown side by side — a claimed level the gate cannot verify would be
downgraded, and that must be visible in the log, not hidden.

Run: python3 docs/board_audit/m17405/insert_verdict.py [--dry-run]
"""
from __future__ import annotations

import argparse
import importlib.util
import json
import subprocess
from pathlib import Path

HERE = Path(__file__).resolve().parent
EVIDENCE = HERE / "evidence.json"
REPO = "TensorCAD"
SCOPE = "m17357 GAP-CAD-D3-GEOM-CONTRACT"
MISSION = 17405          # this audit task (not the subject)
LEVEL = 2                # independent functional implementation; L3/L4 refused
BY = "board-audit"


def psql(sql: str) -> str:
    r = subprocess.run(["runuser", "-u", "postgres", "--", "psql", "-d", "tensorplatform",
                        "-tA", "-v", "ON_ERROR_STOP=1", "-c", sql],
                       capture_output=True, text=True)
    if r.returncode:
        raise SystemExit(f"psql failed: {r.stderr.strip()}")
    return r.stdout.strip()


def gate_preview(ev: dict) -> tuple[int, list[str]]:
    """Call the real maturity-gate judge on the evidence bundle."""
    spec = importlib.util.spec_from_file_location("mg", "/root/tp-agent/maturity_gate.py")
    mg = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mg)          # module level only defines helpers
    return mg.judge(LEVEL, ev, REPO, MISSION)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--dry-run", action="store_true")
    a = ap.parse_args()
    ev = json.loads(EVIDENCE.read_text())

    gl, why = gate_preview(ev)
    print(f"gate preview: claimed L{LEVEL} -> gate_level L{gl} "
          f"({'confirmed' if gl >= 2 else 'REOPEN'})")
    for w in why:
        print("  ", w)

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
        print(sql[:400], "...")
        return 0
    print("row:", psql(sql))
    print(psql("select id, repo, scope, mission_id, level, verdict_by, gate_state, "
               "jsonb_array_length(coalesce(evidence->'checks_not_run','[]'::jsonb)) "
               "from tensor_maturity_verdicts "
               f"where repo='{REPO}' and scope='{SCOPE}' and verdict_by='{BY}';"))
    print("evidence keys:", psql("select string_agg(k, ',' order by k) from "
                                 f"(select jsonb_object_keys(evidence) k from tensor_maturity_verdicts "
                                 f"where repo='{REPO}' and scope='{SCOPE}' and verdict_by='{BY}') s;"))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())