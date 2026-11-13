#!/usr/bin/env python3
"""Probe m17317 mission result + reports + related verdicts."""
import json, subprocess

def psql(sql):
    r = subprocess.run(["runuser", "-u", "postgres", "--", "psql", "-d", "tensorplatform",
                        "-tA", "-F", "\t", "-q", "-v", "ON_ERROR_STOP=1", "-c", sql],
                       capture_output=True, text=True)
    if r.returncode:
        return f"ERR: {r.stderr.strip()[:300]}"
    return r.stdout.strip()

print("=== m17317 full result (json) ===")
print(psql("select result from tensor_dev_missions where id=17317;"))
print()
print("=== m17317 changed_files ===")
print(psql("select jsonb_pretty(result->'changed_files') from tensor_dev_missions where id=17317;"))
print()
print("=== m17317 pr_gate ===")
print(psql("select result->'pr_gate' from tensor_dev_missions where id=17317;"))
print()
print("=== m17252 mission (original) ===")
print(psql("select id, repo, status, created_at::text, left(goal,300) from tensor_dev_missions where id=17252;"))
print()
print("=== m17252 result ===")
print(psql("select result from tensor_dev_missions where id=17252;"))
print()
print("=== all verdicts for GAP-CAD-T2-PARAM-SERIES scope ===")
print(psql("select repo, scope, mission_id, level, gate_level, gate_state, verdict_by, created_at::text "
           "from tensor_maturity_verdicts where scope like '%GAP-CAD-T2-PARAM-SERIES%' order by created_at;"))
print()
print("=== reports for mission 17317 / 17252 (tensor_reports) ===")
print(psql("select id, mission_id, left(title,80), created_at::text from tensor_reports "
           "where mission_id in (17317,17252) order by created_at;"))
