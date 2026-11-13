#!/usr/bin/env python3
"""Probe the TensorPlatform DB for the m17317 cross-audit subject."""
import json, subprocess, os, sys

def psql(sql):
    r = subprocess.run(["runuser", "-u", "postgres", "--", "psql", "-d", "tensorplatform",
                        "-tA", "-F", "\t", "-q", "-v", "ON_ERROR_STOP=1", "-c", sql],
                       capture_output=True, text=True)
    if r.returncode:
        return f"ERR: {r.stderr.strip()[:300]}"
    return r.stdout.strip()

print("=== mission m17317 ===")
print(psql("select id, repo, status, created_at::text, left(goal, 400) from tensor_dev_missions where id=17317;"))
print()
print("=== verdicts for scope containing m17317 / GAP-CAD-T2-PARAM-SERIES ===")
print(psql("select repo, scope, mission_id, level, gate_level, gate_state, verdict_by, created_at::text "
           "from tensor_maturity_verdicts where scope like '%GAP-CAD-T2-PARAM-SERIES%' or mission_id=17317 "
           "order by created_at;"))
print()
print("=== benchmark cases for anchor GAP-CAD-T2-PARAM-SERIES ===")
print(psql("select case_name, metric, tolerance, anchor_slug from tensor_benchmark_cases "
           "where anchor_slug='GAP-CAD-T2-PARAM-SERIES';"))
print()
print("=== all benchmark cases (count + sample) ===")
print(psql("select count(*) from tensor_benchmark_cases;"))
print(psql("select case_name, metric, tolerance, anchor_slug from tensor_benchmark_cases limit 20;"))
print()
print("=== tensorbench_results for mission 17317 ===")
print(psql("select count(*) from tensorbench_results where mission_id=17317;"))
print()
print("=== feature_ledger for GAP-CAD-T2-PARAM-SERIES ===")
print(psql("select slug, domain, status, title from feature_ledger where slug like '%PARAM-SERIES%' or title like '%Series 60%';"))
