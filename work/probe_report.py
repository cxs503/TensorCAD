#!/usr/bin/env python3
"""Verify m17317 report + gate power + audit diff."""
import json, subprocess

def psql(sql):
    r = subprocess.run(["runuser", "-u", "postgres", "--", "psql", "-d", "tensorplatform",
                        "-tA", "-q", "-v", "ON_ERROR_STOP=1", "-c", sql],
                       capture_output=True, text=True)
    if r.returncode:
        return f"ERR: {r.stderr.strip()[:300]}"
    return r.stdout.strip()

print("=== report 16314 (m17317) ===")
print(psql("select id, mission_id, title, summary, scenario_key, tags, created_at::text from tensor_reports where id=16314;"))
print()
print("=== report 16314 markdown (first 2000) ===")
print(psql("select left(markdown, 2000) from tensor_reports where id=16314;"))
