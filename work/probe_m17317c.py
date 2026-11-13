#!/usr/bin/env python3
"""Dump full m17317 result JSON."""
import json, subprocess

def psql(sql):
    r = subprocess.run(["runuser", "-u", "postgres", "--", "psql", "-d", "tensorplatform",
                        "-tA", "-q", "-v", "ON_ERROR_STOP=1", "-c", sql],
                       capture_output=True, text=True)
    if r.returncode:
        return f"ERR: {r.stderr.strip()[:300]}"
    return r.stdout.strip()

raw = psql("select result from tensor_dev_missions where id=17317;")
try:
    d = json.loads(raw)
    print("=== top-level keys ===")
    print(list(d.keys()))
    print()
    print("=== changed_files ===")
    print(json.dumps(d.get("changed_files"), indent=1))
    print()
    print("=== on_disk_changes ===")
    print(json.dumps(d.get("on_disk_changes"), indent=1))
    print()
    print("=== evidence ===")
    print(json.dumps(d.get("evidence"), indent=1)[:3000])
    print()
    print("=== checks ===")
    print(json.dumps(d.get("checks"), indent=1))
    print()
    print("=== maturity_gate ===")
    print(json.dumps(d.get("maturity_gate"), indent=1))
    print()
    print("=== pr_gate ===")
    print(json.dumps(d.get("pr_gate"), indent=1))
    print()
    print("=== headline ===")
    print(d.get("headline"))
except Exception as e:
    print("parse err", e)
    print(raw[:2000])
