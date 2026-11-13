#!/usr/bin/env python3
"""Read-only audit helper: query tensor platform for benchmark cases / mission record.

Usage: python3 tools/audit_q.py cases <anchor_slug>
       python3 tools/audit_q.py mission <id>
       python3 tools/audit_q.py sql "<select ...>"
"""
import json
import os
import sys
import urllib.parse
import urllib.request

BASE = os.environ.get("TENSOR_API_BASE", "http://127.0.0.1:8010/api/v1")
USER = os.environ.get("TENSOR_AGENT_USER", "eco_admin")
PASS = os.environ.get("TENSOR_AGENT_PASS", "")


def token():
    req = urllib.request.Request(
        BASE.replace("/api/v1", "") + "/api/v1/auth/login",
        data=json.dumps({"username": USER, "password": PASS}).encode(),
        headers={"Content-Type": "application/json"},
        method="POST",
    )
    with urllib.request.urlopen(req, timeout=30) as r:
        d = json.load(r)
    return d.get("access_token") or d.get("token")


def get(path, tok):
    req = urllib.request.Request(BASE + path, headers={"Authorization": "Bearer " + tok})
    with urllib.request.urlopen(req, timeout=60) as r:
        return json.load(r)


def main():
    cmd = sys.argv[1]
    tok = token()
    if cmd == "cases":
        slug = urllib.parse.quote(sys.argv[2])
        print(json.dumps(get(f"/tensor/benchmark_cases?anchor_slug={slug}", tok),
                         ensure_ascii=False, indent=2)[:8000])
    elif cmd == "mission":
        print(json.dumps(get(f"/tensor/missions/{sys.argv[2]}", tok),
                         ensure_ascii=False, indent=2)[:12000])
    elif cmd == "get":
        print(json.dumps(get(sys.argv[2], tok), ensure_ascii=False, indent=2)[:8000])
    elif cmd == "sql":
        sql = sys.argv[2]
        url = (BASE + "/tensor/sql?q=" + urllib.parse.quote(sql))
        try:
            print(json.dumps(get(url[len(BASE):], tok), ensure_ascii=False, indent=2)[:10000])
        except Exception as e:
            print("sql endpoint failed:", e)
    else:
        print(__doc__)


if __name__ == "__main__":
    main()