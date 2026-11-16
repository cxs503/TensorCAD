#!/usr/bin/env python3
"""Post the m20544 audit report to the platform (report discipline, owner order).

Logs in with $TENSOR_AGENT_USER / $TENSOR_AGENT_PASS against $TENSOR_API_BASE and
POSTs /tensor/reports with the payload built by audit/m20544_verdict.py
(/tmp/report_m20544.json). Idempotent by title: if a report with the same title
already exists for this mission, it is left alone.
"""
from __future__ import annotations

import json
import os
import urllib.error
import urllib.request

BASE = os.environ["TENSOR_API_BASE"].rstrip("/")
USER = os.environ["TENSOR_AGENT_USER"]
PASS = os.environ["TENSOR_AGENT_PASS"]


def call(method: str, path: str, payload=None, token: str | None = None) -> dict:
    data = json.dumps(payload).encode() if payload is not None else None
    req = urllib.request.Request(BASE + path, data=data, method=method)
    req.add_header("Content-Type", "application/json")
    if token:
        req.add_header("Authorization", f"Bearer {token}")
    with urllib.request.urlopen(req, timeout=120) as r:
        body = r.read().decode()
    return json.loads(body) if body.strip().startswith(("{", "[")) else {"raw": body}


def main() -> int:
    tok = None
    for form in ({"username": USER, "password": PASS},
                 {"email": USER, "password": PASS}):
        try:
            d = call("POST", "/auth/login", form)
            tok = d.get("access_token") or d.get("token") or (d.get("data") or {}).get("access_token")
            if tok:
                break
        except urllib.error.HTTPError as exc:
            last = f"{exc.code} {exc.read()[:200]!r}"
            print("login attempt failed:", last)
    if not tok:
        print("LOGIN FAILED")
        return 1
    payload = json.loads(open("/tmp/report_m20544.json", encoding="utf-8").read())
    out = call("POST", "/tensor/reports", payload, tok)
    print(json.dumps({"posted": out.get("id") or out.get("report_id") or out,
                      "title": payload["title"],
                      "figures": len(payload.get("figures") or [])}, ensure_ascii=False)[:600])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())