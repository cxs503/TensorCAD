#!/usr/bin/env python3
"""Publish the PC6 report payload to the platform report endpoint."""
from __future__ import annotations

import json
import os
import sys
import urllib.error
import urllib.request

BASE = os.environ.get("TENSOR_API_BASE", "http://127.0.0.1:8010/api/v1").rstrip("/")
USER = os.environ.get("TENSOR_AGENT_USER", "")
PASSWORD = os.environ.get("TENSOR_AGENT_PASS", "")
PAYLOAD = "/tmp/pc6_report_payload.json"
MISSION_ID = "4929"


def post(path: str, body: dict, token: str | None = None) -> tuple[int, str]:
    headers = {"Content-Type": "application/json"}
    if token:
        headers["Authorization"] = f"Bearer {token}"
    request = urllib.request.Request(
        BASE + path, data=json.dumps(body).encode(), headers=headers, method="POST"
    )
    try:
        with urllib.request.urlopen(request, timeout=90) as response:
            return response.status, response.read().decode("utf-8", "replace")
    except urllib.error.HTTPError as error:
        return error.code, error.read().decode("utf-8", "replace")


def login() -> str:
    for path in ("/auth/login", "/tensor/auth/login", "/login"):
        status, text = post(path, {"username": USER, "password": PASSWORD})
        if status == 200:
            data = json.loads(text)
            token = data.get("access_token") or data.get("token") or ""
            if token:
                print(f"login via {path} ok")
                return str(token)
        print(f"login via {path} -> {status} {text[:120]}")
    raise SystemExit("could not log in")


def main() -> int:
    token = login()
    with open(PAYLOAD, encoding="utf-8") as handle:
        payload = json.load(handle)
    status, text = post("/tensor/reports", payload, token)
    print("POST /tensor/reports ->", status, text[:600])
    # read back, which is the only proof the write landed
    request = urllib.request.Request(
        f"{BASE}/tensor/reports?mission_id={MISSION_ID}",
        headers={"Authorization": f"Bearer {token}"},
        method="GET",
    )
    try:
        with urllib.request.urlopen(request, timeout=60) as response:
            body = response.read().decode("utf-8", "replace")
        data = json.loads(body)
        items = data if isinstance(data, list) else data.get("items", [])
        print(f"read back {len(items)} report(s) for mission {MISSION_ID}")
        for item in items[:2]:
            print(" title:", str(item.get("title"))[:110])
            print(" md5:", item.get("md5") or item.get("content_md5"))
    except urllib.error.HTTPError as error:
        print("read back failed", error.code, error.read()[:200])
    return 0 if status in (200, 201) else 1


if __name__ == "__main__":
    sys.exit(main())
