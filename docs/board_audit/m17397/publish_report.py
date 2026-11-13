#!/usr/bin/env python
"""Publish the m17397 audit report to the platform report API."""
import json
import os
import urllib.request

base = os.environ.get("TENSOR_API_BASE", "http://127.0.0.1:8010/api/v1").rstrip("/")
user = os.environ["TENSOR_AGENT_USER"]
password = os.environ["TENSOR_AGENT_PASS"]


def post(path, payload, token=None):
    headers = {"Content-Type": "application/json"}
    if token:
        headers["Authorization"] = f"Bearer {token}"
    request = urllib.request.Request(base + path, data=json.dumps(payload).encode(), headers=headers)
    with urllib.request.urlopen(request, timeout=120) as response:
        return json.loads(response.read() or b"{}")


auth = post("/auth/login", {"username": user, "password": password})
token = auth.get("access_token") or auth.get("token")
body = json.load(open("/root/tp-agent/work/m17397/docs/board_audit/m17397/report.json"))
result = post("/tensor/reports", body, token)
print(json.dumps({k: result[k] for k in list(result)[:6]}, default=str))
