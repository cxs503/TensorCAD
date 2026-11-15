"""Post work/m21856_report.json to the platform reports endpoint."""
from __future__ import annotations

import json
import os
import urllib.request
from pathlib import Path

ROOT = Path("/root/tp-agent/work/m21856")
BASE = os.environ["TENSOR_API_BASE"]
USER = os.environ["TENSOR_AGENT_USER"]
PASS = os.environ["TENSOR_AGENT_PASS"]


def post(url: str, payload: dict, token: str | None = None) -> dict:
    data = json.dumps(payload).encode()
    req = urllib.request.Request(url, data=data, method="POST",
                                 headers={"Content-Type": "application/json"})
    if token:
        req.add_header("Authorization", f"Bearer {token}")
    with urllib.request.urlopen(req, timeout=300) as resp:
        return json.loads(resp.read().decode())


login = post(f"{BASE}/auth/login", {"username": USER, "password": PASS})
token = login.get("access_token") or login.get("token")
if not token:
    raise SystemExit(f"login failed: {list(login)}")

payload = json.loads((ROOT / "work" / "m21856_report.json").read_text(encoding="utf-8"))
res = post(f"{BASE}/tensor/reports", payload, token)
print(json.dumps(res, ensure_ascii=False)[:800])