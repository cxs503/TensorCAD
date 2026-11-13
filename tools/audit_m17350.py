#!/usr/bin/env python3
"""m17402 交差审计辅助脚本：登录 API + 查询基准/verdict 表 + 发布报告。"""
import json
import os
import sys
import urllib.request

BASE = os.environ["TENSOR_API_BASE"].rstrip("/")
USER = os.environ["TENSOR_AGENT_USER"]
PASS = os.environ["TENSOR_AGENT_PASS"]


def req(method, path, body=None, token=None):
    url = BASE + path
    data = json.dumps(body).encode() if body is not None else None
    r = urllib.request.Request(url, data=data, method=method)
    r.add_header("Content-Type", "application/json")
    if token:
        r.add_header("Authorization", "Bearer " + token)
    try:
        with urllib.request.urlopen(r, timeout=60) as resp:
            return resp.status, json.loads(resp.read().decode() or "{}")
    except urllib.error.HTTPError as e:
        return e.code, json.loads(e.read().decode() or "{}")


def login():
    st, d = req("POST", "/auth/login", {"username": USER, "password": PASS})
    if st != 200:
        sys.exit(f"login failed {st}: {d}")
    return d.get("access_token") or d.get("token")


if __name__ == "__main__":
    tok = login()
    cmd = sys.argv[1] if len(sys.argv) > 1 else "paths"
    if cmd == "paths":
        st, d = req("GET", "/openapi.json", token=tok)
        if st != 200:
            sys.exit(f"openapi {st}")
        for p in d.get("paths", {}):
            if any(k in p for k in ("maturity", "benchmark", "verdict", "report", "mission")):
                print(p)
    elif cmd == "bench":
        st, d = req("GET", "/tensor/benchmark-cases?anchor_slug=GAP-CAD-B2-MESH-FORMAT-MATRIX", token=tok)
        print(st, json.dumps(d, ensure_ascii=False)[:3000])
    elif cmd == "verdict":
        body = json.load(open(sys.argv[2]))
        st, d = req("POST", "/tensor/maturity-verdicts", body, token=tok)
        print(st, json.dumps(d, ensure_ascii=False)[:2000])
    elif cmd == "report":
        body = json.load(open(sys.argv[2]))
        st, d = req("POST", "/tensor/reports", body, token=tok)
        print(st, json.dumps(d, ensure_ascii=False)[:2000])
