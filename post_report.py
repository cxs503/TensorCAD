#!/usr/bin/env python3
"""Post the acceptance report to the TensorPlatform."""
import json
import os
import urllib.request
import urllib.error

API_BASE = os.environ.get("TENSOR_API_BASE", "http://127.0.0.1:8010/api/v1")
USER = os.environ.get("TENSOR_AGENT_USER", "eco_admin")
PASS = os.environ.get("TENSOR_AGENT_PASS", "")

def login():
    data = json.dumps({"username": USER, "password": PASS}).encode()
    req = urllib.request.Request(f"{API_BASE}/auth/login", data=data, headers={"Content-Type": "application/json"})
    with urllib.request.urlopen(req) as resp:
        return json.loads(resp.read())["access_token"]

def post_report(token, report_path="/root/tp-agent/work/m25806/runs/report.json"):
    with open(report_path) as f:
        report = json.load(f)
    data = json.dumps(report).encode()
    req = urllib.request.Request(f"{API_BASE}/tensor/reports", data=data,
                                  headers={"Content-Type": "application/json",
                                           "Authorization": f"Bearer {token}"})
    try:
        with urllib.request.urlopen(req) as resp:
            result = json.loads(resp.read())
            print(f"Report posted: {result}")
            return result
    except urllib.error.HTTPError as e:
        body = e.read().decode()
        print(f"HTTP {e.code}: {body}")
        return None

if __name__ == "__main__":
    token = login()
    print(f"Logged in, token: {token[:20]}...")
    post_report(token)
