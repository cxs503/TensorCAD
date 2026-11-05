import json
import os
import urllib.request

base = os.environ["TENSOR_API_BASE"]
with open("/tmp/report.json") as handle:
    payload = json.load(handle)
body = json.dumps({"username": os.environ["TENSOR_AGENT_USER"],
                   "password": os.environ["TENSOR_AGENT_PASS"]}).encode()
req = urllib.request.Request(base + "/auth/login", data=body,
                             headers={"Content-Type": "application/json"})
resp = json.load(urllib.request.urlopen(req, timeout=60))
token = resp.get("access_token") or resp.get("token")
print("login keys:", sorted(resp)[:6], "token len", len(token or ""))
data = json.dumps(payload).encode()
req = urllib.request.Request(base + "/tensor/reports", data=data,
                             headers={"Content-Type": "application/json",
                                      "Authorization": "Bearer " + token})
try:
    out = json.load(urllib.request.urlopen(req, timeout=180))
except urllib.error.HTTPError as exc:
    print("HTTP", exc.code, exc.read()[:500].decode(errors="replace"))
else:
    print("report id:", out.get("id") or out.get("report_id") or str(out)[:300])
