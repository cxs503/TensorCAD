import json, os, urllib.request, urllib.error
API = os.environ.get("TENSOR_API_BASE", "http://127.0.0.1:8010/api/v1")
import os
USER = os.environ.get("TENSOR_AGENT_USER", "eco_admin")
PASS = os.environ.get("TENSOR_AGENT_PASS", "")
def call(method, path, body=None, tok=None):
    data = json.dumps(body).encode() if body is not None else None
    headers = {"Content-Type": "application/json"}
    if tok:
        headers["Authorization"] = "Bearer " + tok
    r = urllib.request.Request(API + path, data=data, headers=headers, method=method)
    try:
        with urllib.request.urlopen(r) as resp:
            return resp.status, json.loads(resp.read() or b"null")
    except urllib.error.HTTPError as e:
        return e.code, json.loads(e.read() or b"null")
st, out = call("POST", "/auth/login", {"username": USER, "password": PASS})
tok = out.get("access_token", "") if isinstance(out, dict) else ""
for p in ["/missions?limit=5", "/tensor/missions/24137", "/missions/24137/status", "/tensor/mission/24137", "/api/v1/missions/24137"]:
    st, out = call("GET", p, tok=tok)
    print(p, "->", st, str(out)[:300])
