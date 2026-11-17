import json, os, urllib.request, urllib.error
API = os.environ.get("TENSOR_API_BASE", "http://127.0.0.1:8010/api/v1")
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
st, out = call("GET", "/tensor/missions/24137", tok=tok)
print(json.dumps(out, ensure_ascii=False, indent=1)[:9000])
