import json, os, urllib.request

BASE = os.environ["TENSOR_API_BASE"].rstrip("/")
def req(path, body=None, token=None, method=None):
    data = json.dumps(body).encode() if body is not None else None
    r = urllib.request.Request(BASE + path, data=data, method=method or ("POST" if data else "GET"))
    r.add_header("Content-Type", "application/json")
    if token:
        r.add_header("Authorization", "Bearer " + token)
    with urllib.request.urlopen(r, timeout=60) as f:
        return json.load(f)

tok = req("/auth/login", {"username": os.environ["TENSOR_AGENT_USER"],
                          "password": os.environ["TENSOR_AGENT_PASS"]})
tok = tok.get("access_token") or tok.get("token")
open("/tmp/tok.txt", "w").write(tok)
print("login ok", len(tok))

rs = req("/tensor/reports?limit=25", token=tok)
rs = rs if isinstance(rs, list) else (rs.get("reports") or rs.get("items") or rs.get("data") or [])
print("N reports:", len(rs))
for r in rs:
    print(r.get("id"), "|", r.get("mission_id"), "|", str(r.get("title"))[:90])