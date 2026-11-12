import json, os, urllib.request

base = os.environ["TENSOR_API_BASE"].rstrip("/")


def post(path, body, token=None):
    req = urllib.request.Request(base + path, data=json.dumps(body).encode(), method="POST")
    req.add_header("Content-Type", "application/json")
    if token:
        req.add_header("Authorization", "Bearer " + token)
    with urllib.request.urlopen(req, timeout=60) as r:
        return json.load(r)


def get(path, token):
    req = urllib.request.Request(base + path)
    req.add_header("Authorization", "Bearer " + token)
    with urllib.request.urlopen(req, timeout=60) as r:
        return json.load(r)


d = post("/auth/login", {"username": os.environ["TENSOR_AGENT_USER"],
                         "password": os.environ["TENSOR_AGENT_PASS"]})
tok = d.get("access_token") or d.get("token") or d.get("data", {}).get("access_token")
open("/tmp/tok_m16690", "w").write(tok or "")
reps = get("/tensor/reports?limit=12", tok)
items = reps.get("items") or reps.get("reports") or reps
if isinstance(items, dict):
    items = items.get("items", [])
for r in items[:12]:
    print(r.get("id"), "|", (r.get("title") or "")[:70], "|", str(r.get("created_at"))[:10],
          "|", r.get("mission_id"))