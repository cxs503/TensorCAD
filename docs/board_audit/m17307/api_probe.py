"""Audit probe: read m17251 artifacts, benchmark cases, reports via API/DB."""
import json, subprocess, urllib.request, urllib.parse, os

BASE = "http://127.0.0.1:8010/api/v1"
login = json.load(open("/tmp/login.json"))
tk = login["access_token"]

def api(path, method="GET", body=None):
    url = BASE + "/" + path.lstrip("/")
    data = json.dumps(body).encode() if body is not None else None
    req = urllib.request.Request(url, data=data, method=method)
    req.add_header("Authorization", "Bearer " + tk)
    if data:
        req.add_header("Content-Type", "application/json")
    try:
        with urllib.request.urlopen(req, timeout=60) as r:
            return r.status, json.loads(r.read().decode())
    except Exception as e:
        err = getattr(e, "read", lambda: b"")()
        return getattr(e, "code", -1), err[:2000].decode(errors="replace")

if __name__ == "__main__":
    import sys
    for p in sys.argv[1:]:
        st, d = api(p)
        print("===", p, st)
        print(json.dumps(d, ensure_ascii=False)[:4000] if not isinstance(d, str) else d[:4000])