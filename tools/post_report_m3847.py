"""Publish the m3847 (PR #44 cure) progress report to the platform.

Builds the payload from the committed artifacts — markdown from
``docs/m3847_pr44_merge_cure.md``, figure from ``artifacts/m3847/gz_integration.png`` —
and POSTs it to ``$TENSOR_API_BASE/tensor/reports`` after logging in as the agent user.
The same payload is left in ``artifacts/m3847/report_payload.json`` as evidence.

    python3 tools/post_report_m3847.py
"""

from __future__ import annotations

import base64
import json
import os
import urllib.error
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
BASE = os.environ["TENSOR_API_BASE"].rstrip("/")
USER = os.environ["TENSOR_AGENT_USER"]
PASSWORD = os.environ["TENSOR_AGENT_PASS"]

CAPTION = (
    "PR #44 cure cross-feature evidence. (a) branch GZ solver vs the classical "
    "wall-sided closed form on a 2 x 0.8 m box with 1.2 m freeboard; the dotted line "
    "is the validity limit tan(phi)=2T/B. (b) the merge seam: the branch solver driven "
    "on main's transom hull (tau=0.20) with the branch's freeboard topsides, and on "
    "main's Wigley with 0.3 m topsides. (c) log-scale relative deviation of (a): "
    "<=1.11e-08 inside the window, 5.66 % one step outside it (asserted departure, not "
    "tolerance). (d) test counts of the cured tree vs the main baseline under the same "
    "pytest command: 879 vs 818 passed, 2 vs 2 failed (pre-existing), 11 vs 11 skipped."
)

SUMMARY = (
    "PR #44 (agent/m3347, old head 9f1572f) was dirty against main 2ff80b2. Merged main "
    "into the branch with a true 3-way merge: 1 textual conflict (CHANGELOG 0.21.0 "
    "section) resolved as a verbatim union of both sides' bullets; 3 auto-merged shared "
    "files verified symbol-by-symbol against both parents (125 exports = 112 union 120, "
    "none lost); 37 single-sided files byte-identical to the side that changed them; "
    "single_sided_content_lost = []. Merge commit 5fe2912. Tests on the cured tree: 879 "
    "passed / 2 failed / 11 skipped vs main baseline 818 / 2 / 11 (+61 = the branch's GZ "
    "suite; the same two pre-existing failures). Cross-feature evidence: the branch GZ "
    "solver on main's transom hull and Wigley (GM agreement 0.0, displacement conserved "
    "exactly over 0-90 deg), box GZ vs wall-sided closed form 1.11e-08 inside the "
    "validity window, transom hydrostatics vs closed form 0.785 % at the reference grid. "
    "Push to origin/agent/m3347 was NOT possible from this sandbox (GitHub unreachable, "
    "git push blocked by task policy): the cure ships as commit 5fe2912 plus the "
    "incremental bundle artifacts/m3847/agent_m3347_cured.bundle (prerequisite main "
    "2ff80b2)."
)


def post(path: str, body: dict, token: str | None = None) -> dict:
    headers = {"Content-Type": "application/json"}
    if token:
        headers["Authorization"] = f"Bearer {token}"
    request = urllib.request.Request(
        BASE + path, data=json.dumps(body).encode(), headers=headers, method="POST"
    )
    with urllib.request.urlopen(request, timeout=120) as response:
        return json.loads(response.read().decode() or "{}")


def login() -> str:
    body = {"username": USER, "password": PASSWORD}
    last: Exception | None = None
    for path in (
        "/tensor/auth/login",
        "/tensor/login",
        "/auth/login",
        "/api/login",
    ):
        try:
            reply = post(path, body)
        except urllib.error.HTTPError as exc:  # wrong endpoint: try the next guess
            last = exc
            continue
        for key in ("token", "access_token", "jwt"):
            if reply.get(key):
                return reply[key]
        data = reply.get("data") or {}
        for key in ("token", "access_token", "jwt"):
            if data.get(key):
                return data[key]
        last = RuntimeError(f"no token in reply from {path}: {reply}")
    raise SystemExit(f"login failed: {last}")


def main() -> int:
    payload = {
        "title": "PR #44 cure (m3847): semantic union of agent/m3347 and main in TensorCAD",
        "summary": SUMMARY,
        "scenario_key": "tensorcad/pr44-merge-cure",
        "mission_id": 3847,
        "tags": [
            "TensorCAD",
            "merge-cure",
            "PR44",
            "semantic-union",
            "GZ-curve",
            "transom-stern",
            "m3847",
        ],
        "markdown": (ROOT / "docs" / "m3847_pr44_merge_cure.md").read_text(),
        "figures": [
            {
                "name": "gz_integration.png",
                "caption": CAPTION,
                "png_base64": base64.b64encode(
                    (ROOT / "artifacts" / "m3847" / "gz_integration.png").read_bytes()
                ).decode(),
            }
        ],
    }
    (ROOT / "artifacts" / "m3847" / "report_payload.json").write_text(
        json.dumps(payload, indent=2) + "\n"
    )
    token = login()
    reply = post("/tensor/reports", payload, token)
    print("POST /tensor/reports ->", json.dumps(reply)[:400])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())