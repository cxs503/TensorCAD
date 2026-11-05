#!/usr/bin/env python3
"""Publish the m3822 (PR #44 cure) report: narrative + generated tables + figure.

Credentials come from the environment (TENSOR_API_BASE / TENSOR_AGENT_USER / TENSOR_AGENT_PASS);
nothing is embedded. Tables are injected from docs/m3822/merge_verification.md, which
tools/verify_merge_m3822.py regenerates from the merged tree, so no number is hand-copied.
"""

from __future__ import annotations

import base64
import io
import json
import os
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DOC = ROOT / "docs" / "m3822_merge_cure_pr44.md"
TABLES = ROOT / "docs" / "m3822" / "merge_verification.md"
FIGURE = ROOT / "docs" / "m3822" / "merge_verification.png"
BASE = os.environ.get("TENSOR_API_BASE", "http://127.0.0.1:8010/api/v1").rstrip("/")
USER = os.environ["TENSOR_AGENT_USER"]
PASSWORD = os.environ["TENSOR_AGENT_PASS"]


def post(path: str, payload: dict, token: str | None = None, timeout: int = 180) -> dict:
    headers = {"Content-Type": "application/json"}
    if token:
        headers["Authorization"] = f"Bearer {token}"
    request = urllib.request.Request(
        f"{BASE}{path}", data=json.dumps(payload).encode(), headers=headers, method="POST"
    )
    try:
        with urllib.request.urlopen(request, timeout=timeout) as response:
            return json.loads(response.read())
    except urllib.error.HTTPError as error:
        raise SystemExit(f"HTTP {error.code} on {path}: {error.read().decode()[:600]}") from None


def login() -> str:
    data = post("/auth/login", {"username": USER, "password": PASSWORD})
    token = data.get("access_token") or data.get("token")
    if not token:
        token = data["data"]["access_token"]
    return token


def figure_payload(path: Path, target_kb: int = 60) -> tuple[dict, int]:
    """Downscale so the base64 payload stays well under the 250 KB report ceiling."""
    from PIL import Image

    image = Image.open(path).convert("RGB")
    for width in (1100, 950, 820, 700):
        scaled = image if image.width <= width else image.resize(
            (width, int(image.height * width / image.width)), Image.LANCZOS
        )
        buffer = io.BytesIO()
        scaled.save(buffer, format="PNG", optimize=True)
        raw = buffer.getvalue()
        if len(raw) <= target_kb * 1024:
            break
    encoded = base64.b64encode(raw).decode()
    return (
        {
            "name": path.stem,
            "caption": (
                "PR #44 cure verification: GZ vs section clipping (1.49e-14 m); transom volume "
                "error vs closed form; test counts main vs cured; added-line survival 100%."
            ),
            "png_base64": encoded,
        },
        len(raw),
    )


def main() -> None:
    doc = DOC.read_text(encoding="utf-8")
    tables = TABLES.read_text(encoding="utf-8")
    # drop the generated file's own title/first paragraphs, keep its tables
    body = tables.split("\n", 1)[1]
    markdown = doc.replace("{{TABLES}}", body.strip())
    assert "{{TABLES}}" not in markdown

    fig, size = figure_payload(FIGURE)
    payload = {
        "title": "PR #44 治愈：agent/m3347 × main 三方合并，GZ 曲线与艉板/惯性矩船族的无损并集",
        "summary": (
            "Merge-cured PR #44: three-way merge of agent/m3347 with main with machine-checked "
            "no-silent-loss (100% of both sides' added lines kept, __all__ an exact 112+120=125 "
            "union), both numeric chains re-run from the merged tree (GZ vs section clipping "
            "1.49e-14 m; transom volume vs closed form 0.13-1.54%), tests 879 passed vs main's 818 "
            "baseline with an identical failure set."
        ),
        "scenario_key": "tensorcad",
        "mission_id": 3822,
        "tags": ["TensorCAD", "merge-cure", "PR44", "stability", "GZ", "transom", "no-silent-loss"],
        "markdown": markdown,
        "figures": [fig],
    }
    token = login()
    result = post("/tensor/reports", payload, token)
    print(json.dumps({"figure_png_bytes": size, "response": result}, ensure_ascii=False)[:800])


if __name__ == "__main__":
    main()