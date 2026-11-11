#!/usr/bin/env python3
"""Live survey of open CAD / geometry-kernel repositories for FEAT-OSS-INVENTORY-CAD.

Why this exists
---------------
The census must be reproducible and must not depend on hand-typed star counts or
memorised licence names.  Every field in ``artifacts/oss-inventory-cad-v1/raw/*.json``
is obtained from the repository itself at a pinned commit:

* ``commit`` / ``committed_at``  - a shallow blob-less clone (``--filter=blob:none``)
  when git handshakes fast enough, otherwise the GitHub REST API (metadata + one
  recursive tree call), with the same commit SHA recorded either way;
* ``paths``                      - the real file tree (``git ls-tree -r`` or
  ``git/trees?recursive=1``);
* ``license_file`` / ``license_text_head`` / ``license_tokens`` - the licence blob
  actually stored in the tree, classified by *text* (never by memory);
* per feature point ``candidates`` - paths whose name matches that feature point's
  signature (``tools/oss_inventory_cad/feature_points.json``), so "repo X implements
  feature F" is derived from X's own file names;
* ``pointer``                    - ``path:line`` where the feature probe token actually
  occurs in the blob at that commit (lazy blob fetch via git, or a raw.githubusercontent
  fetch at the pinned SHA), i.e. a citable pointer rather than a guess.

A repository that cannot be reached at all is recorded with ``status: clone_failed`` /
``api_failed`` plus the error text - the negative control in repos.json exists precisely
to prove that failures are recorded rather than silently dropped from the census.

Usage
-----
    python tools/oss_inventory_cad/survey.py                  # incremental
    python tools/oss_inventory_cad/survey.py --budget 80      # stop after ~80s (resumable)
    python tools/oss_inventory_cad/survey.py --refresh        # re-clone everything
    python tools/oss_inventory_cad/survey.py --only Open-Cascade-SAS/OCCT
"""
from __future__ import annotations

import argparse
import json
import os
import re
import shutil
import subprocess
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
TOOL = Path(__file__).resolve().parent
OUT = ROOT / "artifacts" / "oss-inventory-cad-v1" / "raw"
CACHE = Path(os.environ.get("OSS_CAD_CACHE", "/tmp/oss-cad-cache"))

SOURCE_EXT = {
    ".py", ".h", ".hpp", ".hh", ".hxx", ".cpp", ".cxx", ".cc", ".c",
    ".inl", ".ipp", ".txx", ".ts", ".tsx", ".js", ".mjs", ".rs", ".java", ".cs",
    ".f90", ".f", ".f95", ".jl", ".go", ".mm", ".m", ".cu", ".cl", ".glsl",
}
#: Directories whose contents are *users* of a kernel, not the kernel: keeping them
#: would let "tests/test_boolean.py" masquerade as a boolean implementation.
EXCLUDE_SEG = {
    "test", "tests", "testing", "doc", "docs", "documentation", "example",
    "examples", "sample", "samples", "demo", "demos", "benchmark", "benchmarks",
    "third_party", "thirdparty", "3rdparty", "extern", "external", "submodules",
    "node_modules", "vendor", ".git", "ci", ".github",
}

LICENSE_NAMES = re.compile(
    r"^(licen[cs]e|copying|copyright|licence|unlicense|notice)(\.(md|txt|rst|html))?"
    r"(\.[a-z0-9]+)?$", re.I)

#: Ordered most-restrictive-first; the licence *class* is the strictest token found,
#: which is the conservative reading for a dual-licensed tree.
LICENSE_TOKENS = [
    ("AGPL-3.0", re.compile(r"GNU AFFERO GENERAL PUBLIC LICENSE", re.I)),
    ("GPL-3.0", re.compile(r"GNU GENERAL PUBLIC LICENSE\s*(?:Version|v\.?)\s*3", re.I)),
    ("GPL-2.0", re.compile(r"GNU GENERAL PUBLIC LICENSE\s*(?:Version|v\.?)\s*2", re.I)),
    ("LGPL-3.0", re.compile(r"GNU LESSER GENERAL PUBLIC LICENSE\s*(?:Version|v\.?)\s*3", re.I)),
    ("LGPL-2.1", re.compile(r"GNU LESSER GENERAL PUBLIC LICENSE\s*(?:Version|v\.?)\s*2", re.I)),
    ("MPL-2.0", re.compile(r"Mozilla Public License Version 2\.0", re.I)),
    ("EPL-2.0", re.compile(r"Eclipse Public License", re.I)),
    ("Apache-2.0", re.compile(r"Apache License,?\s*Version 2\.0", re.I)),
    ("BSD-3-Clause", re.compile(r"Redistribution and use in source and binary forms.{0,600}?"
                               r"Neither the name", re.I | re.S)),
    ("BSD-2-Clause", re.compile(r"Redistribution and use in source and binary forms", re.I)),
    ("MIT", re.compile(r"Permission is hereby granted, free of charge", re.I)),
    ("NOSA-1.3", re.compile(r"NASA Open Source Agreement", re.I)),
    ("CC0-1.0", re.compile(r"CC0 1\.0 Universal", re.I)),
    ("Unlicense", re.compile(r"This is free and unencumbered software", re.I)),
]

#: Licence *policy* for reuse inside TensorCAD (see docs/oss_inventory_cad/ REPORT §4).
CLASS_POLICY = {
    "permissive": "fits",
    "copyleft-weak": "adapt",
    "copyleft-strong": "reference",
    "no-license": "reference",
    "restricted": "reference",
}


SPDX_RE = re.compile(
    r"\b(AGPL-3\.0|GPL-[23]\.0(?:-only|-or-later)?|LGPL-[23]\.[01]|MPL-2\.0|"
    r"Apache-2\.0|BSL-1\.0|MIT|BSD-[23]-Clause)\b")


def classify_license(text: str) -> tuple[list[str], str]:
    """Return (tokens found, class). Class is the strictest token present."""
    tokens = [name for name, rx in LICENSE_TOKENS if rx.search(text)]
    if not tokens:
        # Some trees ship a *summary* file (CGAL: LICENSE.md lists per-package SPDX
        # ids).  Fall back to SPDX identifiers before declaring "no licence".
        spdx = sorted(set(m.group(1) for m in SPDX_RE.finditer(text)))
        if spdx:
            tokens = spdx
        else:
            # Last resort: prose mention (CGAL's LICENSE.md describes the licence in
            # running text - "GNU Lesser General Public License ... v3" - without an
            # SPDX id in the first kilobytes read).
            low = text.lower()
            if "lgpl" in low or "lesser general public license" in low:
                tokens = ["LGPL(prose)"]
            elif "gpl" in low or "general public license" in low:
                tokens = ["GPL(prose)"]
    if not tokens:
        return [], "no-license"
    if any(t.startswith("AGPL") or t.startswith("GPL") for t in tokens):
        return tokens, "copyleft-strong"
    if any(t.startswith(("LGPL", "MPL", "EPL")) for t in tokens):
        return tokens, "copyleft-weak"
    if any(t.startswith("NOSA") for t in tokens):
        return tokens, "restricted"
    return tokens, "permissive"


# --------------------------------------------------------------------------- backends
class GitRepo:
    """Shallow blob-less clone; blobs are lazily fetched per pointer request."""

    kind = "git"

    def __init__(self, repo: str, dest: Path, timeout: int = 55):
        self.dest, self.timeout = dest, timeout
        shutil.rmtree(dest, ignore_errors=True)
        p = subprocess.run(
            ["git", "clone", "--depth", "1", "--filter=blob:none", "--no-checkout", "-q",
             f"https://github.com/{repo}.git", str(dest)],
            capture_output=True, text=True, timeout=timeout)
        if p.returncode != 0:
            err = (p.stderr or p.stdout or "").strip().splitlines()
            raise RuntimeError("git clone failed: " + (" | ".join(err[-2:]) if err else "exit"))
        g = self.git(["log", "-1", "--format=%H%n%cI"])
        sha, _, date = g.stdout.strip().partition("\n")
        if not sha:
            raise RuntimeError("git clone produced no HEAD")
        self.commit, self.committed_at = sha, date
        self._paths = None

    def git(self, args, timeout=180):
        return subprocess.run(["git", *args], cwd=str(self.dest), capture_output=True,
                              text=True, timeout=timeout)

    @property
    def paths(self) -> list[str]:
        if self._paths is None:
            p = self.git(["ls-tree", "-r", "--name-only", "HEAD"], timeout=180)
            if p.returncode != 0:
                raise RuntimeError(f"ls-tree failed: {p.stderr.strip()[:200]}")
            self._paths = [ln for ln in p.stdout.splitlines() if ln]
        return self._paths

    def show(self, path: str) -> str | None:
        p = self.git(["show", f"HEAD:{path}"], timeout=180)
        return p.stdout[:400_000] if p.returncode == 0 else None


class RawRepo:
    """API metadata + one recursive tree call; blobs from raw.githubusercontent.com.

    Used only when the git handshake does not finish inside the per-repo budget -
    the platform's egress path stalls on a few very large repositories.  The commit
    SHA is pinned identically, so pointers stay citable.
    """

    kind = "api+raw"

    def __init__(self, repo: str, timeout: int = 30):
        meta = self._json(f"https://api.github.com/repos/{repo}", timeout)
        self.branch = meta.get("default_branch") or "main"
        tree = self._json(
            f"https://api.github.com/repos/{repo}/git/trees/{self.branch}?recursive=1",
            timeout)
        self._paths = [n["path"] for n in tree.get("tree", []) if n.get("type") == "blob"]
        if not self._paths:
            raise RuntimeError("api tree empty")
        self.commit = tree.get("sha") or meta.get("default_branch")
        self.committed_at = meta.get("pushed_at") or ""
        self.api_meta = {
            "stars": meta.get("stargazers_count"),
            "size_kb": meta.get("size"),
            "archived": meta.get("archived"),
            "spdx": (meta.get("license") or {}).get("spdx_id"),
            "html_language": meta.get("language"),
        }
        self.repo = repo

    @staticmethod
    def _json(url: str, timeout: int) -> dict:
        req = urllib.request.Request(url, headers={"User-Agent": "tensorcad-oss-inventory/1"})
        try:
            with urllib.request.urlopen(req, timeout=timeout) as r:
                return json.loads(r.read().decode())
        except urllib.error.HTTPError as e:
            raise RuntimeError(f"HTTP {e.code} for {url}") from None
        except urllib.error.URLError as e:
            raise RuntimeError(f"URLError {e.reason} for {url}") from None

    @property
    def paths(self) -> list[str]:
        return self._paths

    def show(self, path: str) -> str | None:
        q = urllib.parse.quote(path, safe="/")
        url = f"https://raw.githubusercontent.com/{self.repo}/{self.commit}/{q}"
        req = urllib.request.Request(url, headers={"User-Agent": "tensorcad-oss-inventory/1"})
        try:
            with urllib.request.urlopen(req, timeout=30) as r:
                return r.read(400_000).decode("utf-8", "replace")
        except Exception:  # noqa: BLE001 - a missing blob is data, not a crash
            return None


# --------------------------------------------------------------------------- helpers
def is_source(path: str) -> bool:
    if Path(path).suffix.lower() not in SOURCE_EXT:
        return False
    segs = {s.lower() for s in Path(path).parts[:-1]}
    return not (segs & EXCLUDE_SEG)


def score(path: str, sigs: list[str]) -> int:
    low, name = path.lower(), Path(path).name.lower()
    s = 0
    for sig in sigs:
        if sig in name:
            s += 6
        elif sig in low:
            s += 2
    return s


def pick_candidates(paths: list[str], sigs: list[str], limit: int = 3) -> list[str]:
    scored = [(-score(p, sigs), len(p), p) for p in paths if score(p, sigs) > 0]
    scored.sort()
    return [p for _, _, p in scored[:limit]]


def find_pointer(source, path: str, probes: list[str]) -> dict:
    """Locate the first line in ``path`` carrying a feature probe token."""
    blob = source.show(path)
    if blob is None:
        return {"pointer": None, "pointer_status": "blob_fetch_failed"}
    lines = blob.splitlines()
    for tok in probes:
        rx = re.compile(re.escape(tok))
        for i, line in enumerate(lines, 1):
            if rx.search(line):
                hi = min(len(lines), i + 7)
                return {"pointer": f"{path}:{i}-{hi}", "pointer_status": "blob_verified",
                        "pointer_token": tok, "pointer_line_text": line.strip()[:180]}
    return {"pointer": None, "pointer_status": "keyword_absent_in_blob",
            "pointer_n_lines": len(lines)}


def license_of(source, paths: list[str]) -> dict:
    cands = [p for p in paths if LICENSE_NAMES.match(Path(p).name)]
    # Rank by *name kind* first: a bare "copyright" notice is a file-header summary,
    # not the licence grant (CGAL's root "copyright" lists authors only, so picking it
    # silently turned a dual GPL/LGPL tree into "no-license").
    def rank(p: str) -> tuple:
        base = Path(p).name.lower().split(".")[0]
        kind = 0 if base.startswith(("license", "licence")) else (
            1 if base.startswith("copying") else 2)
        return (p.count("/"), kind, len(p))
    cands.sort(key=rank)
    if not cands:
        return {"license_file": None, "license_tokens": [], "license_class": "no-license",
                "license_text_head": None, "license_status": "no_license_file_in_tree"}
    text = None
    for cand in cands[:4]:
        text = source.show(cand)
        if text and len(text.strip()) > 40:
            break
    if not text:
        return {"license_file": cands[0], "license_tokens": [], "license_class": "unknown",
                "license_text_head": None, "license_status": "blob_unavailable"}
    tokens, klass = classify_license(text)
    if not tokens:
        # A licence file exists but carries no recognisable grant in the bytes read:
        # "unknown" is honest, "no-license" would wrongly imply all-rights-reserved.
        klass = "unknown"
    return {"license_file": cands[0], "license_tokens": tokens, "license_class": klass,
            "license_text_head": " ".join(text.split())[:220],
            "license_status": "classified", "license_candidates": cands[:4]}


def survey_repo(entry: dict, features: list[dict], refresh: bool, budget: int = 55) -> dict:
    repo, dest = entry["repo"], CACHE / re.sub(r"[^A-Za-z0-9._-]+", "__", entry["repo"])
    t0 = time.time()
    rec = {"repo": repo, "tier": entry["tier"], "focus": entry.get("focus", []),
           "why_included": entry.get("why", ""), "url": f"https://github.com/{repo}",
           "surveyed_at": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())}

    source, git_err = None, None
    if dest.is_dir() and not refresh:
        try:
            source = GitRepo.__new__(GitRepo)          # reuse cache without re-cloning
            source.dest, source.timeout, source._paths = dest, budget, None
            g = source.git(["log", "-1", "--format=%H%n%cI"])
            sha, _, date = g.stdout.strip().partition("\n")
            if sha:
                source.commit, source.committed_at = sha, date
            else:
                source = None
        except Exception:  # noqa: BLE001
            source = None
    if source is None:
        try:
            source = GitRepo(repo, dest, timeout=budget)
        except Exception as exc:  # noqa: BLE001 - fall back to API+raw
            git_err = f"{type(exc).__name__}: {str(exc)[:200]}"
            try:
                source = RawRepo(repo, timeout=25)
            except Exception as exc2:  # noqa: BLE001
                rec.update({"status": "clone_failed", "clone_error": git_err,
                            "api_error": f"{type(exc2).__name__}: {str(exc2)[:200]}",
                            "elapsed_s": round(time.time() - t0, 2)})
                return rec

    paths = source.paths
    src = [p for p in paths if is_source(p)]
    exts = {Path(p).suffix.lower() for p in src}
    rec.update({
        "status": "surveyed", "backend": source.kind,
        "clone_error": git_err,
        "commit": source.commit, "committed_at": source.committed_at,
        "n_paths": len(paths), "n_source_paths": len(src),
        "dominant_ext": (max(exts, key=lambda e: sum(1 for p in src if p.endswith(e)))
                         if src else None),
    })
    if getattr(source, "api_meta", None):
        rec["api_meta"] = source.api_meta
    rec.update(license_of(source, paths))

    per_feature = {}
    for f in features:
        cands = pick_candidates(src, f["sig"])
        row = {"candidates": cands, "n_candidates": len(cands)}
        row.update(find_pointer(source, cands[0], f["probe"]) if cands
                   else {"pointer": None, "pointer_status": "no_path_signature"})
        per_feature[f["id"]] = row
    rec["features"] = per_feature
    rec["elapsed_s"] = round(time.time() - t0, 2)
    return rec


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--refresh", action="store_true")
    ap.add_argument("--only", default=None)
    ap.add_argument("--limit", type=int, default=0)
    ap.add_argument("--budget", type=int, default=80,
                    help="soft wall-clock budget for this invocation (resumable)")
    ap.add_argument("--per-repo", type=int, default=55, help="per-repo git clone budget (s)")
    args = ap.parse_args(argv)

    repos = json.loads((TOOL / "repos.json").read_text())["repos"]
    features = json.loads((TOOL / "feature_points.json").read_text())["feature_points"]
    if args.only:
        repos = [r for r in repos if r["repo"] == args.only]
    if args.limit:
        repos = repos[: args.limit]
    OUT.mkdir(parents=True, exist_ok=True)

    t_start = time.time()
    for i, entry in enumerate(repos, 1):
        path = OUT / (re.sub(r"[^A-Za-z0-9._-]+", "__", entry["repo"]) + ".json")
        if path.exists() and not args.refresh:
            prev = json.loads(path.read_text())
            if prev.get("status") == "surveyed":
                print(f"[{i}/{len(repos)}] {entry['repo']}: cached "
                      f"({str(prev.get('commit'))[:10]})", flush=True)
                continue
        if time.time() - t_start > args.budget:
            print(f"[budget] stopping before {entry['repo']} "
                  f"({round(time.time() - t_start)}s > {args.budget}s)", flush=True)
            break
        print(f"[{i}/{len(repos)}] {entry['repo']}: surveying ...", flush=True)
        try:
            rec = survey_repo(entry, features, args.refresh, budget=args.per_repo)
        except Exception as exc:  # noqa: BLE001 - a survey failure is data, not a crash
            rec = {"repo": entry["repo"], "tier": entry["tier"],
                   "focus": entry.get("focus", []), "why_included": entry.get("why", ""),
                   "url": f"https://github.com/{entry['repo']}", "status": "survey_error",
                   "clone_error": f"{type(exc).__name__}: {str(exc)[:200]}"}
        path.write_text(json.dumps(rec, ensure_ascii=False, indent=1) + "\n")
        print(f"    -> {rec['status']} via {rec.get('backend', '-')} "
              f"({rec.get('elapsed_s', '?')}s) {str(rec.get('commit', ''))[:10]}", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())