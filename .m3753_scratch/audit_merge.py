#!/usr/bin/env python3
"""Audit a merge commit for silent content loss against both parents."""
import subprocess, sys, re

REPO = "/root/tp-agent/work/m3753"
BASE, MAIN, BRANCH, MERGE = "ca0dde3", "311f45a", "ad86f4f", "546e001"


def blob(rev, path):
    r = subprocess.run(["git", "-C", REPO, "show", f"{rev}:{path}"],
                       capture_output=True, text=True)
    return r.stdout if r.returncode == 0 else None


def diff_added(rev, path):
    r = subprocess.run(["git", "-C", REPO, "diff", "-U0", BASE, rev, "--", path],
                       capture_output=True, text=True)
    return [ln[1:] for ln in r.stdout.splitlines()
            if ln.startswith("+") and not ln.startswith("+++")]


def norm(s):
    return re.sub(r"\s+", " ", s).strip()


paths = [p.strip() for p in open("/tmp/conflicts.txt") if p.strip()]
bad = 0
for p in paths:
    merged = blob(MERGE, p)
    if merged is None:
        print(f"[MISSING-IN-MERGE] {p}")
        bad += 1
        continue
    mset = {norm(l) for l in merged.splitlines()}
    mtext = "".join(merged.split())
    for side, rev in (("main", MAIN), ("branch", BRANCH)):
        lost = [l for l in diff_added(rev, p) if norm(l) and norm(l) not in mset]
        lost2 = [l for l in lost if "".join(l.split()) not in mtext]
        if lost2:
            bad += 1
            print(f"[LOST-FROM-{side}] {p}: {len(lost2)} line(s)")
            for l in lost2[:25]:
                print("    -", l[:150])
print("files:", len(paths), "problems:", bad)
sys.exit(0 if bad == 0 else 1)