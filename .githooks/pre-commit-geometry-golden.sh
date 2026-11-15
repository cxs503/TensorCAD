#!/usr/bin/env bash
# pre-commit geometry golden gate (GAP-CAD-D2-GOLDEN-REGRESSION), dual-derived arm.
#
# PROVENANCE.  Drop-in authored in TensorPapers at
# papers/gap-cad-d2-golden-regression/ci/pre-commit-geometry-golden.sh, installed here as
# a husky-style chain step (see the tail of .githooks/pre-commit) instead of a private
# .git/hooks/pre-commit, because the repository already routes hooks through
# core.hooksPath=.githooks (tools/golden/install_git_hooks.sh) and a checkout of a Git
# directory is not a reviewable change.
#
# Behaviour: if any staged path can move geometry/hydrostatics, run the independent
# golden suite and the mutation arm from TensorPapers and block the commit on a non-zero
# exit.  Everything else commits untouched, so the hook costs nothing on documentation
# changes.  ~35 s when it does run (the repository's own gate in .githooks/pre-commit is
# the cheap ~5 s first line of defence; this is the independent second opinion).
#
# Environment:
#   GOLDEN_DIR   override the suite location (missing explicit path = hard error)
#   GOLDEN_STAGED_FILE  read staged paths from a file (test harness; never used by git)
#   PYTHON       interpreter (default: python3)
#
# A checkout with no TensorPapers next to it warns and lets the commit through: locally
# unenforceable is not the same as locally failed, and the CI workflow
# .github/workflows/geometry-golden.yml does fail closed.
set -euo pipefail

GOLDEN_STAGED_FILE="${GOLDEN_STAGED_FILE:-}"
# repo_root = the checkout this hook ships in (source of the gate's judging code);
# TENSORCAD_PATH = the library under test, which may be pointed elsewhere.
repo_root=$(git rev-parse --show-toplevel)
TENSORCAD_PATH="${TENSORCAD_PATH:-$repo_root}"
export TENSORCAD_PATH

# Paths whose change can move a golden number.
WATCHED='^(tensorcad/(hull|wigley|series60|dtmb5415|hydrostatics|hydrostatics_analytic|mesh|topology|exchange|geometry_audit|voxel)\.py)'

if [ -n "$GOLDEN_STAGED_FILE" ]; then
  staged=$(cat "$GOLDEN_STAGED_FILE")
else
  staged=$(git diff --cached --name-only --diff-filter=ACMR)
fi
if ! printf '%s\n' "$staged" | grep -Eq "$WATCHED"; then
  exit 0
fi

echo "[golden] geometry-affecting change detected:"
printf '%s\n' "$staged" | grep -E "$WATCHED" | sed 's/^/  /'

if [ -n "${GOLDEN_DIR:-}" ]; then
  # An explicitly named suite that is not there is a misconfiguration, not a reason to
  # silently downgrade the gate.
  if [ ! -f "$GOLDEN_DIR/suite.py" ]; then
    echo "[golden] GOLDEN_DIR=$GOLDEN_DIR holds no suite.py -- refusing to guess" >&2
    exit 1
  fi
  suite_dir="$GOLDEN_DIR"
else
  suite_dir=""
  for candidate in \
      "$TENSORCAD_PATH/../TensorPapers/papers/gap-cad-d2-golden-regression" \
      "$HOME/TensorPapers/papers/gap-cad-d2-golden-regression"; do
    if [ -f "$candidate/suite.py" ]; then
      suite_dir="$candidate"
      break
    fi
  done
  if [ -z "$suite_dir" ]; then
    echo "[golden] suite not found next to this checkout (looked for ../TensorPapers and" >&2
    echo "[golden] \$HOME/TensorPapers); this commit is not locally gated.  Set GOLDEN_DIR," >&2
    echo "[golden] or rely on .github/workflows/geometry-golden.yml, which fails closed." >&2
    exit 0
  fi
fi

python_bin=${PYTHON:-python3}
if ! command -v "$python_bin" >/dev/null 2>&1; then
  echo "[golden] no $python_bin on PATH -- this commit is not locally gated" >&2
  exit 0
fi

export PYTHONPATH="$TENSORCAD_PATH${PYTHONPATH:+:$PYTHONPATH}"

assert_report="$repo_root/tools/golden/gap_d2_assert_report.py"
if [ ! -f "$assert_report" ]; then
  echo "[golden] $assert_report is missing -- cannot judge the run, refusing" >&2
  exit 1
fi

# Exit codes are not the contract: suite.py exits 0 when a case stops being measured
# (SKIP/UNAVAILABLE), and the mutation arm's contract is its rates.  Both are judged
# from the reports by the same script CI uses, so the two carriers cannot drift.
"$python_bin" "$suite_dir/suite.py"   || { echo "[golden] BLOCKED: golden suite regressed" >&2; exit 1; }
"$python_bin" "$assert_report" "$suite_dir" --suite-only || {
  echo "[golden] BLOCKED: a golden case was not measured" >&2; exit 1; }
"$python_bin" "$suite_dir/mutants.py" || { echo "[golden] BLOCKED: detection power lost" >&2; exit 1; }
"$python_bin" "$assert_report" "$suite_dir" --mutants-only || {
  echo "[golden] BLOCKED: detection power below contract" >&2; exit 1; }

echo "[golden] ok"
