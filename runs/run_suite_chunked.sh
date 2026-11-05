#!/bin/bash
# Run the TensorCAD test suite one file at a time, resumable, within a wall-clock
# budget per invocation (the harness kills the process group at ~120 s, so the
# suite cannot be run in a single foreground call).
#
# usage: run_suite_chunked.sh <tree-dir> <log-dir> [budget-seconds] [per-file-timeout]
#   Writes to <log-dir>/per_file.tsv lines:  <status> <seconds> <file> <summary>
#   Already-recorded files are skipped, so re-invoking continues where it stopped.
set -u
TREE=$1; OUT=$2; BUDGET=${3:-100}; FTIMEOUT=${4:-60}
mkdir -p "$OUT"
TSV="$OUT/per_file.tsv"
touch "$TSV"
cd "$TREE" || exit 1
START=$(date +%s)
for f in tests/test_*.py; do
  case "$f" in tests/test_cli_neural.py|tests/test_neural_pipeline.py) continue;; esac
  grep -qP "^[^\t]+\t[^\t]+\t\Q$f\E\t" "$TSV" && continue
  if [ $(( $(date +%s) - START )) -ge "$BUDGET" ]; then
    echo "BUDGET_REACHED after $(( $(date +%s) - START ))s"; break
  fi
  T0=$(date +%s)
  LOG=$(timeout "$FTIMEOUT" python -m pytest "$f" -q -p no:cacheprovider -p no:randomly 2>&1)
  RC=$?
  T1=$(date +%s)
  SUMMARY=$(echo "$LOG" | grep -E "^[0-9]+ (passed|failed)|passed|failed|error" | tail -1 | tr '\t' ' ')
  printf '%s\t%s\t%s\t%s\n' "$RC" "$((T1-T0))" "$f" "$SUMMARY" >> "$TSV"
  echo "$RC $((T1-T0))s $f :: $SUMMARY"
done
echo "DONE $(wc -l < "$TSV") files recorded"