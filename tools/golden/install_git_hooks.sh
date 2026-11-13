#!/bin/sh
# Install the repository's git hooks (currently the golden geometry gate).
#
# This is the only place that touches git configuration, and it only sets the hook
# directory for this checkout: `git config core.hooksPath .githooks`. It never changes
# remotes, user identity or anything else. Run it once per clone/worktree:
#
#     tools/golden/install_git_hooks.sh
#     git config --get core.hooksPath      # -> .githooks
#
# Opting out for one commit: TENSORCAD_GOLDEN_SKIP=1 git commit ...
# The CI workflow `.github/workflows/golden-geometry.yml` enforces the same gate, so a
# checkout without the hook installed is still covered.
set -eu

repo_root=$(git rev-parse --show-toplevel)
cd "$repo_root"

if [ ! -x ".githooks/pre-commit" ]; then
    chmod +x ".githooks/pre-commit"
fi

git config core.hooksPath .githooks
printf 'installed: core.hooksPath=%s\n' "$(git config --get core.hooksPath)"
printf 'hook: %s\n' "$repo_root/.githooks/pre-commit"
printf 'self-test (docs-only staging must skip, geometry staging must run):\n'
printf 'docs/baseline_fleet.md\n' > /tmp/golden-staged-docs.txt
PYTHON=${PYTHON:-python3} TENSORCAD_GOLDEN_STAGED_FILE=/tmp/golden-staged-docs.txt \
    "$repo_root/.githooks/pre-commit"