# PR #40 cure — publication refs

The cure is a two-parent merge commit that puts `main` back inside `agent/m3169`:

| ref | SHA |
|---|---|
| `agent/m3169` before the cure (PR #40 head) | `fcf15adccd0dcbae671868dcf51e064c5e398946` |
| cure merge commit (parents `fcf15ad` and `2ff80b2` = `main`) | `1b8fd067f16744d77b6e52b3f21cbc3886f06238` |
| cure branch tip with the evidence commit (`cure/m3169`) | `a05702a` |
| evidence mirror on `agent/m3777` (this branch) | `bc2e482` |

`fcf15ad` is an ancestor of every commit above, so updating `agent/m3169` to `a05702a`
is a plain fast-forward: no force, and no other branch moves. Publication is the
platform's push step (the agent shell is not permitted to push), so the merge is left
committed and verified in the worktree.

Reproduce the evidence (must run on the cured branch, it identifies the merge by its
second parent being `origin/main`):

```sh
python docs/merge_cure/run_cure_evidence.py
```