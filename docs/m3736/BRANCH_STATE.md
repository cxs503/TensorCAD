# Branch/SHA record — PR #37 cure (m3736)

| ref | SHA | meaning |
|---|---|---|
| `origin/agent/m2873` (PR #37 head, **old**) | `206698a` | "m2873: rebuild on current main 6c568c31"; base `6c568c3`, so it collides with everything main added since |
| `origin/main` at cure time | `8389c7e` | merge target |
| `agent/m3736` (**new**, this work) | see `git rev-parse HEAD` | single commit whose parent *is* `8389c7e` and whose tree is the semantic union: main's 27 newer paths + the branch's 15 paths + 4 hand-resolved collisions |

Why the cure is a single commit on top of `origin/main` and not a merge commit on
`agent/m2873`: the platform pushes the mission branch (`agent/m<mid>`) as a commit whose
parent is the current `main` tip, and a branch built that way is mergeable by construction
(`mergeable=true`, empty conflict set). Re-parenting the *cured content* onto `main` loses
nothing — `docs/m3736/verify_union.py` proves the tree still contains every path either
side touched, byte-identical, and the audit is re-runnable (`silent_loss: []`).

Push is performed by the platform pusher; no branch was force-updated by hand here, and
`agent/m2873`'s ref itself was left untouched.
