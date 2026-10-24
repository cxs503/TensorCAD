# Commands used for this run

Environment: 2 vCPU / 3 GB RAM / no GPU, CPython 3.11.2, `torch 2.2.2+cu121`
(the CUDA runtime is a wheel dependency only), `safetensors 0.8.0`, `numpy 2.4.6`.

```bash
python -m venv .venv
.venv/bin/pip install -e ".[test,dev]"
# PyPI's current torch wheel pulls the whole CUDA toolkit; pin an older CPU-runnable
# wheel instead. torch 2.2.2 needs its NVIDIA runtime deps present to import.
.venv/bin/pip install "torch==2.2.2" safetensors numpy
```

## Smoke test (20 + 20 steps; expected to roll back)

```bash
python examples/neural_geometry/edit_sphere_dent.py \
  --base-steps 20 --edit-steps 20 --seed 503 \
  --output /tmp/tc-exp/smoke.ngo
# exit 1, transaction=rolled_back, no bundle directory
```

## Baseline sphere (revision 1)

```bash
python examples/neural_geometry/train_sphere.py \
  --steps 1800 --seed 503 --output results/neural-sphere.ngo
```

## Formal dent experiment (revision 2)

```bash
python examples/neural_geometry/edit_sphere_dent.py \
  --base-steps 1200 --edit-steps 1000 --seed 503 \
  --output results/sphere-dent-seed-503.ngo
```

## Determinism repeat

```bash
python examples/neural_geometry/edit_sphere_dent.py \
  --base-steps 1200 --edit-steps 1000 --seed 503 \
  --output /tmp/tc-exp/dent-repeat.ngo
sha256sum results/sphere-dent-seed-503.ngo/weights.safetensors \
          /tmp/tc-exp/dent-repeat.ngo/weights.safetensors
```

## Deliberate rollback cases

```bash
# Candidate surface leaves the 0.36-0.58 m radial search window: zero dent rays.
python examples/neural_geometry/edit_sphere_dent.py \
  --base-steps 1200 --edit-steps 1000 --seed 503 \
  --depth 0.25 --sigma 0.15 --output /tmp/tc-exp/rollback1.ngo

python examples/neural_geometry/edit_sphere_dent.py \
  --base-steps 1200 --edit-steps 1000 --seed 503 \
  --depth 0.45 --sigma 0.12 --output /tmp/tc-exp/rollback2.ngo

# Rollback targeting a path that already holds a committed revision 2.
python examples/neural_geometry/edit_sphere_dent.py \
  --base-steps 1200 --edit-steps 1000 --seed 503 \
  --depth 0.45 --sigma 0.12 --output results/sphere-dent-seed-503.ngo
sha256sum results/sphere-dent-seed-503.ngo/weights.safetensors   # unchanged

# The JSON must stay parseable even when metrics are inf/NaN.
python - <<'PY'
import json
for name in ("rollback1", "rollback2", "rollback-existing"):
    text = open(f"/tmp/tc-exp/{name}.stdout.json").read()
    json.loads(text, parse_constant=lambda c: (_ for _ in ()).throw(ValueError(c)))
    print(name, "valid JSON, no NaN/Infinity tokens")
PY
```

## Independent verification

```bash
cd docs/experiments/2026-10-24-neural-sphere-dent
python verification/verify_profile.py --parent-out /tmp/parent.safetensors
python verification/protected_stats.py --parent /tmp/parent.safetensors
```

`verify_profile.py` replays the base training from seed 503 / 1200 steps, probes the
surface radius at fixed polar angles with a dense scan + bisection, and writes
`verification/angular_profile.json`. `protected_stats.py` samples 4096 fixed-seed
directions and reports mean/p95/p99/max displacement plus the baseline's own
radial form error, writing `verification/protected_region_stats.json`.

## Integrity checks

```bash
python - <<'PY'
import hashlib, json
from pathlib import Path
from tensorcad.neural_geometry import load_bundle

root = Path("docs/experiments/2026-10-24-neural-sphere-dent/bundles")
for name in ("neural-sphere.ngo", "sphere-dent-seed-503.ngo"):
    bundle = load_bundle(root / name)          # verifies SHA-256 and metadata paths
    actual = hashlib.sha256((root / name / "weights.safetensors").read_bytes()).hexdigest()
    print(name, "revision", bundle.revision, "sha256", actual)
PY
```
