# Benchmarks / baseline_fleet — published geometry vs reconstructed displacement

Case: `GAP-CAD-T1-BASELINE-FLEET`. Documentation: `docs/baseline_fleet_build.md`.

## What is compared

| hull | ref (published ∇, m³) | got (mesh ∇, m³) | rel. error | watertight | verdict |
|---|---|---|---|---|---|
| KCS | 52 030 | 52 187.98 | **+0.304 %** | yes | VERIFIED |
| JBC | 178 369.9 | — | — | — | WITHHELD (IGES control-net layout guard fired) |
| ONRT | 8 500.8 (from published Cb·Lwl·Bwl·T) | — | — | — | WITHHELD (model-scale coordinates, no scale factor obtained) |
| KVLCC2 | no published table reachable | — | — | — | WITHHELD (no reference) |
| DTC | — | — | — | — | NOT SUPPORTED (no official geometry located) |
| Series 60 | — | — | — | — | NOT SUPPORTED (no official IGES; only the analytic generator in `tensorcad/series60.py`) |

Reference definition: the publisher's own displacement from the workshop page cached
under `artifacts/baseline_fleet/raw/pages/`; every value carries the page SHA-256 and the
matched snippet in `artifacts/baseline_fleet/fleet_specs.json` (nothing is hand-typed).

Secondary published checks on KCS (same page): Cb 0.6505 (got 0.6516, +0.17 %),
LCB −1.48 % Lpp (got −1.48 %), LWL 232.5 m (cut at the design waterline: 230.25 m).

## Reproduce

```bash
python tools/baseline_fleet/fetch_sources.py --out artifacts/baseline_fleet   # network, records SHA-256
python tools/baseline_fleet/build_fleet.py --samples-per-metre 6              # rebuild + benchmark
python -m pytest tests/test_iges_baseline.py -q                               # 4 passed, 1 xfail
```

`result.json` is the machine-readable copy of `artifacts/baseline_fleet/benchmark.json`
produced by the run above (2026-11-13, `--samples-per-metre 6`, dx = 0.5 m, dz = 0.1 m,
240 ring points).
