"""Run the geometry-representation benchmark and write results.json + report.md.

What this runner is for (issue #19): put every representation TensorCAD can produce
through one harness, publish the numbers, and let the numbers decide which representation
suits which job. It is deliberately *not* a leaderboard: a row that cannot be measured here
is reported as `not_supported` with the command that would measure it, and a gate that
fails stays failed.

```bash
# what CI runs: both cases, small sampling, all gates
python examples/benchmarks/run_benchmark.py --quick --output results/bench

# the reference run
python examples/benchmarks/run_benchmark.py --output results/bench

# include a trained neural object and the optional kernel path
python examples/neural_geometry/edit_sphere_dent.py --output results/sphere-dent.ngo   # needs .[neural]
python examples/benchmarks/run_benchmark.py --neural-bundle results/sphere-dent.ngo \
    --output results/bench                                              # STEP row needs .[ocp]
```

Exit code is 1 when any gate fails, so the command is usable as a CI step. The artifacts
written are `results.json` (machine readable), `report.md` and `plot.svg` (both generated
from the same data) and, per case, small OBJ/msh files when `--artifacts` is passed.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from tensorcad.benchmark import run_benchmark


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--quick", action="store_true", help="small sampling for CI; gates are unchanged")
    parser.add_argument("--output", default="results/bench", help="directory for results.json and report.md")
    parser.add_argument("--neural-bundle", default=None, help="a .ngo bundle to import neural metrics from")
    parser.add_argument("--ocp-artifact", default=None, help="where to write the STEP file, when OCP is installed")
    parser.add_argument("--artifacts", action="store_true", help="also write the meshes and gmsh/STEP files")
    args = parser.parse_args(argv)

    output = Path(args.output)
    result = run_benchmark(
        quick=args.quick,
        neural_bundle=args.neural_bundle,
        ocp_artifact=args.ocp_artifact,
        write_artifacts=args.artifacts,
        artifact_dir=output if args.artifacts else None,
        command="python examples/benchmarks/run_benchmark.py " + " ".join(argv if argv is not None else sys.argv[1:]),
    )
    json_path = result.write_json(output / "results.json")
    report_path = result.write_markdown(output / "report.md")
    plot_path = result.write_svg(output / "plot.svg")

    environment = result.environment
    print(
        f"tensorcad {environment['tensorcad']} | python {environment['python']} | "
        f"numpy {environment['numpy']} | torch {environment['torch']} | ocp {environment['ocp']}"
    )
    for case in result.cases:
        print(f"\n{case.identifier}  ({case.seconds:.2f} s)")
        print(f"  {case.question}")
        for representation in case.representations:
            passed = sum(1 for metric in representation.metrics.values() if metric.passed is True)
            gated = sum(1 for metric in representation.metrics.values() if metric.passed is not None)
            print(
                f"    {representation.identifier:18s} {representation.status:14s} "
                f"{representation.seconds:6.2f} s   gates {passed}/{gated}"
            )
            for name, metric in representation.metrics.items():
                if metric.value is None:
                    continue
                marker = "" if metric.passed is None else ("  ok" if metric.passed else "  FAIL")
                print(f"        {name:36s} {metric.value:14.6g} {metric.unit:8s} {metric.status}{marker}")
        for gate in case.gates:
            print(f"    [{'PASS' if gate.passed else 'FAIL'}] {gate.name}: {gate.detail}")

    print(f"\nwrote {json_path}")
    print(f"wrote {report_path}")
    print(f"wrote {plot_path}")
    if result.passed:
        print(f"all gates passed in {result.seconds:.2f} s")
        return 0
    print(f"FAILED ({len(result.failures)}): " + "; ".join(result.failures))
    return 1


if __name__ == "__main__":
    raise SystemExit(main())
