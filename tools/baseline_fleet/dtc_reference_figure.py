"""Figure for the DTC reference-binding record (``docs/dtc_reference_binding.md``).

Two panels, both fed from bytes rather than typed numbers:

* **left** -- the sectional-area curve ``A(x)`` of the DTC hull, integrated from the *publisher*
  IGES archive (``igs_ordner.zip``) re-fetched by :mod:`tools.baseline_fleet.dtc_reference_bind`;
  its own integral is compared with the published displacement that the reference record uses.
* **right** -- the three geometric checks of
  ``benchmarks/verified/fleet_refs/dtc_hull_measured.json`` (beam, block coefficient,
  displacement) as signed deviations from the published values, against the ±1 % acceptance band.

Usage::

    python tools/baseline_fleet/dtc_reference_figure.py --out docs/figs/dtc_reference_bind.png
"""

from __future__ import annotations

import argparse
import json
import sys
import tempfile
import zipfile
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[2]
if str(REPO) not in sys.path:  # pragma: no cover - import plumbing
    sys.path.insert(0, str(REPO))

from tensorcad.fleet_baseline import (  # noqa: E402
    field_from_points,
    field_integrals,
    partition_surfaces,
    sample_patches,
)
from tensorcad.iges import read_iges  # noqa: E402

DEFAULT_ARCHIVE = Path("/root/.cache/hull_sources/dtc_igs_ordner.zip")
DEFAULT_MEASURED = REPO / "benchmarks" / "verified" / "fleet_refs" / "dtc_hull_measured.json"
DEFAULT_OUT = REPO / "docs" / "figs" / "dtc_reference_bind.png"


def hull_points(archive: Path, spacing: float = 1.0) -> np.ndarray:
    """Sample the hull member of the publisher archive into an ``(n, 3)`` array."""
    with zipfile.ZipFile(archive) as bundle:
        for info in bundle.infolist():
            if not info.filename.replace("\\", "/").endswith("DTC_hull.igs"):
                continue
            clone = zipfile.ZipInfo(info.filename.replace("\\", "/"))
            for attribute in ("compress_type", "header_offset", "CRC", "compress_size",
                              "file_size", "flag_bits", "date_time", "external_attr",
                              "create_system", "create_version", "extract_version"):
                setattr(clone, attribute, getattr(info, attribute))
            payload = bundle.read(clone)
            break
        else:  # pragma: no cover - the archive always carries the hull
            raise SystemExit(f"DTC_hull.igs not found inside {archive}")
    with tempfile.TemporaryDirectory() as tmp:
        path = Path(tmp) / "DTC_hull.igs"
        path.write_bytes(payload)
        model = read_iges(str(path))
    usable, _ = partition_surfaces(model.surfaces())
    return np.asarray(list(sample_patches(usable, spacing=spacing)), dtype=float)


def sectional_area(points: np.ndarray, draft: float, dx: float, dz: float) -> tuple:
    """Station areas ``A(x)`` from the repo's half-breadth field.

    The field is the same interpolated ``b[i][k] = max |y|`` grid the fleet build measures with
    (:func:`tensorcad.fleet_baseline.field_from_points`), so the curve plotted here integrates to
    the same displacement the measurement record reports.  Binning the raw sample points
    directly -- the first attempt -- only catches patch vertices and under-reads the area by
    about a third, which is why the field, not the point cloud, is integrated.
    """
    field = field_from_points(points, draft, dx=dx, dz=dz)
    nz = field.nz
    areas = []
    for i in range(field.nx):
        column = [field.b[i][k] for k in range(nz)]
        areas.append(2.0 * dz * (0.5 * column[0] + sum(column[1:-1]) + 0.5 * column[-1]))
    xs = [field.x(i) for i in range(field.nx)]
    volume = sum(0.5 * (areas[i] + areas[i + 1]) * field.dx for i in range(field.nx - 1))
    return np.asarray(xs), np.asarray(areas), volume, field_integrals(field)["volume_m3"]


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--archive", default=str(DEFAULT_ARCHIVE))
    parser.add_argument("--measured", default=str(DEFAULT_MEASURED))
    parser.add_argument("--out", default=str(DEFAULT_OUT))
    parser.add_argument("--dx", type=float, default=0.5)
    parser.add_argument("--dz", type=float, default=0.1)
    args = parser.parse_args(argv)

    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    measured = json.loads(Path(args.measured).read_text())
    reference = measured["published_reference"]
    points = hull_points(Path(args.archive))
    x, areas, volume, field_volume = sectional_area(points, reference["draft"], args.dx,
                                                   args.dz)
    volume = float(volume)
    volume_error = 100.0 * (volume - reference["displacement_volume_m3"]) / reference[
        "displacement_volume_m3"]

    figure, (left, right) = plt.subplots(1, 2, figsize=(12.4, 4.6))
    left.plot(x, areas, color="#1f4e79", lw=1.8)
    left.fill_between(x, 0.0, areas, color="#1f4e79", alpha=0.12)
    left.set_xlabel("x [m] from the publisher's origin")
    left.set_ylabel("sectional area A(x) [m²]")
    left.set_title("DTC sectional area integrated from igs_ordner.zip\n"
                   f"∫A dx = {volume:,.0f} m³ vs published "
                   f"{reference['displacement_volume_m3']:,.0f} m³ "
                   f"({volume_error:+.2f} %)", fontsize=10)
    left.grid(alpha=0.3)

    checks = measured["measured"]["checks"]
    labels, deviations = [], []
    for key, label in (("beam", "B"), ("block_coefficient", "C_B"),
                       ("displacement_volume", "∇")):
        labels.append(label)
        deviations.append(checks[key]["rel_error_percent"])
    positions = np.arange(len(labels))
    bars = right.bar(positions, deviations, color="#c00000", alpha=0.85, width=0.55)
    right.axhline(0.0, color="black", lw=0.8)
    for bound in (-1.0, 1.0):
        right.axhline(bound, color="grey", ls="--", lw=1.0)
    right.text(len(labels) - 0.45, 1.02, "±1 % acceptance", color="grey", fontsize=8, ha="right")
    right.bar_label(bars, fmt="%+.3f %%", fontsize=9, padding=2)
    right.set_xticks(positions, labels)
    right.set_ylim(-1.4, 1.4)
    right.set_ylabel("deviation from published [%]")
    right.set_title("Hull measured from the same archive\n"
                    "(dtc_hull_measured.json, dx 0.5 m / dz 0.1 m)", fontsize=10)
    right.grid(alpha=0.3, axis="y")

    figure.tight_layout()
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(out, dpi=160)
    print(f"wrote {out}")
    print(f"sectional-area integral {volume:.3f} m3 ({volume_error:+.3f} %) "
          f"vs published {reference['displacement_volume_m3']}")
    print(f"field_integrals cross-check {field_volume:.3f} m3 "
          f"({100.0 * (field_volume - volume) / volume:+.3f} % vs the trapezoid above)")
    print("checks:", {label: round(value, 4) for label, value in zip(labels, deviations)})
    return 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())