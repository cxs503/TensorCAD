#!/usr/bin/env python3
"""Build the Verdict oracle shared library used by ``tensorcad.meshquality``.

What it builds
--------------
``tools/verdict_oracle/verdict_oracle.cpp`` plus the Verdict element-metric sources from
VTK's ``ThirdParty/verdict/vtkverdict`` directory, into ``build/libverdict_oracle.so``.
The wrapper adds no metric algebra (see the .cpp header); the numbers are Verdict's, which
are also the numbers ``vtkMeshQuality`` -- and therefore ``meshmagick.print_quality()`` --
reports.

Two ways to supply the Verdict sources:

* ``--source DIR``  an existing ``.../ThirdParty/verdict/vtkverdict`` checkout;
* ``--fetch`` (default when ``--source`` is absent) a sparse, blob-filtered clone of the
  VTK repository at :data:`VTK_REF`, so only ``ThirdParty/verdict`` is transferred.

The build writes ``build/provenance.json``: the VTK reference and commit the sources came
from, the compiler, the flags, the source file digests and the digest of the library.  That
file is what the documentation quotes, so a reader can tell which Verdict they are looking
at instead of "some Verdict".

Usage::

    python tools/verdict_oracle/build.py                     # fetch VTK_REF, build
    python tools/verdict_oracle/build.py --source /path/to/vtkverdict
    python tools/verdict_oracle/build.py --check             # is a library present?

Exit status is 0 when a library exists afterwards, 1 otherwise (CI-friendly: a missing
compiler is reported, not swallowed).
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import platform
import shutil
import subprocess
import sys
from pathlib import Path

#: VTK release the Verdict sources are taken from (see tools/verdict_oracle/verdict_oracle.cpp).
VTK_REF = "v9.3.1"
VTK_URL = "https://github.com/Kitware/VTK.git"
VTK_SUBDIR = "ThirdParty/verdict/vtkverdict"

HERE = Path(__file__).resolve().parent
REPO_ROOT = HERE.parent.parent
BUILD_DIR = HERE / "build"
LIBRARY_NAME = "libverdict_oracle.so"

#: Verdict sources; every element kind the wrapper exposes needs its own translation unit.
VERDICT_SOURCES = (
    "V_EdgeMetric.cpp",
    "V_GaussIntegration.cpp",
    "V_HexMetric.cpp",
    "V_KnifeMetric.cpp",
    "V_PyramidMetric.cpp",
    "V_QuadMetric.cpp",
    "V_TetMetric.cpp",
    "V_TriMetric.cpp",
    "V_WedgeMetric.cpp",
    "VerdictVector.cpp",
)

CXXFLAGS = ("-O2", "-std=c++17", "-fPIC", "-shared", "-fvisibility=default")


def sha256_of(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def fetch_sources(cache_dir: Path, ref: str) -> Path:
    """Sparse, blob-filtered clone of VTK at ``ref``; returns the vtkverdict directory."""
    source = cache_dir / f"VTK-{ref}" / VTK_SUBDIR
    if (source / "verdict.h").is_file():
        print(f"sources: reusing {source}")
        return source
    cache_dir.mkdir(parents=True, exist_ok=True)
    checkout = cache_dir / f"VTK-{ref}"
    if checkout.exists():
        shutil.rmtree(checkout)
    _run(
        ["git", "clone", "--depth", "1", "--branch", ref, "--filter=blob:none", "--sparse",
         VTK_URL, str(checkout)]
    )
    _run(["git", "sparse-checkout", "set", VTK_SUBDIR], cwd=checkout)
    if not (source / "verdict.h").is_file():
        raise SystemExit(f"fetched {ref} but {source / 'verdict.h'} is missing")
    return source


def _run(cmd: list[str], cwd: Path | None = None) -> None:
    print("+", " ".join(cmd))
    subprocess.run(cmd, cwd=cwd, check=True)


def write_config(source: Path, build_dir: Path) -> Path:
    """Instantiate ``verdict_config.h.in`` without CMake.

    The template only needs ``VERDICT_VERSION``; the two ``#cmakedefine`` blocks stay
    undefined, which is the static, unmangled configuration (the wrapper wants the plain
    ``verdict::`` names).
    """
    template = (source / "verdict_config.h.in").read_text()
    body = template.replace("@verdict_VERSION_FLAT@", "0")
    body = body.replace("#cmakedefine BUILD_SHARED_LIBS", "/* #undef BUILD_SHARED_LIBS */")
    body = body.replace("#cmakedefine VERDICT_MANGLE", "/* #undef VERDICT_MANGLE */")
    target = build_dir / "verdict_config.h"
    target.write_text(body)
    return target


def compiler_version(cxx: str) -> str:
    try:
        out = subprocess.run([cxx, "--version"], capture_output=True, text=True, check=True)
    except (OSError, subprocess.CalledProcessError) as exc:  # pragma: no cover - host specifics
        return f"unavailable ({exc})"
    return out.stdout.splitlines()[0].strip()


def git_revision(checkout: Path) -> str | None:
    try:
        out = subprocess.run(
            ["git", "rev-parse", "HEAD"], cwd=checkout, capture_output=True, text=True, check=True
        )
    except (OSError, subprocess.CalledProcessError):
        return None
    return out.stdout.strip() or None


def build(source: Path, cxx: str, build_dir: Path = BUILD_DIR) -> Path:
    build_dir.mkdir(parents=True, exist_ok=True)
    config = write_config(source, build_dir)
    missing = [name for name in VERDICT_SOURCES if not (source / name).is_file()]
    if missing:
        raise SystemExit(f"source tree {source} is missing {', '.join(missing)}")
    wrapper = HERE / "verdict_oracle.cpp"
    library = build_dir / LIBRARY_NAME
    cmd = [cxx, *CXXFLAGS, "-I", str(source), "-I", str(build_dir), "-o", str(library), str(wrapper)]
    cmd += [str(source / name) for name in VERDICT_SOURCES]
    _run(cmd)

    sources = {name: sha256_of(source / name) for name in VERDICT_SOURCES}
    checkout = source.parents[1] if (source.parents[1] / ".git").exists() else None
    provenance = {
        "library": str(library.relative_to(REPO_ROOT)),
        "library_sha256": sha256_of(library),
        "wrapper": str(wrapper.relative_to(REPO_ROOT)),
        "wrapper_sha256": sha256_of(wrapper),
        "verdict_source_dir": str(source),
        "vtk_ref": VTK_REF if VTK_REF in str(source) else None,
        "vtk_commit": git_revision(checkout) if checkout else None,
        "verdict_sources": sources,
        "verdict_config": str(config),
        "compiler": compiler_version(cxx),
        "cxxflags": list(CXXFLAGS),
        "platform": f"{platform.system()} {platform.machine()}",
    }
    (build_dir / "provenance.json").write_text(json.dumps(provenance, indent=2) + "\n")
    print(f"built {library} ({os.path.getsize(library)} bytes)")
    print(f"provenance: {build_dir / 'provenance.json'}")
    return library


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--source", type=Path, default=None,
                        help=f"existing {VTK_SUBDIR} directory (skips the clone)")
    parser.add_argument("--cache", type=Path, default=Path.home() / ".cache" / "tensorcad-vtk",
                        help="where the sparse VTK checkout is kept")
    parser.add_argument("--ref", default=VTK_REF, help=f"VTK tag/branch to fetch (default {VTK_REF})")
    parser.add_argument("--cxx", default=os.environ.get("CXX", "g++"))
    parser.add_argument("--check", action="store_true", help="only report whether a library exists")
    args = parser.parse_args(argv)

    library = BUILD_DIR / LIBRARY_NAME
    if args.check:
        print(f"{library}: {'present' if library.is_file() else 'absent'}")
        return 0 if library.is_file() else 1
    if shutil.which(args.cxx) is None:
        print(f"error: C++ compiler {args.cxx!r} not found; the oracle needs one", file=sys.stderr)
        return 1
    source = args.source
    if source is None:
        try:
            source = fetch_sources(args.cache, args.ref)
        except (OSError, subprocess.CalledProcessError) as exc:
            print(f"error: could not fetch VTK {args.ref}: {exc}", file=sys.stderr)
            return 1
    if not (source / "verdict.h").is_file():
        print(f"error: {source} does not look like {VTK_SUBDIR}", file=sys.stderr)
        return 1
    build(source, args.cxx)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())