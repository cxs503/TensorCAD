"""Optional solver, meshing and CAD-kernel adapters.

Nothing in the core package may import from here: the geometry contract in
:mod:`tensorcad.geometry` has to stay usable with only the standard library
installed. Adapters translate that contract to a specific consumer -- a volume
mesher, a BEM/LBM/FEM solver's geometry import, or an OCP/OpenCASCADE B-rep
session -- and each of them lives behind an extra in ``pyproject.toml``.

Planned layout (see issue #4 for the OCP/STEP work)::

    tensorcad.adapters.ocp     # B-rep build, kernel validity check, STEP round-trip
    tensorcad.adapters.gmsh    # surface -> volume mesh handoff
"""
from __future__ import annotations

__all__: list[str] = []
