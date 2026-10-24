"""Optional solver, meshing and CAD-kernel adapters.

Nothing in the core package may import from here: the geometry contract in
:mod:`tensorcad.geometry` has to stay usable with only the standard library
installed. Adapters translate that contract to a specific consumer -- a volume
mesher, a BEM/LBM/FEM solver's geometry import, or an OCP/OpenCASCADE B-rep
session -- and each of them lives behind an extra in ``pyproject.toml``.

Layout::

    tensorcad.adapters.fem     # node/element tables + gmsh 2.2 / JSON export (no extra)
    tensorcad.adapters.ocp     # sew -> solid -> STEP, behind the "ocp" extra

Two different promises live here, and conflating them is the usual way "we export
STEP" becomes a claim nobody can check:

* :mod:`tensorcad.adapters.fem` is a *format* adapter. It needs nothing beyond the
  standard library plus NumPy, because the receiving side is a mesh format (gmsh 2.2
  ASCII, JSON), not a kernel.
* :mod:`tensorcad.adapters.ocp` is a *kernel* adapter. It sews the triangles into a
  closed shell, makes a solid, asks OpenCASCADE's own checker whether the solid is
  valid, and writes STEP. The result is a valid faceted solid -- thousands of planar
  faces, one per triangle -- which is what a solver-facing exchange file needs and
  what a designer does *not* want (no fillets, no analytic surfaces, no feature
  tree). The module reports the facts it measured (free edges, shell closure,
  validity, volume, re-imported volume) and makes no claim about being a B-rep model
  of the hull.

Neither module is imported by the core: ``import tensorcad`` works with the standard
library alone, and every adapter is reached explicitly.
"""
from __future__ import annotations

__all__: list[str] = []
