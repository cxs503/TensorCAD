"""High-order (p-version) geometry for the Tensor ecosystem.

Three pieces, each with a measurable accuracy claim rather than a promised one:

* :mod:`tensorcad.highorder.quadrature` -- Gauss-Legendre tensor-product rules over
  NURBS knot spans, the measured convergence order of the integrals they produce, the
  isoparametric curved element they feed to a solver, and the *faceted* integrator that
  serves as the negative control.
* :mod:`tensorcad.highorder.nurbs_geom` -- first and second fundamental forms, unit
  normals, principal curvatures on a NURBS patch (via exact Bezier-span extraction, so
  second derivatives are analytic rather than finite-differenced), plus the closed-form
  curvature of an ellipsoid used as the analytic truth.
* :mod:`tensorcad.highorder.elliptic` -- AGM/Carlson elliptic integrals and the Legendre
  formula for a triaxial ellipsoid's area: an implementation that shares no code with the
  NURBS route, which is what makes the cross-check in
  ``benchmarks/verified/cad_highorder_geometry`` evidence instead of a tautology.
* :mod:`tensorcad.highorder.contract` -- the ``(x, w, n, |J|)`` element contract with a
  content hash, written as NPZ/HDF5/JSON so TensorDG/TensorFEM can consume it verbatim.

Only the standard library is required by the first three; NumPy/h5py are imported lazily
inside the contract writers so that the core package keeps installing with nothing but
CPython.
"""

from __future__ import annotations

from . import contract, elliptic, nurbs_geom, quadrature

__all__ = ["contract", "elliptic", "nurbs_geom", "quadrature"]
