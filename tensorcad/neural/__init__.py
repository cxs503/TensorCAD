"""Optional neural geometry runtime (requires the ``neural`` extra).

Two submodules with different dependency footprints:

* :mod:`tensorcad.neural.extract` needs only NumPy — it turns any implicit-field
  callable into a closed triangle mesh;
* :mod:`tensorcad.neural.field` and :mod:`tensorcad.neural.objects` need torch.

This ``__init__`` deliberately imports neither, so a consumer that only needs
extraction can install NumPy alone and still ``import tensorcad.neural.extract``.
Importing :mod:`tensorcad` never imports this package at all: the core geometry
layer stays standard-library only.
"""
from __future__ import annotations

__all__: list[str] = []
