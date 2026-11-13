"""The ``(x, w, n, |J|)`` element contract TensorDG/TensorFEM consume.

The contract is deliberately boring: flat arrays of float64 plus a JSON header, with a
SHA-256 over the canonical byte layout. Boring is the point -- a downstream solver needs
to be able to (a) read the arrays without importing this package, (b) check that the file
it read is the file that was written, and (c) recompute the surface integrals from the
arrays and get the number this package reported, bit for bit.

Round-trip identity is therefore a *gate*, not a nicety: the writer stores float64 and
the reader returns float64, so every value must return bit-identical and a recomputed
integral must match to ``0.0`` (not "about 1e-14"). The tests assert ``0.0``; the
campaign budget of ``1e-14`` is what a *format conversion* through text would need, and is
reported alongside as the looser bound.

Formats
-------
* ``.npz`` -- the primary format. NumPy is imported lazily so that the core package keeps
  installing with no third-party dependency at all.
* ``.h5`` -- written only if ``h5py`` is installed; the same arrays under a ``/elements``
  group with the header JSON as an attribute plus a ``header`` dataset.
* ``.json`` -- the portable fallback. Floats are written with ``repr`` (shortest
  round-trip form), so it is also bit-exact for float64, just bigger.

Named limitations
-----------------
* Only the tensor-product element set this package produces: quadrilateral elements in a
  rectangular patch. A solver that needs triangles gets them from the tessellation, not
  from here, and the negative control in the campaign is exactly the statement that
  triangles lose the p-order.
* The content hash covers the *payload* (arrays + header), not the git commit of the
  generator; ``provenance`` carries the commit/version fields and is itself hashed, so a
  hand-edited provenance does change the hash.
* No compression is applied to the NPZ so that the byte layout -- and therefore the file
  digest recorded in the evidence -- is reproducible across zlib versions.
"""

from __future__ import annotations

import hashlib
import json
import struct
from collections.abc import Iterable, Sequence
from pathlib import Path
from typing import Any

from .quadrature import CurvedElement

__all__ = [
    "ELEMENTS_SCHEMA",
    "build_payload",
    "content_hash",
    "read_contract",
    "recompute_integrals",
    "verify_roundtrip",
    "write_contract",
]

#: Schema tag written into every file; bump only on a breaking layout change.
ELEMENTS_SCHEMA = "tensorcad.highorder.elements/1.0"

#: Field order inside the canonical hash payload. Frozen: changing it changes every hash.
_HASH_FIELDS = ("points", "weights", "normals", "jacobian_det", "uv", "element_id")


def _flatten(elements: Sequence[CurvedElement]) -> dict[str, list[Any]]:
    points: list[float] = []
    weights: list[float] = []
    normals: list[float] = []
    jacobians: list[float] = []
    uv: list[float] = []
    element_ids: list[float] = []
    for element in elements:
        element_weights = element.weights
        if not (
            len(element.points) == len(element_weights) == len(element.normals) == len(element.uv)
        ):
            raise ValueError(f"element {element.element_id} has inconsistent point counts")
        for point, weight, normal, jacobian, parameter in zip(
            element.points,
            element_weights,
            element.normals,
            element.jacobian_det,
            element.uv,
            strict=True,
        ):
            points.extend(point)
            weights.append(weight)
            normals.extend(normal)
            jacobians.append(jacobian)
            uv.extend(parameter)
            element_ids.append(float(element.element_id))
    return {
        "points": points,
        "weights": weights,
        "normals": normals,
        "jacobian_det": jacobians,
        "uv": uv,
        "element_id": element_ids,
    }


def build_payload(
    elements: Sequence[CurvedElement],
    *,
    provenance: dict[str, Any] | None = None,
    quantity: str = "surface_elements",
) -> dict[str, Any]:
    """The in-memory contract: flat arrays plus a header describing them."""
    if not elements:
        raise ValueError("the contract needs at least one element")
    degrees = {tuple(element.degree) for element in elements}
    orders = {element.quadrature_order for element in elements}
    if len(degrees) != 1:
        raise ValueError(f"mixed polynomial degrees in one contract: {sorted(degrees)}")
    if len(orders) != 1:
        raise ValueError(f"mixed quadrature orders in one contract: {sorted(orders)}")
    arrays = _flatten(elements)
    degree = next(iter(degrees))
    header: dict[str, Any] = {
        "schema": ELEMENTS_SCHEMA,
        "quantity": quantity,
        "rows": len(arrays["weights"]),
        "elements": len(elements),
        "degree": [degree[0], degree[1]],
        "quadrature_order": next(iter(orders)),
        "layout": {
            "points": ["row", 3],
            "weights": ["row"],
            "normals": ["row", 3],
            "jacobian_det": ["row"],
            "uv": ["row", 2],
            "element_id": ["row"],
        },
        "conventions": {
            "weights": "w_ref * |S_u x S_v| * metric_scale; summing them gives the area",
            "normals": "unit S_u x S_v, outward for an outward-oriented patch",
            "jacobian_det": "|S_u x S_v| (not multiplied by the reference weight)",
            "element_id": "0-based index into the element set the writer received",
        },
        "provenance": dict(provenance or {}),
    }
    header["content_hash"] = content_hash(arrays, header)
    return {"header": header, "arrays": arrays}


def _canonical_bytes(arrays: dict[str, list[Any]], header: dict[str, Any]) -> bytes:
    digest = hashlib.sha256()
    digest.update(b"tensorcad.highorder.elements/1.0\x00")
    digest.update(json.dumps(header, sort_keys=True, separators=(",", ":")).encode("utf-8"))
    for name in _HASH_FIELDS:
        values = arrays[name]
        digest.update(name.encode("ascii"))
        digest.update(struct.pack(f"<{len(values)}d", *[float(value) for value in values]))
    return digest.digest()


def content_hash(arrays: dict[str, list[Any]], header: dict[str, Any]) -> str:
    """SHA-256 over the canonical layout; stable across runs and machines."""
    return _canonical_bytes(arrays, header).hex()


def _without_hash(header: dict[str, Any]) -> dict[str, Any]:
    return {key: value for key, value in header.items() if key != "content_hash"}


def write_contract(
    path: str | Path,
    elements: Sequence[CurvedElement],
    *,
    provenance: dict[str, Any] | None = None,
    quantity: str = "surface_elements",
) -> dict[str, Any]:
    """Write the element set and return its header (including the content hash)."""
    payload = build_payload(elements, provenance=provenance, quantity=quantity)
    target = Path(path)
    header: dict[str, Any] = payload["header"]
    arrays = payload["arrays"]
    target.parent.mkdir(parents=True, exist_ok=True)
    suffix = target.suffix.lower()
    if suffix == ".npz":
        import numpy as np  # imported lazily: the core package stays stdlib-only

        np.savez(
            target,
            meta=np.array(json.dumps(header, sort_keys=True)),
            points=np.array(arrays["points"], dtype=np.float64).reshape(-1, 3),
            weights=np.array(arrays["weights"], dtype=np.float64),
            normals=np.array(arrays["normals"], dtype=np.float64).reshape(-1, 3),
            jacobian_det=np.array(arrays["jacobian_det"], dtype=np.float64),
            uv=np.array(arrays["uv"], dtype=np.float64).reshape(-1, 2),
            element_id=np.array(arrays["element_id"], dtype=np.float64),
        )
    elif suffix == ".json":
        target.write_text(
            json.dumps({"header": header, "arrays": arrays}, sort_keys=True),
            encoding="utf-8",
        )
    elif suffix in (".h5", ".hdf5"):
        import h5py  # type: ignore[import-untyped]
        import numpy as np

        with h5py.File(target, "w") as handle:
            group = handle.create_group("elements")
            group.attrs["header"] = json.dumps(header, sort_keys=True)
            group.create_dataset("meta", data=json.dumps(header, sort_keys=True))
            group.create_dataset(
                "points", data=np.array(arrays["points"], dtype=np.float64).reshape(-1, 3)
            )
            group.create_dataset("weights", data=np.array(arrays["weights"], dtype=np.float64))
            group.create_dataset(
                "normals", data=np.array(arrays["normals"], dtype=np.float64).reshape(-1, 3)
            )
            group.create_dataset(
                "jacobian_det", data=np.array(arrays["jacobian_det"], dtype=np.float64)
            )
            group.create_dataset("uv", data=np.array(arrays["uv"], dtype=np.float64).reshape(-1, 2))
            group.create_dataset(
                "element_id", data=np.array(arrays["element_id"], dtype=np.float64)
            )
    else:
        raise ValueError(f"unsupported contract extension {suffix!r}: use .npz, .json or .h5")
    return header


def _read_json(path: Path) -> dict[str, Any]:
    document: dict[str, Any] = json.loads(path.read_text(encoding="utf-8"))
    return {"header": document["header"], "arrays": document["arrays"]}


def _read_npz(path: Path) -> dict[str, Any]:
    import numpy as np

    with np.load(path, allow_pickle=False) as handle:
        header: dict[str, Any] = json.loads(str(handle["meta"]))
        arrays = {
            "points": [float(value) for value in np.asarray(handle["points"]).reshape(-1)],
            "weights": [float(value) for value in np.asarray(handle["weights"]).reshape(-1)],
            "normals": [float(value) for value in np.asarray(handle["normals"]).reshape(-1)],
            "jacobian_det": [float(value) for value in np.asarray(handle["jacobian_det"]).reshape(-1)],
            "uv": [float(value) for value in np.asarray(handle["uv"]).reshape(-1)],
            "element_id": [float(value) for value in np.asarray(handle["element_id"]).reshape(-1)],
        }
    return {"header": header, "arrays": arrays}


def _text(value: Any) -> str:
    """HDF5 hands back bytes for a JSON string dataset; decode it explicitly."""
    return value.decode("utf-8") if isinstance(value, bytes) else str(value)


def _read_hdf5(path: Path) -> dict[str, Any]:
    import h5py
    import numpy as np

    with h5py.File(path, "r") as handle:
        group = handle["elements"]
        header: dict[str, Any] = json.loads(_text(group["meta"][()]))
        arrays = {
            "points": [float(value) for value in np.asarray(group["points"]).reshape(-1)],
            "weights": [float(value) for value in np.asarray(group["weights"]).reshape(-1)],
            "normals": [float(value) for value in np.asarray(group["normals"]).reshape(-1)],
            "jacobian_det": [float(value) for value in np.asarray(group["jacobian_det"]).reshape(-1)],
            "uv": [float(value) for value in np.asarray(group["uv"]).reshape(-1)],
            "element_id": [float(value) for value in np.asarray(group["element_id"]).reshape(-1)],
        }
    return {"header": header, "arrays": arrays}


def read_contract(path: str | Path) -> dict[str, Any]:
    """Read any of the three formats back into the in-memory contract."""
    target = Path(path)
    suffix = target.suffix.lower()
    if suffix == ".npz":
        return _read_npz(target)
    if suffix == ".json":
        return _read_json(target)
    if suffix in (".h5", ".hdf5"):
        return _read_hdf5(target)
    raise ValueError(f"unsupported contract extension {suffix!r}: use .npz, .json or .h5")


def recompute_integrals(arrays: dict[str, list[Any]]) -> dict[str, float]:
    """Area and divergence-theorem volume implied by the arrays alone.

    This is what a solver does with the file: no surface, no patch, just points, weights
    and normals. Reproducing :func:`tensorcad.highorder.quadrature.surface_area` from the
    arrays is the operational meaning of "the contract round-trips".
    """
    weights = arrays["weights"]
    points = arrays["points"]
    normals = arrays["normals"]
    area = sum(float(value) for value in weights)
    volume = 0.0
    for index, weight in enumerate(weights):
        px, py, pz = points[3 * index], points[3 * index + 1], points[3 * index + 2]
        nx, ny, nz = normals[3 * index], normals[3 * index + 1], normals[3 * index + 2]
        volume += float(weight) * (px * nx + py * ny + pz * nz) / 3.0
    return {"area": area, "volume": volume}


def verify_roundtrip(path: str | Path, original: dict[str, Any]) -> dict[str, Any]:
    """Write-then-read comparison: max absolute deviation per array and the hash verdict."""
    reread = read_contract(path)
    fields: dict[str, dict[str, float]] = {}
    for name in _HASH_FIELDS:
        left = [float(value) for value in original["arrays"][name]]
        right = [float(value) for value in reread["arrays"][name]]
        if len(left) != len(right):
            raise ValueError(f"array {name!r} changed length across the round trip")
        fields[name] = {
            "max_abs_difference": max(
                (abs(a - b) for a, b in zip(left, right, strict=True)), default=0.0
            ),
            "bit_identical": left == right,
        }
    original_hash = original["header"]["content_hash"]
    reread_hash = content_hash(reread["arrays"], _without_hash(reread["header"]))
    first = recompute_integrals(original["arrays"])
    second = recompute_integrals(reread["arrays"])
    return {
        "path": str(path),
        "schema": reread["header"]["schema"],
        "rows": reread["header"]["rows"],
        "fields": fields,
        "hash_match": original_hash == reread_hash == reread["header"]["content_hash"],
        "content_hash": original_hash,
        "recomputed": first,
        "recomputed_reread": second,
        "recomputed_max_abs_difference": max(
            abs(first[key] - second[key]) for key in first
        ),
        "all_bit_identical": all(entry["bit_identical"] for entry in fields.values()),
    }


def iter_arrays(payload: dict[str, Any]) -> Iterable[tuple[str, list[Any]]]:
    """Iterate the payload arrays in canonical order (convenience for writers)."""
    for name in _HASH_FIELDS:
        yield name, payload["arrays"][name]
