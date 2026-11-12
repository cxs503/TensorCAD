"""Smoke probe for tensorcad.mesh_repair (not part of the suite)."""
from __future__ import annotations

from tensorcad.geometry import Mesh
from tensorcad.mesh_repair import (
    fill_holes,
    manifoldize,
    measure_mesh,
    orient_consistently,
    repair_mesh,
)


def box_mesh(n: int = 2, size=(1.0, 1.0, 1.0)) -> Mesh:
    """A closed box, each face split into n x n quads (2 n^2 triangles)."""
    lx, ly, lz = size
    vertices: list[tuple[float, float, float]] = []
    lookup: dict[tuple[int, int, int], int] = {}

    def vid(i, j, k):
        key = (i, j, k)
        if key not in lookup:
            lookup[key] = len(vertices)
            vertices.append((lx * i / n, ly * j / n, lz * k / n))
        return lookup[key]

    faces = []
    for axis in range(3):
        for sign in (0, 1):
            u, v = [a for a in range(3) if a != axis]
            for a in range(n):
                for b in range(n):
                    corners = []
                    for da, db in ((0, 0), (1, 0), (1, 1), (0, 1)):
                        index = [0, 0, 0]
                        index[axis] = sign * n
                        index[u] = a + da
                        index[v] = b + db
                        corners.append(vid(*index))
                    # two triangles, wound so the normal points outward
                    if sign:
                        faces.append((corners[0], corners[1], corners[2]))
                        faces.append((corners[0], corners[2], corners[3]))
                    else:
                        faces.append((corners[0], corners[2], corners[1]))
                        faces.append((corners[0], corners[3], corners[2]))
    return Mesh(vertices=tuple(vertices), faces=tuple(faces))


def main():
    box = box_mesh(2)
    print("box faces", len(box.faces), "volume", round(measure_mesh(box).signed_volume, 12))
    print("clean:", measure_mesh(box).to_dict())

    # remove one triangle -> hole of 3 edges
    holed = Mesh(vertices=box.vertices, faces=box.faces[:-1])
    m = measure_mesh(holed)
    print("\nhole:", {k: m.to_dict()[k] for k in ("open_edges", "faces", "closed")})
    fill = fill_holes(holed)
    print("fill:", fill.to_dict()["filled"], [dict(h) for h in fill.holes])
    print("after fill volume", round(measure_mesh(fill.mesh).signed_volume, 12),
          "open", measure_mesh(fill.mesh).open_edges)

    # remove a full quad (2 triangles) -> 4-edge loop
    holed4 = Mesh(vertices=box.vertices, faces=box.faces[:2] + box.faces[4:])
    fill4 = fill_holes(holed4)
    print("\nquad hole:", fill4.to_dict()["filled"], [dict(h) for h in fill4.holes],
          "volume", round(measure_mesh(fill4.mesh).signed_volume, 12))

    # reversed face -> orientation conflict
    bad = Mesh(vertices=box.vertices, faces=tuple([(f[0], f[2], f[1])] + list(box.faces[1:])))
    print("\nreversed face conflicts:", measure_mesh(bad).winding_conflicts)
    fixed, step = orient_consistently(bad)
    print("orient step:", step.to_dict(), "volume", round(measure_mesh(fixed).signed_volume, 12))

    # whole mesh inward
    inward = Mesh(vertices=box.vertices, faces=tuple((f[0], f[2], f[1]) for f in box.faces))
    fixed2, step2 = orient_consistently(inward)
    print("\ninward mesh volume before", round(measure_mesh(inward).signed_volume, 12),
          "step", step2.to_dict(), "after", round(measure_mesh(fixed2).signed_volume, 12))

    # duplicate face -> non-manifold
    dup = Mesh(vertices=box.vertices, faces=box.faces + (box.faces[0],))
    print("\nduplicate:", {k: measure_mesh(dup).to_dict()[k] for k in
                            ("nonmanifold_edges", "open_edges", "duplicate_faces")})
    rep = repair_mesh(dup, reference_volume=1.0)
    print("repair dup accepted", rep.accepted, "|", rep.reason)
    print("  after", {k: rep.after.to_dict()[k] for k in ("open_edges", "nonmanifold_edges", "faces")},
          "vol", round(rep.after.signed_volume, 12), "rel change", rep.volume_relative_change)

    # full chain on a mesh with a hole
    rep2 = repair_mesh(holed4, reference_volume=1.0)
    print("\nrepair holed accepted", rep2.accepted, "|", rep2.reason)
    print("  volume rel change", rep2.volume_relative_change, "faces",
          rep2.before.faces, "->", rep2.after.faces)


if __name__ == "__main__":
    main()
