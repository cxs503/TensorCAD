import io, re, sys

root = "/root/tp-agent/work/m3832/"

# ---------- CHANGELOG.md : semantic union under 0.21.0 ----------
p = root + "CHANGELOG.md"
lines = open(p).read().split("\n")
i1 = lines.index("<<<<<<< HEAD")
i2 = lines.index("=======")
i3 = lines.index(">>>>>>> origin/main")
head_block = lines[i1 + 1:i2]
main_block = lines[i2 + 1:i3]
assert head_block and main_block, ("empty side", len(head_block), len(main_block))
assert head_block[0].startswith("- **Righting-arm curve"), head_block[0]
assert main_block[0] == "### Added", main_block[0]
new = lines[:i1] + head_block + [""] + main_block + lines[i3 + 1:]
open(p, "w").write("\n".join(new))

# ---------- tests/test_public_surface.py : union of the module list ----------
p = root + "tests/test_public_surface.py"
lines = open(p).read().split("\n")
i1 = lines.index("<<<<<<< HEAD")
i2 = lines.index("=======")
i3 = lines.index(">>>>>>> origin/main")
head_block = lines[i1 + 1:i2]
main_block = lines[i2 + 1:i3]
assert len(head_block) == 2 and len(main_block) == 1, (head_block, main_block)
assert '"righting_arm"' in head_block[0] and '"bonjean"' in head_block[0], head_block
assert '"transom"' in main_block[0], main_block
new = lines[:i1] + head_block + [main_block[0]] + lines[i3 + 1:]
open(p, "w").write("\n".join(new))

print("CHANGELOG head block lines:", len(head_block), "main block lines:", len(main_block))
print("merged module list line:", main_block[0].strip())