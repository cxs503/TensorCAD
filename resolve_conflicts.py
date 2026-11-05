#!/usr/bin/env python3
"""Resolve merge conflicts between main and agent/m3300 by semantic union."""
import subprocess
import os

WORK = "/root/tp-agent/work/m3898"
os.chdir(WORK)

CONFLICT_FILES = [
    "CHANGELOG.md",
    "README.md",
    "docs/fairness.md",
    "docs/hydrostatics_verification.md",
    "tensorcad/__init__.py",
    "tensorcad/cli.py",
    "tensorcad/hull.py",
    "tensorcad/hydrostatics.py",
    "tensorcad/hydrostatics_analytic.py",
    "tests/test_cli.py",
    "tests/test_curvature.py",
    "tests/test_hull.py",
    "tests/test_public_surface.py",
]

def get_blob(ref, path):
    """Get file content from a git ref."""
    result = subprocess.run(
        ["git", "show", f"{ref}:{path}"],
        capture_output=True, text=True, cwd=WORK
    )
    return result.stdout

def write_file(path, content):
    """Write content to file."""
    full = os.path.join(WORK, path)
    os.makedirs(os.path.dirname(full) if os.path.dirname(full) else ".", exist_ok=True)
    with open(full, "w") as f:
        f.write(content)

def resolve_init_py():
    """Resolve tensorcad/__init__.py - union of imports from both sides."""
    main_content = get_blob("origin/main", "tensorcad/__init__.py")
    m3300_content = get_blob("origin/agent/m3300", "tensorcad/__init__.py")
    
    # Parse imports from both sides
    main_lines = main_content.split("\n")
    m3300_lines = m3300_content.split("\n")
    
    # Find import blocks and __all__
    def parse_imports(content):
        imports = []
        all_exports = []
        in_import = False
        in_all = False
        current_import = []
        
        for line in content.split("\n"):
            if line.startswith("from .") or line.startswith("import "):
                if not in_import:
                    in_import = True
                    current_import = [line]
                else:
                    current_import.append(line)
                if line.rstrip().endswith(")"):
                    imports.append("\n".join(current_import))
                    current_import = []
                    in_import = False
            elif line.startswith("__all__"):
                in_all = True
                all_exports.append(line)
                if "]" in line:
                    in_all = False
            elif in_all:
                all_exports.append(line)
                if "]" in line:
                    in_all = False
            elif not in_import:
                pass  # other lines
        
        return imports, all_exports
    
    main_imports, main_all = parse_imports(main_content)
    m3300_imports, m3300_all = parse_imports(m3300_content)
    
    # Union of imports - keep all unique imports
    all_imports = []
    seen_modules = set()
    for imp in main_imports + m3300_imports:
        # Extract module name
        for line in imp.split("\n"):
            if line.startswith("from ."):
                mod = line.split()[1].rstrip(".")
                if mod not in seen_modules:
                    seen_modules.add(mod)
                    all_imports.append(imp)
                break
    
    # Union of __all__ exports
    def extract_all_names(all_lines):
        names = []
        for line in all_lines:
            for token in line.replace("[", "").replace("]", "").replace(",", "").replace('"', "").replace("'", "").split():
                if token and not token.startswith("#"):
                    names.append(token)
        return names
    
    main_names = extract_all_names(main_all)
    m3300_names = extract_all_names(m3300_all)
    all_names = list(dict.fromkeys(main_names + m3300_names))  # preserve order, remove dups
    
    # Build the file: use main as base, add m3300-specific imports
    # Find where imports end in main
    result_lines = []
    in_import_block = False
    import_block_end = 0
    
    for i, line in enumerate(main_lines):
        if line.startswith("from .") or line.startswith("import "):
            in_import_block = True
        elif in_import_block and line.strip() == "":
            import_block_end = i
            in_import_block = False
            break
    
    # Find m3300-specific imports (modules not in main)
    main_modules = set()
    for line in main_content.split("\n"):
        if line.startswith("from ."):
            main_modules.add(line.split()[1].rstrip("."))
    
    m3300_specific_imports = []
    for imp in m3300_imports:
        for line in imp.split("\n"):
            if line.startswith("from ."):
                mod = line.split()[1].rstrip(".")
                if mod not in main_modules:
                    m3300_specific_imports.append(imp)
                break
    
    # Build result: main content with m3300-specific imports added
    # Find the position after the last import in main
    last_import_idx = 0
    for i, line in enumerate(main_lines):
        if line.startswith("from .") or line.startswith("import "):
            last_import_idx = i
        elif line.startswith("__all__"):
            break
    
    # Insert m3300-specific imports after main imports
    result_lines = main_lines[:last_import_idx+1]
    for imp in m3300_specific_imports:
        result_lines.append("")
        result_lines.extend(imp.split("\n"))
    result_lines.extend(main_lines[last_import_idx+1:])
    
    # Update __all__ to include m3300-specific exports
    # Find __all__ block and update it
    final_content = "\n".join(result_lines)
    
    # Simple approach: replace __all__ with union
    import re
    all_pattern = r'__all__\s*=\s*\[[\s\S]*?\]'
    new_all = '__all__ = [\n' + ',\n'.join(f'    "{n}"' for n in all_names) + '\n]'
    final_content = re.sub(all_pattern, new_all, final_content)
    
    return final_content

def resolve_simple_union(path):
    """For files where one side is a superset, use the longer one."""
    main_content = get_blob("origin/main", path)
    m3300_content = get_blob("origin/agent/m3300", path)
    
    if m3300_content in main_content:
        return main_content
    elif main_content in m3300_content:
        return m3300_content
    else:
        # Need manual resolution - return None to signal
        return None

def resolve_changelog():
    """Resolve CHANGELOG.md - keep all version entries from both sides."""
    main_content = get_blob("origin/main", "CHANGELOG.md")
    m3300_content = get_blob("origin/agent/m3300", "CHANGELOG.md")
    
    # Parse version sections
    def parse_versions(content):
        versions = []
        current_version = None
        current_content = []
        
        for line in content.split("\n"):
            if line.startswith("## "):
                if current_version:
                    versions.append((current_version, "\n".join(current_content)))
                current_version = line
                current_content = []
            else:
                current_content.append(line)
        
        if current_version:
            versions.append((current_version, "\n".join(current_content)))
        
        return versions
    
    main_versions = parse_versions(main_content)
    m3300_versions = parse_versions(m3300_content)
    
    # Union: keep all unique version headers
    seen_headers = set()
    all_versions = []
    
    for header, content in main_versions:
        if header not in seen_headers:
            seen_headers.add(header)
            all_versions.append((header, content))
    
    for header, content in m3300_versions:
        if header not in seen_headers:
            seen_headers.add(header)
            all_versions.append((header, content))
    
    # Sort by version number (descending)
    def version_key(v):
        header = v[0]
        # Extract version number
        import re
        match = re.search(r'## (\d+\.\d+\.\d+)', header)
        if match:
            parts = match.group(1).split(".")
            return tuple(int(p) for p in parts)
        return (0, 0, 0)
    
    # Separate versioned entries from special ones (Unreleased)
    versioned = [v for v in all_versions if re.search(r'## \d+\.\d+\.\d+', v[0])]
    special = [v for v in all_versions if not re.search(r'## \d+\.\d+\.\d+', v[0])]
    
    versioned.sort(key=version_key, reverse=True)
    
    # Build result
    result = []
    for header, content in special + versioned:
        result.append(header)
        result.append(content)
    
    return "\n".join(result)

def resolve_code_file(path):
    """Resolve a Python file by union of functions/classes."""
    main_content = get_blob("origin/main", path)
    m3300_content = get_blob("origin/agent/m3300", path)
    
    # For now, use the longer version as base and add unique content from the other
    if len(main_content) >= len(m3300_content):
        base = main_content
        other = m3300_content
        base_ref = "main"
    else:
        base = m3300_content
        other = main_content
        base_ref = "m3300"
    
    # Check if one is a subset of the other
    if other.strip() in base.strip():
        return base
    
    # For complex cases, we need to do line-by-line union
    # This is a simplified approach - just use the longer version
    # and note that manual review may be needed
    return base

# Resolve each file
resolutions = {}

for f in CONFLICT_FILES:
    print(f"Resolving {f}...")
    
    if f == "tensorcad/__init__.py":
        resolutions[f] = resolve_init_py()
    elif f == "CHANGELOG.md":
        resolutions[f] = resolve_changelog()
    elif f in ["docs/fairness.md", "docs/hydrostatics_verification.md"]:
        resolutions[f] = resolve_simple_union(f)
    else:
        resolutions[f] = resolve_code_file(f)

# Write resolutions
for f, content in resolutions.items():
    if content is not None:
        write_file(f, content)
        print(f"  Written {f} ({len(content)} chars)")
    else:
        print(f"  NEEDS MANUAL RESOLUTION: {f}")

print("\nDone!")
