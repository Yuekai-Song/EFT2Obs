#!/usr/bin/env python3
"""Generate clangd commands for Rivet plugins using the local Rivet installation.

Run inside the compiler environment used to build the plugins.
"""

import json
import os
from pathlib import Path
import shlex
import shutil
import subprocess


def main():
    root = Path(__file__).resolve().parents[1]
    plugin_dir = root / "RivetPlugins"
    compiler = shlex.split(os.environ.get("CXX") or "g++")
    compiler_path = shutil.which(compiler[0])
    if compiler_path is None:
        raise SystemExit("C++ compiler unavailable; initialize the build environment first")
    compiler[0] = compiler_path
    flags = shlex.split(subprocess.check_output(
        [str(root / "local/bin/rivet-config"), "--cppflags", "--cxxflags"],
        text=True,
    ))
    # rivet-build also defines NDEBUG and builds position-independent code.
    flags += ["-DNDEBUG", "-fPIC"]
    for variable in ("CPPFLAGS", "CXXFLAGS"):
        flags += shlex.split(os.environ.get(variable, ""))
    sources = sorted(plugin_dir.glob("*.cc"))
    if not sources:
        raise SystemExit("No Rivet plugin sources found")
    entries = [{
        "directory": str(plugin_dir),
        "file": str(source),
        "arguments": compiler + flags + ["-c", str(source)],
    } for source in sources]
    destination = root / "compile_commands.json"
    destination.write_text(json.dumps(entries, indent=2) + "\n")
    print(f"Wrote {len(entries)} plugin commands to {destination}")
    print(f"Compiler: {compiler_path}")


if __name__ == "__main__":
    main()
