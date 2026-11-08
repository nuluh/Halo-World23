#!/usr/bin/env python3
"""Verify that supported headers do not pull experimental/support APIs transitively."""
from __future__ import annotations

import argparse
import pathlib
import subprocess
import tempfile

FORBIDDEN = (
    "openvfife/material/experimental.hpp",
    "openvfife/section/experimental.hpp",
    "openvfife/section/grid_section.hpp",
    "openvfife/element/experimental.hpp",
    "headers/config.hpp",
    "examples/support/config.hpp",
)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--compiler", required=True)
    parser.add_argument("--source-dir", required=True)
    parser.add_argument("--header", required=True)
    args = parser.parse_args()

    source_dir = pathlib.Path(args.source_dir).resolve()
    with tempfile.TemporaryDirectory() as tmp:
        source = pathlib.Path(tmp) / "surface.cpp"
        source.write_text(f"#include <{args.header}>\n", encoding="utf-8")
        command = [
            args.compiler,
            "-std=c++20",
            f"-I{source_dir / 'include'}",
            f"-I{source_dir / 'headers'}",
            "-M",
            str(source),
        ]
        result = subprocess.run(command, check=True, capture_output=True, text=True)

    dependencies = result.stdout.replace("\\\n", " ")
    leaks = [entry for entry in FORBIDDEN if entry in dependencies]
    if leaks:
        raise SystemExit(
            f"{args.header} leaks unsupported headers transitively: " + ", ".join(leaks)
        )
    print(f"supported include surface ok: {args.header}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
