#!/usr/bin/env python3
"""Write the selected-file and nonrecursive SHA-256 delivery manifests."""
from __future__ import annotations

import argparse
from hashlib import sha256
from pathlib import Path


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", type=Path, default=Path(__file__).resolve().parents[1])
    ns = ap.parse_args()
    root = ns.root.resolve()
    selected_name = "CERTIFICATE_SELECTED_PATHS-r034.txt"
    sums_name = "DELIVERY_SHA256SUMS-r034.txt"
    excluded = {selected_name, sums_name}
    paths = sorted(p.relative_to(root).as_posix() for p in root.rglob("*") if p.is_file() and p.name not in excluded)
    # The selected-path manifest lists both manifests; the checksum manifest
    # deliberately omits only itself to avoid a recursive hash.
    paths.extend((selected_name, sums_name))
    (root / selected_name).write_text("\n".join(paths) + "\n")
    sums = []
    for rel in paths:
        if rel == sums_name:
            continue
        p = root / rel
        sums.append(f"{sha256(p.read_bytes()).hexdigest()}  {rel}")
    (root / sums_name).write_text("\n".join(sums) + "\n")


if __name__ == "__main__":
    main()
