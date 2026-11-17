#!/usr/bin/env python3
"""Create the deterministic package-file checksum manifest (excluding itself)."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "SHA256SUMS.json"


def sha(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def main() -> None:
    rows = []
    for p in sorted(ROOT.rglob("*")):
        if not p.is_file() or p == OUT or p.name == ".DS_Store":
            continue
        rows.append({"path": p.relative_to(ROOT).as_posix(), "bytes": p.stat().st_size, "sha256": sha(p)})
    payload = {
        "format": "sha256/byte manifest for round046-core6 staging package; manifest is excluded from its own file list",
        "files": rows,
    }
    OUT.write_text(json.dumps(payload, indent=2) + "\n")
    print(f"wrote {OUT}: {len(rows)} files")


if __name__ == "__main__":
    main()
