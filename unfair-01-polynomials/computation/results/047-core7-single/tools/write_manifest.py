#!/usr/bin/env python3
"""Write a deterministic SHA-256 manifest for this staging package."""

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
        if p.is_file() and p != OUT and p.name != ".DS_Store":
            rows.append({"path": p.relative_to(ROOT).as_posix(), "bytes": p.stat().st_size, "sha256": sha(p)})
    OUT.write_text(json.dumps({"format": "SHA-256 package manifest; manifest excludes itself", "files": rows}, indent=2) + "\n")
    print(json.dumps({"files": len(rows), "manifest": str(OUT)}, indent=2))


if __name__ == "__main__":
    main()
