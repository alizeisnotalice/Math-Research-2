#!/usr/bin/env python3
"""Check byte sizes and hashes listed in SHA256SUMS.json."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def sha(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def main() -> None:
    rows = json.loads((ROOT / "SHA256SUMS.json").read_text())["files"]
    bad = [r["path"] for r in rows if not (ROOT / r["path"]).is_file()
           or (ROOT / r["path"]).stat().st_size != r["bytes"]
           or sha(ROOT / r["path"]) != r["sha256"]]
    if bad:
        raise SystemExit("bad package files: " + ", ".join(bad))
    print(json.dumps({"status": "passed", "files": len(rows)}, indent=2))


if __name__ == "__main__":
    main()
