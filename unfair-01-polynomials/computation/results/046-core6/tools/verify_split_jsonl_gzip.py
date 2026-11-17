#!/usr/bin/env python3
"""Verify the line-preserving terminal-journal split without needing its raw gzip."""

from __future__ import annotations

import gzip
import hashlib
import json
from pathlib import Path


def sha(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def main() -> None:
    here = Path(__file__).resolve().parents[1] / "data/single"
    meta = json.loads((here / "single-terminal-witnesses-split.json").read_text())
    plain = hashlib.sha256()
    lines = 0
    for row in meta["parts"]:
        path = here / row["name"]
        assert path.stat().st_size == row["compressed_bytes"]
        assert sha(path) == row["compressed_sha256"]
        count = 0
        bytes_in_part = 0
        with gzip.open(path, "rb") as f:
            for line in f:
                plain.update(line)
                lines += 1
                count += 1
                bytes_in_part += len(line)
        assert bytes_in_part == row["uncompressed_bytes"]
    assert lines == meta["line_count"]
    assert plain.hexdigest() == meta["source_uncompressed_sha256"]
    print(json.dumps({"status": "passed", "parts": len(meta["parts"]), "line_count": lines}, indent=2))


if __name__ == "__main__":
    main()
