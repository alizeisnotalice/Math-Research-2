#!/usr/bin/env python3
"""Line-preserving archival merge for the 680 independent AB record files."""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
from pathlib import Path


def sha(raw: bytes) -> str:
    return hashlib.sha256(raw).hexdigest()


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("records", type=Path)
    p.add_argument("out", type=Path)
    args = p.parse_args()
    paths = sorted(args.records.glob("*.json"))
    if len(paths) != 680:
        raise SystemExit(f"expected 680 records, found {len(paths)}")
    args.out.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(args.out, "wt", encoding="utf-8", newline="") as out:
        for pth in paths:
            raw = pth.read_bytes()
            rec = json.loads(raw)
            row = {
                "index": rec["index"],
                "original_name": pth.name,
                "original_sha256": sha(raw),
                "original_utf8": raw.decode("utf-8"),
            }
            out.write(json.dumps(row, separators=(",", ":"), ensure_ascii=False) + "\n")
    print(json.dumps({"records": len(paths), "out_sha256": sha(args.out.read_bytes())}, indent=2))


if __name__ == "__main__":
    main()
