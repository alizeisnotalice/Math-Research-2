#!/usr/bin/env python3
"""Restore exact original AB record files from the archived merged gzip."""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "data/independent/s1-ab-records-680.jsonl.gz"


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("out", type=Path, help="absent or empty output directory")
    args = p.parse_args()
    if args.out.exists() and any(args.out.iterdir()):
        raise SystemExit("output directory must be absent or empty")
    args.out.mkdir(parents=True, exist_ok=True)
    count = 0
    with gzip.open(SOURCE, "rt", encoding="utf-8") as f:
        for line in f:
            row = json.loads(line)
            raw = row["original_utf8"].encode("utf-8")
            assert hashlib.sha256(raw).hexdigest() == row["original_sha256"]
            target = args.out / row["original_name"]
            target.write_bytes(raw)
            count += 1
    assert count == 680
    print(json.dumps({"status": "extracted", "records": count, "out": str(args.out)}, indent=2))


if __name__ == "__main__":
    main()
