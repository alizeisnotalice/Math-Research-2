#!/usr/bin/env python3
"""Verify every raw independent AB record preserved in the merged gzip."""

from __future__ import annotations

import gzip
import hashlib
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
PATH = ROOT / "data/independent/s1-ab-records-680.jsonl.gz"


def main() -> None:
    count = 0
    source_hash = None
    with gzip.open(PATH, "rt", encoding="utf-8") as f:
        for line in f:
            row = json.loads(line)
            raw = row["original_utf8"].encode("utf-8")
            assert hashlib.sha256(raw).hexdigest() == row["original_sha256"]
            rec = json.loads(raw)
            assert rec["index"] == row["index"] == count
            assert row["original_name"] == f"{count:04}.json"
            assert rec["status"] == "completed" and rec["returncode"] == 0
            source_hash = rec["source_sha256"] if source_hash is None else source_hash
            assert rec["source_sha256"] == source_hash
            count += 1
    assert count == 680
    print(json.dumps({"status": "passed", "records": count, "search_source_sha256": source_hash}, indent=2))


if __name__ == "__main__":
    main()
