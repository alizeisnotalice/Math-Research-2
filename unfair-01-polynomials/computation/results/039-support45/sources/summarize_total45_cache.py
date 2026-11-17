#!/usr/bin/env python3
"""Create a compact, independently hash-bound summary of total45 cache rows."""
from __future__ import annotations

import argparse
import hashlib
import json
import sys
from collections import Counter
from pathlib import Path

from verify_total45 import cache_binding, cache_key, build_manifest, sha


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("run", type=Path)
    ap.add_argument("cache", type=Path)
    ap.add_argument("--out", type=Path, required=True)
    args = ap.parse_args()
    manifest = build_manifest(args.run)
    if manifest["errors"]:
        raise SystemExit("manifest rejected")
    binding = cache_binding(manifest)
    if json.loads((args.cache / "binding.json").read_text()) != binding:
        raise SystemExit("cache binding is stale")
    totals = Counter(raw=0, canonical=0, excluded=0, fair=0, unknown=0)
    h = hashlib.sha256(); names = []
    for event in sorted(manifest["events"], key=cache_key):
        path = args.cache / (cache_key(event) + ".json")
        row = json.loads(path.read_text())
        if row.get("binding") != binding or row.get("event") != {"gzip_sha256": event["gzip_sha256"], "raw_sha256": event["raw_sha256"]}:
            raise SystemExit(f"bad cache binding: {path.name}")
        # Current gzip contents are rehashed, so the summary cannot cite an
        # event digest for a changed on-disk ledger.
        packed = args.run / f"m{event['m']:02d}n{event['n']:02d}" / f"i{event['index']:06d}.jsonl.gz"
        if sha(packed) != event["gzip_sha256"]:
            raise SystemExit(f"current gzip mismatch: {path.name}")
        result = row["result"]
        if result["errors"]:
            raise SystemExit(f"cached semantic failure: {path.name}")
        totals.update(result["counts"])
        b = path.read_bytes(); h.update(path.name.encode() + b"\0" + hashlib.sha256(b).digest()); names.append(path.name)
    coverage = [{k: x[k] for k in ("m", "n", "expected_I_count", "seen_I_count", "missing")} for x in manifest["coverage"]]
    out = {"status": "passed_complete_total45_slice_only", "total_degree": 45, "cap_each": 6,
           "chunks": len(names), "counts": dict(totals), "cache_binding": binding,
           "cache_records_merkle_input_sha256": h.hexdigest(), "coverage": coverage,
           "scope": "Each chunk is separately recorded in the uncommitted cache. This compact report does not certify any degree-46 record or the interrupted parent 45/46 run.",
           "interpreter": sys.version}
    args.out.write_text(json.dumps(out, indent=2, sort_keys=True) + "\n")
    print(json.dumps({k: out[k] for k in ("status", "chunks", "counts", "cache_records_merkle_input_sha256")}, sort_keys=True))


if __name__ == "__main__":
    main()
