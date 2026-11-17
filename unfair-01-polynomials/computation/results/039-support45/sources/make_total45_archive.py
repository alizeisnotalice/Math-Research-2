#!/usr/bin/env python3
"""Make small archival indexes from a completed, hash-bound total-45 audit.

The large per-chunk cache and full semantic report stay in CodexTemp.  This
writer reads every cache row and every current gzip again; it does not trust
the earlier aggregate report merely because it exists.
"""
from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import os
import time
from collections import Counter, defaultdict
from pathlib import Path

from verify_total45 import build_manifest, cache_binding, cache_key, sha


def write_json(path: Path, value) -> str:
    data = json.dumps(value, indent=2, sort_keys=True).encode() + b"\n"
    path.write_bytes(data)
    return hashlib.sha256(data).hexdigest()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("run", type=Path)
    ap.add_argument("cache", type=Path)
    ap.add_argument("report", type=Path)
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--prior-runner-pid", type=int, required=True)
    args = ap.parse_args()
    assert Path("/Volumes/data").is_mount()
    args.out.mkdir(parents=True, exist_ok=True)
    report = json.loads(args.report.read_text())
    assert report["status"] == "passed_complete_total45_slice_only"
    manifest = build_manifest(args.run)
    assert not manifest["errors"] and all(not x["missing"] for x in manifest["coverage"])
    binding = cache_binding(manifest)
    assert json.loads((args.cache / "binding.json").read_text()) == binding
    totals = Counter(raw=0, canonical=0, excluded=0, fair=0, unknown=0)
    pairs = defaultdict(lambda: Counter(raw=0, canonical=0, excluded=0, fair=0, unknown=0))
    rows = []
    cache_merkle = hashlib.sha256()
    for event in sorted(manifest["events"], key=cache_key):
        cache_path = args.cache / (cache_key(event) + ".json")
        cache_bytes = cache_path.read_bytes()
        cache = json.loads(cache_bytes)
        event_binding = {"gzip_sha256": event["gzip_sha256"], "raw_sha256": event["raw_sha256"]}
        assert cache["binding"] == binding and cache["event"] == event_binding
        result = cache["result"]
        assert not result["errors"] and result["event_gzip_sha256"] == event["gzip_sha256"] and result["event_raw_sha256"] == event["raw_sha256"]
        packed = args.run / f"m{event['m']:02d}n{event['n']:02d}" / f"i{event['index']:06d}.jsonl.gz"
        assert sha(packed) == event["gzip_sha256"]
        counts = result["counts"]
        assert counts == event["counts"]
        totals.update(counts); pairs[(event["m"], event["n"])].update(counts)
        cache_merkle.update(cache_path.name.encode() + b"\0" + hashlib.sha256(cache_bytes).digest())
        rows.append({"m": event["m"], "n": event["n"], "index": event["index"], "I": event["I"],
                     "counts": counts, **event_binding, "cache_sha256": hashlib.sha256(cache_bytes).hexdigest()})
    assert len(rows) == 35442
    assert dict(totals) == report["counts"] == {"raw": 606304193, "canonical": 303168708, "excluded": 301344396, "fair": 1824312, "unknown": 0}
    coverage = {(x["m"], x["n"]): x for x in manifest["coverage"]}
    degree_pairs = []
    for pair in sorted(pairs):
        c = coverage[pair]
        assert c["seen_I_count"] == c["expected_I_count"] and not c["missing"]
        degree_pairs.append({"m": pair[0], "n": pair[1], "I_count": c["seen_I_count"], "counts": dict(pairs[pair])})
    index = {"format": "a30-total45-chunk-binding-index-v1", "scope": "only total-degree 45; neither total-degree 46 nor the reboot-interrupted 45/46 parent batch is certified",
             "binding": binding, "chunks": len(rows), "rows": rows, "cache_records_merkle_input_sha256": cache_merkle.hexdigest()}
    compressed = gzip.compress(json.dumps(index, separators=(",", ":")).encode(), compresslevel=9, mtime=0)
    index_path = args.out / "total45-chunk-bindings.json.gz"
    assert len(compressed) < 10 * 1024 * 1024
    index_path.write_bytes(compressed)
    pairs_doc = {"status": "passed_complete_total45_slice_only", "total_degree": 45, "cap_each": 6,
                 "degree_pair_count": len(degree_pairs), "degree_pairs": degree_pairs, "totals": dict(totals),
                 "cache_records_merkle_input_sha256": cache_merkle.hexdigest(), "binding": binding,
                 "scope": index["scope"]}
    pairs_sha = write_json(args.out / "total45-degree-pairs.json", pairs_doc)
    try:
        os.kill(args.prior_runner_pid, 0)
        runner_absent = False
    except ProcessLookupError:
        runner_absent = True
    receipt = {"status": "completed_total45_slice_only", "semantic_audit_report_sha256": sha(args.report),
               "semantic_audit_report_bytes": args.report.stat().st_size, "semantic_audit_status": report["status"],
               "prior_runner_pid": args.prior_runner_pid, "prior_runner_absent_at_receipt": runner_absent,
               "chunk_binding_index": {"path": index_path.name, "bytes": len(compressed), "sha256": hashlib.sha256(compressed).hexdigest(), "chunks": len(rows)},
               "degree_pair_aggregate": {"path": "total45-degree-pairs.json", "sha256": pairs_sha},
               "cache_merkle_input_sha256": cache_merkle.hexdigest(),
               "claim_boundary": index["scope"], "generated_unix": time.time()}
    receipt_sha = write_json(args.out / "total45-exit-receipt.json", receipt)
    print(json.dumps({"index_bytes": len(compressed), "index_sha256": receipt["chunk_binding_index"]["sha256"],
                      "degree_pair_sha256": pairs_sha, "receipt_sha256": receipt_sha, "runner_absent": runner_absent}, sort_keys=True))


if __name__ == "__main__":
    main()
