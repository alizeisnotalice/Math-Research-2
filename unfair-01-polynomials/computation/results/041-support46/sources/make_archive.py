#!/usr/bin/env python3
"""Create compact, independently inspectable artifacts after a passed total-46 report."""
from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import os
from collections import Counter, defaultdict
from pathlib import Path


def sha(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1 << 20), b""): h.update(block)
    return h.hexdigest()


def write_gzip(path: Path, obj) -> tuple[str, str]:
    raw = (json.dumps(obj, sort_keys=True, separators=(",", ":")) + "\n").encode()
    with gzip.GzipFile(path, "wb", compresslevel=9, mtime=0) as f: f.write(raw)
    return sha(path), hashlib.sha256(raw).hexdigest()


def main() -> None:
    ap = argparse.ArgumentParser(); ap.add_argument("report", type=Path); ap.add_argument("--out", type=Path, required=True)
    a = ap.parse_args()
    if a.out.exists(): raise RuntimeError("archive exists")
    report = json.loads(a.report.read_text())
    if report.get("status") != "passed_complete_combined_total46": raise RuntimeError("report not a complete pass")
    rows = report["chunk_results"]; pairs = defaultdict(Counter); bindings = []
    for row in rows:
        source, job = row["source"], row["job"]
        m, n = int(job[1:3]), int(job[4:6])
        pairs[(m, n)]["chunks"] += 1; pairs[(m, n)].update(row["counts"])
        bindings.append({"source": source, "job": job, "gzip_sha256": row["event_gzip_sha256"], "raw_sha256": row["event_raw_sha256"], "counts": row["counts"], "verified_rows": row["verified_rows"]})
    a.out.mkdir(parents=True)
    full_gz, full_raw = write_gzip(a.out / "combined46-full-report.json.gz", report)
    bind_gz, bind_raw = write_gzip(a.out / "combined46-chunk-bindings.json.gz", bindings)
    degree = [{"m": m, "n": n, **dict(pairs[(m, n)])} for m, n in sorted(pairs)]
    (a.out / "combined46-degree-pairs.json").write_text(json.dumps(degree, indent=2, sort_keys=True) + "\n")
    receipt = {"format": "a30-combined46-archive-receipt-v1", "report_status": report["status"], "chunks": len(rows),
               "report_sha256": sha(a.report), "full_gzip_sha256": full_gz, "full_uncompressed_sha256": full_raw,
               "bindings_gzip_sha256": bind_gz, "bindings_uncompressed_sha256": bind_raw,
               "degree_pairs_sha256": sha(a.out / "combined46-degree-pairs.json"),
               "scope": report["scope"]}
    (a.out / "combined46-exit-receipt.json").write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n")
    print(json.dumps(receipt, sort_keys=True))


if __name__ == "__main__": main()
