#!/usr/bin/env python3
"""Smoke and adversarial tests for the total-45 slice verifier."""
from __future__ import annotations

import copy
import gzip
import hashlib
import json
import shutil
import subprocess
import sys
from pathlib import Path

from verify_total45 import build_manifest, cache_binding, cache_key, load_audited, valid_cached

HERE = Path(__file__).resolve().parent
SOURCE = Path("/Volumes/data/CodexTemp/T2/unfair-01-polynomials/round038-macbook-support46")
RUN = SOURCE / "run"
TMP = Path("/Volumes/data/CodexTemp/T2/unfair-01-polynomials/agents/verify_slice45_39/test-fixtures")


def h(b: bytes) -> str:
    return hashlib.sha256(b).hexdigest()


def fixture(name: str) -> Path:
    target = TMP / name
    if target.exists():
        shutil.rmtree(target)
    (target / "run").mkdir(parents=True)
    for name0 in ("frozen_support.py", "remote.py"):
        shutil.copy2(SOURCE / name0, target / name0)
    for name0 in ("events.jsonl", "payload.py"):
        shutil.copy2(RUN / name0, target / "run" / name0)
    return target


def event_rows(path: Path):
    return [json.loads(x) for x in path.read_text().splitlines() if x]


def write_events(path: Path, rows):
    path.write_text("".join(json.dumps(x, separators=(",", ":")) + "\n" for x in rows))


def main():
    # Existing independent 20-record endpoint is executed unmodified.
    smoke = Path("/Volumes/data/CodexTemp/T2/unfair-01-polynomials/round038-macbook-support46/smoke")
    smoke_report = HERE / "smoke-20-report.json"
    command = [sys.executable, str(HERE / "audited_chunk_verifier.py"), str(smoke), "--report", str(smoke_report), "--require-complete"]
    subprocess.run(command, check=True)
    assert json.loads(smoke_report.read_text())["verified_rows"] == 20

    good = fixture("good")
    manifest = build_manifest(good / "run")
    assert not manifest["errors"] and all(not row["missing"] for row in manifest["coverage"])

    deleted = fixture("deleted")
    rows = event_rows(deleted / "run" / "events.jsonl")
    drop = next(i for i, row in enumerate(rows) if row.get("event") == "chunk" and row["m"] + row["n"] == 45)
    rows.pop(drop); write_events(deleted / "run" / "events.jsonl", rows)
    assert any(e["kind"] == "incomplete_total45_coverage" for e in build_manifest(deleted / "run")["errors"])

    repeated = fixture("repeated")
    rows = event_rows(repeated / "run" / "events.jsonl")
    row = next(row for row in rows if row.get("event") == "chunk" and row["m"] + row["n"] == 45)
    rows.append(copy.deepcopy(row)); write_events(repeated / "run" / "events.jsonl", rows)
    assert any(e["kind"] == "duplicate_total45_I_index" for e in build_manifest(repeated / "run")["errors"])

    # Rebind both hashes after changing a witness; semantic recomputation must
    # still reject it, so this is stronger than a digest-only mutation test.
    audited = load_audited()
    rows = event_rows(RUN / "events.jsonl")
    event = next(row for row in rows if row.get("event") == "chunk" and row["m"] + row["n"] == 45 and row["counts"]["excluded"])
    packed = (RUN / f"m{event['m']:02d}n{event['n']:02d}" / f"i{event['index']:06d}.jsonl.gz").read_bytes()
    records = [json.loads(x) for x in gzip.decompress(packed).splitlines()]
    position = next(i for i, record in enumerate(records) if record[2] == "excluded")
    records[position][3] = [0, 0, 0]
    raw = b"".join(json.dumps(x, separators=(",", ":")).encode() + b"\n" for x in records)
    changed = copy.deepcopy(event); changed["raw_sha256"] = h(raw); changed["gzip_sha256"] = h(gzip.compress(raw, mtime=0))
    witness_errors = []
    # Place changed bytes into a dedicated root at the same expected path.
    witness_root = fixture("witness") / "run"
    dest = witness_root / f"m{changed['m']:02d}n{changed['n']:02d}" / f"i{changed['index']:06d}.jsonl.gz"
    dest.parent.mkdir(); dest.write_bytes(gzip.compress(raw, mtime=0))
    audited.verify_chunk(witness_root, changed, 6, witness_errors)
    assert any(e["kind"] == "record_sequence_or_certificate_mismatch" for e in witness_errors)

    # A cached success becomes unusable as soon as the current gzip bytes no
    # longer match the event binding, even if the cache JSON itself is intact.
    cache_file = TMP / "cache.json"
    cache_row = {"binding": cache_binding(manifest), "event": {"gzip_sha256": event["gzip_sha256"], "raw_sha256": event["raw_sha256"]},
                 "result": {"errors": [], "counts": {}, "job": "synthetic"}}
    cache_file.write_text(json.dumps(cache_row))
    assert valid_cached(RUN, event, cache_file, cache_binding(manifest)) is not None
    cache_row["event"]["gzip_sha256"] = "0" * 64
    cache_file.write_text(json.dumps(cache_row))
    assert valid_cached(RUN, event, cache_file, cache_binding(manifest)) is None
    result = {"smoke_rows": 20, "deleted_I": "rejected", "repeated_I": "rejected", "rehash_bound_witness_change": "rejected", "cache_wrong_gzip_binding": "rejected"}
    (HERE / "test-result.json").write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, sort_keys=True))


if __name__ == "__main__":
    main()
