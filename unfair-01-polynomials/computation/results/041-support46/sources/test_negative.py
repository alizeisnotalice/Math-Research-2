#!/usr/bin/env python3
"""Negative tests: partition deletion/duplication and hash-rebound witness change."""
from __future__ import annotations

import copy
import gzip
import hashlib
import importlib.util
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent


def load(name):
    spec = importlib.util.spec_from_file_location(name, HERE / f"{name}.py")
    module = importlib.util.module_from_spec(spec); spec.loader.exec_module(module); return module


def main() -> None:
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("old", type=Path); ap.add_argument("new", type=Path); ap.add_argument("plan", type=Path)
    ap.add_argument("frozen", type=Path); ap.add_argument("remote", type=Path); ap.add_argument("out", type=Path)
    a = ap.parse_args()
    if a.out.exists(): raise RuntimeError("negative-test output exists")
    v = load("verify_combined46"); audited = load("audited_chunk_verifier")
    manifest = v.build_manifest(a.old, a.new, a.plan, a.frozen, a.remote)
    if manifest["errors"]: raise RuntimeError(manifest["errors"][:2])
    plan = json.loads(a.plan.read_text())
    deleted = v.partition_errors(manifest["table"], manifest["old_chunks"][:-1], manifest["new_chunks"], plan)
    repeated = v.partition_errors(manifest["table"], manifest["old_chunks"] + [manifest["old_chunks"][0]], manifest["new_chunks"], plan)
    if not any(x["kind"] == "old_events_not_exact_plan_received" for x in deleted): raise RuntimeError("deletion accepted")
    if not any(x["kind"] == "duplicate_chunk" for x in repeated): raise RuntimeError("duplicate accepted")
    event = next(x for x in manifest["old_chunks"] if x["counts"]["excluded"] > 0)
    source = a.old / f"m{event['m']:02d}n{event['n']:02d}" / f"i{event['index']:06d}.jsonl.gz"
    rows = [json.loads(x) for x in gzip.decompress(source.read_bytes()).decode().splitlines()]
    row = next(x for x in rows if x[2] == "excluded")
    row[3] = [999, 999, 999]  # semantic witness change, then rebind both hashes
    raw = ("\n".join(json.dumps(x, separators=(",", ":")) for x in rows) + "\n").encode()
    packed = gzip.compress(raw, compresslevel=1, mtime=0)
    changed = copy.deepcopy(event); changed["raw_sha256"] = hashlib.sha256(raw).hexdigest(); changed["gzip_sha256"] = hashlib.sha256(packed).hexdigest()
    target = a.out / f"m{event['m']:02d}n{event['n']:02d}"; target.mkdir(parents=True)
    (target / f"i{event['index']:06d}.jsonl.gz").write_bytes(packed)
    errors = []; audited.verify_chunk(a.out, changed, 6, errors)
    if not any(x["kind"] == "record_sequence_or_certificate_mismatch" for x in errors): raise RuntimeError("hash-rebound witness change accepted")
    result = {"format": "a30-combined46-negative-v1", "deleted_chunk_rejected": True, "duplicate_chunk_rejected": True,
              "hash_rebound_witness_change_rejected": True, "source_job": v.key(event), "semantic_errors": errors}
    (a.out / "result.json").write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps({k: v for k, v in result.items() if k != "semantic_errors"}, default=str, sort_keys=True))


if __name__ == "__main__": main()
