#!/usr/bin/env python3
"""Verify the complete total-degree-45 slice of an interrupted 45/46 run.

The source run has no terminal event because the host rebooted.  This program
therefore does *not* accept or synthesize a completion receipt for the full
45/46 run.  It validates only the independently enumerable total-45 subset.
"""
from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import os
import time
from collections import Counter, defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

HERE = Path(__file__).resolve().parent
AUDITED = HERE / "audited_chunk_verifier.py"
EXPECTED_AUDITED_SHA256 = "a557bf120c6bed03adbd8b2ce153a88bf76e678ed94c6bc8815d14dc33edbe4b"
TOTAL = 45
CAP = 6


def sha(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def atomic_json(path: Path, value) -> None:
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n")
    temporary.replace(path)


def load_audited():
    if sha(AUDITED) != EXPECTED_AUDITED_SHA256:
        raise RuntimeError("audited verifier bytes changed")
    spec = importlib.util.spec_from_file_location("audited_chunk_verifier", AUDITED)
    assert spec and spec.loader
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def canonical_events(root: Path):
    rows = [json.loads(line) for line in (root / "events.jsonl").read_text().splitlines() if line]
    assert rows and rows[0].get("event") == "started"
    start = rows[0]
    assert (start.get("total_min"), start.get("total_max"), start.get("cap"), start.get("smoke"), start.get("workers")) == (45, 46, 6, False, 15)
    source = root.parent / "frozen_support.py"
    remote = root.parent / "remote.py"
    assert start["source_sha256"] == sha(source)
    # The bridge payload must be exactly reproducible even though its missing
    # exit receipt cannot certify the interrupted parent batch.
    payload = root / "payload.py"
    assert payload.is_file()
    audited = load_audited()
    assert payload.read_bytes() == audited.expected_payload(source.read_text(), remote.read_text(), False)
    selected = []
    for row in rows[1:]:
        if row.get("event") == "chunk" and int(row["m"]) + int(row["n"]) == TOTAL:
            selected.append(row)
    return start, selected


def expected_by_pair(audited):
    return {(m, TOTAL - m): {index: support for index, support in enumerate(audited.independent_supports(m, CAP))}
            for m in range(1, TOTAL // 2 + 1)}


def build_manifest(root: Path):
    audited = load_audited()
    start, rows = canonical_events(root)
    expected = expected_by_pair(audited)
    seen = defaultdict(set)
    errors = []
    jobs = []
    for row in rows:
        pair = (int(row["m"]), int(row["n"]))
        index = int(row["index"])
        support = tuple(row["I"])
        if pair not in expected or index not in expected[pair]:
            errors.append({"kind": "outside_total45_expected_index", "pair": pair, "index": index})
            continue
        if expected[pair][index] != support:
            errors.append({"kind": "wrong_I_for_index", "pair": pair, "index": index})
            continue
        if index in seen[pair]:
            errors.append({"kind": "duplicate_total45_I_index", "pair": pair, "index": index})
            continue
        seen[pair].add(index)
        jobs.append(row)
    coverage = []
    for pair, values in sorted(expected.items()):
        missing = sorted(set(values) - seen[pair])
        coverage.append({"m": pair[0], "n": pair[1], "expected_I_count": len(values), "seen_I_count": len(seen[pair]), "missing": missing})
    if any(x["missing"] for x in coverage):
        errors.append({"kind": "incomplete_total45_coverage"})
    return {"start": start, "events": jobs, "coverage": coverage, "errors": errors,
            "source_sha256": sha(root.parent / "frozen_support.py"), "remote_sha256": sha(root.parent / "remote.py"),
            "events_sha256": sha(root / "events.jsonl"), "audited_verifier_sha256": sha(AUDITED)}


def worker(root_text: str, event: dict):
    root = Path(root_text)
    audited = load_audited()
    errors = []
    result = audited.verify_chunk(root, event, CAP, errors)
    return {"job": result["job"], "errors": errors, "verified_rows": result["verified_rows"],
            "counts": dict(result["counts"]), "event_gzip_sha256": event["gzip_sha256"],
            "event_raw_sha256": event["raw_sha256"]}


def cache_key(event: dict) -> str:
    return f"m{int(event['m']):02d}n{int(event['n']):02d}-i{int(event['index']):06d}"


def cache_binding(manifest: dict) -> dict:
    return {key: manifest[key] for key in ("source_sha256", "remote_sha256", "events_sha256", "audited_verifier_sha256")}


def valid_cached(root: Path, event: dict, cache_file: Path, binding: dict) -> dict | None:
    """Reuse a result only if all source/event/current-gzip bindings agree."""
    try:
        row = json.loads(cache_file.read_text())
        if row.get("binding") != binding or row.get("event") != {
            "gzip_sha256": event["gzip_sha256"], "raw_sha256": event["raw_sha256"]
        }:
            return None
        packed = root / f"m{int(event['m']):02d}n{int(event['n']):02d}" / f"i{int(event['index']):06d}.jsonl.gz"
        if not packed.is_file() or sha(packed) != event["gzip_sha256"]:
            return None
        result = row.get("result")
        if not isinstance(result, dict) or result.get("errors"):
            return None
        return result
    except (OSError, ValueError, TypeError, KeyError):
        return None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("run", type=Path)
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--workers", type=int, default=4)
    ap.add_argument("--manifest-only", action="store_true")
    ap.add_argument("--cache", type=Path, default=HERE / "total45-cache")
    args = ap.parse_args()
    assert Path("/Volumes/data").is_mount()
    assert args.run.is_dir() and 1 <= args.workers <= 4
    started = time.monotonic()
    manifest = build_manifest(args.run)
    report = {"status": "partial_parent_batch_total45_slice_under_verification", "parent_batch": {"total_min": 45, "total_max": 46, "terminal_receipt": "absent due to reboot; no completion inferred"},
              "total_degree": TOTAL, "cap_each": CAP, "manifest": manifest, "started_seconds": started,
              "workers": args.workers, "chunk_results": [], "counts": {}, "cache": {}}
    if manifest["errors"]:
        report["status"] = "rejected_manifest"
    elif args.manifest_only:
        report["status"] = "total45_manifest_complete_but_contents_unchecked"
    else:
        args.cache.mkdir(parents=True, exist_ok=True)
        binding = cache_binding(manifest)
        atomic_json(args.cache / "binding.json", binding)
        totals = Counter(raw=0, canonical=0, excluded=0, fair=0, unknown=0)
        pending = []
        reused = 0
        for event in manifest["events"]:
            cached = valid_cached(args.run, event, args.cache / (cache_key(event) + ".json"), binding)
            if cached is None:
                pending.append(event)
            else:
                report["chunk_results"].append(cached)
                totals.update(cached["counts"])
                reused += 1
        report["cache"] = {"directory": str(args.cache), "binding": binding, "reused_chunks": reused,
                           "scheduled_chunks": len(pending), "current_gzip_sha256_rechecked_before_reuse": True}
        with ProcessPoolExecutor(max_workers=args.workers) as pool:
            futures = {pool.submit(worker, str(args.run), row): row for row in pending}
            for future in as_completed(futures):
                row = future.result()
                report["chunk_results"].append(row)
                totals.update(row["counts"])
                event = futures[future]
                atomic_json(args.cache / (cache_key(event) + ".json"), {
                    "binding": binding,
                    "event": {"gzip_sha256": event["gzip_sha256"], "raw_sha256": event["raw_sha256"]},
                    "result": row,
                })
                if row["errors"]:
                    report["status"] = "rejected_chunk_content"
        report["chunk_results"].sort(key=lambda x: x["job"])
        report["counts"] = dict(totals)
        if report["status"] == "partial_parent_batch_total45_slice_under_verification":
            report["status"] = "passed_complete_total45_slice_only"
    report["elapsed_seconds"] = time.monotonic() - started
    report["scope"] = ("Every total-degree-45 chunk has independently recomputed its canonical supports, singleton-unit witnesses, hashes and counts. "
                       "The total-degree-46 portion and the interrupted parent 45/46 batch remain incomplete and unclaimed.")
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"status": report["status"], "chunks": len(report["chunk_results"]), "counts": report["counts"], "elapsed_seconds": report["elapsed_seconds"]}, sort_keys=True))
    if report["status"] != "passed_complete_total45_slice_only":
        raise SystemExit(2)


if __name__ == "__main__":
    main()
