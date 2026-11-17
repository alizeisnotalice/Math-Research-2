#!/usr/bin/env python3
"""Independent semantic verification of the old+resume total-degree-46 ledger."""
from __future__ import annotations

import argparse
import gzip
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
AUDITED_SHA = "a557bf120c6bed03adbd8b2ce153a88bf76e678ed94c6bc8815d14dc33edbe4b"
TOTAL, CAP = 46, 6


def sha(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def load_audited():
    if sha(AUDITED) != AUDITED_SHA:
        raise RuntimeError("audited verifier bytes changed")
    spec = importlib.util.spec_from_file_location("audited_chunk_verifier", AUDITED)
    module = importlib.util.module_from_spec(spec)
    assert spec and spec.loader
    spec.loader.exec_module(module)
    return module


def read_jsonl(path: Path) -> list[dict]:
    return [json.loads(line) for line in path.read_text().splitlines() if line.strip()]


def key(row: dict) -> tuple[int, int, int, tuple[int, ...]]:
    return int(row["m"]), int(row["n"]), int(row["index"]), tuple(row["I"])


def payload_bytes(frozen: bytes, remote: bytes, missing: list[dict], plan_sha: str) -> bytes:
    text = frozen.decode()
    prefix, marker, _ = text.partition("\ndef main():")
    if not marker:
        raise ValueError("frozen source missing main marker")
    source_sha = hashlib.sha256(frozen).hexdigest()
    payload = ("SOURCE_SHA=" + repr(source_sha) + "\n" + prefix + "\n" + remote.decode()
               + "\nrun_selected(jobs=" + json.dumps(missing, separators=(",", ":"))
               + ", plan_sha256=" + repr(plan_sha)
               + ", source_sha256=SOURCE_SHA, workers=15, seconds=14400)\n")
    return payload.encode()


def expected_table(audited) -> dict[tuple[int, int], dict[int, tuple[int, ...]]]:
    return {(m, TOTAL - m): {i: support for i, support in enumerate(audited.independent_supports(m, CAP))}
            for m in range(1, TOTAL // 2 + 1)}


def partition_errors(table, old_chunks: list[dict], new_chunks: list[dict], plan: dict) -> list[dict]:
    errors = []
    plan_received = {key(x): x for x in plan.get("received", [])}
    plan_missing = {key(x): x for x in plan.get("missing", [])}
    old_seen, new_seen = set(), set()
    for label, rows, seen, expected_set in (("old", old_chunks, old_seen, plan_received),
                                             ("new", new_chunks, new_seen, plan_missing)):
        for row in rows:
            try:
                item = key(row); pair = item[:2]
            except Exception as exc:
                errors.append({"kind": "malformed_chunk", "source": label, "error": repr(exc)}); continue
            if pair not in table or item[2] not in table[pair] or table[pair][item[2]] != item[3]:
                errors.append({"kind": "index_or_I_outside_independent_table", "source": label, "job": repr(item)})
                continue
            if item in seen:
                errors.append({"kind": "duplicate_chunk", "source": label, "job": repr(item)}); continue
            seen.add(item)
            if item not in expected_set:
                errors.append({"kind": "chunk_not_in_bound_plan_side", "source": label, "job": repr(item)})
                continue
            bound = expected_set[item]
            for name in ("gzip_sha256", "raw_sha256"):
                # Old received rows are hash-bound by the resume plan.  The
                # plan intentionally contains no precomputed bytes for jobs
                # that did not yet exist; those new hashes are instead bound
                # by the immutable payload/source and independently checked
                # against their local gzip and decompressed ledger below.
                if name in bound and row.get(name) != bound.get(name):
                    errors.append({"kind": "event_hash_disagrees_with_plan", "source": label, "job": repr(item), "field": name})
    if set(plan_received) != old_seen:
        errors.append({"kind": "old_events_not_exact_plan_received", "expected": len(plan_received), "seen": len(old_seen)})
    if set(plan_missing) != new_seen:
        errors.append({"kind": "new_events_not_exact_plan_missing", "expected": len(plan_missing), "seen": len(new_seen)})
    if old_seen & new_seen:
        errors.append({"kind": "old_new_overlap", "count": len(old_seen & new_seen)})
    all_expected = {(m, n, i, support) for (m, n), values in table.items() for i, support in values.items()}
    if old_seen | new_seen != all_expected:
        errors.append({"kind": "combined_partition_not_full_independent_coverage",
                       "expected": len(all_expected), "seen": len(old_seen | new_seen)})
    return errors


def build_manifest(old_root: Path, new_root: Path, plan_path: Path, frozen_path: Path, remote_path: Path) -> dict:
    audited = load_audited(); errors = []
    plan = json.loads(plan_path.read_text())
    old_events = read_jsonl(old_root / "events.jsonl")
    new_events = read_jsonl(new_root / "events.jsonl")
    old_start = old_events[0] if old_events else {}
    if (old_start.get("event"), old_start.get("total_min"), old_start.get("total_max"), old_start.get("cap"), old_start.get("smoke"), old_start.get("workers")) != ("started", 45, 46, 6, False, 15):
        errors.append({"kind": "old_start_contract"})
    old_frozen = old_root.parent / "frozen_support.py"; old_remote = old_root.parent / "remote.py"
    if not old_frozen.is_file() or old_start.get("source_sha256") != sha(old_frozen):
        errors.append({"kind": "old_frozen_binding"})
    old_payload = old_root / "payload.py"
    if not old_payload.is_file() or not old_remote.is_file() or old_payload.read_bytes() != audited.expected_payload(old_frozen.read_text(), old_remote.read_text(), False):
        errors.append({"kind": "old_payload_binding"})
    if plan.get("format") != "a30-resume-total46-v1" or plan.get("parameters") != {"total": 46, "cap": 6, "remote_workers": 15, "remote_seconds_limit": 14400}:
        errors.append({"kind": "plan_contract"})
    original = plan.get("original", {})
    if (original.get("events_sha256") != sha(old_root / "events.jsonl")
            or original.get("frozen_support_sha256") != sha(old_frozen)
            or original.get("remote_sha256") != sha(old_remote)):
        errors.append({"kind": "plan_old_input_binding"})
    plan_sha = sha(plan_path)
    manifest_path, payload_path, receipt_path = new_root / "launch-manifest.json", new_root / "payload.py", new_root / "exit-receipt.json"
    launch = json.loads(manifest_path.read_text()) if manifest_path.is_file() else {}
    frozen, remote = frozen_path.read_bytes(), remote_path.read_bytes()
    if hashlib.sha256(frozen).hexdigest() != original.get("frozen_support_sha256"):
        errors.append({"kind": "resume_frozen_not_bound_to_old_classifier"})
    required_launch = {"format": "a30-resume46-launch-v1", "plan_sha256": plan_sha,
                       "frozen_support_sha256": hashlib.sha256(frozen).hexdigest(),
                       "remote_resume46_sha256": hashlib.sha256(remote).hexdigest(),
                       "selected_jobs": len(plan.get("missing", [])), "workers": 15, "seconds_limit": 14400,
                       "original_events_sha256": sha(old_root / "events.jsonl")}
    if any(launch.get(k) != v for k, v in required_launch.items()):
        errors.append({"kind": "new_launch_manifest_binding"})
    rebuilt = payload_bytes(frozen, remote, plan.get("missing", []), plan_sha)
    if not payload_path.is_file() or payload_path.read_bytes() != rebuilt or launch.get("payload_sha256") != hashlib.sha256(rebuilt).hexdigest():
        errors.append({"kind": "new_payload_binding"})
    receipt = json.loads(receipt_path.read_text()) if receipt_path.is_file() else {}
    if not (receipt.get("returncode") == 0 and receipt.get("selected_set_complete") is True
            and receipt.get("chunks_persisted") == len(plan.get("missing", []))
            and receipt.get("terminal_event") == "resume46_selected_complete"
            and receipt.get("not_parent_45_46_complete") is True):
        errors.append({"kind": "new_exit_receipt"})
    if len(new_events) < 4 or new_events[0].get("event") != "resume46_local_started":
        errors.append({"kind": "new_event_start"})
    elif new_events[0].get("plan_sha256") != plan_sha or new_events[0].get("payload_sha256") != hashlib.sha256(rebuilt).hexdigest():
        errors.append({"kind": "new_local_start_binding"})
    remote_start = next((x for x in new_events if x.get("event") == "resume46_remote_started"), None)
    worker_event = next((x for x in new_events if x.get("event") == "resume46_remote_workers"), None)
    terminal = new_events[-1] if new_events else {}
    if not remote_start or remote_start.get("frozen_source_sha256") != hashlib.sha256(frozen).hexdigest() or remote_start.get("plan_sha256") != plan_sha or remote_start.get("workers") != 15:
        errors.append({"kind": "new_remote_start_binding"})
    if not worker_event or len(worker_event.get("worker_pids", [])) != 15:
        errors.append({"kind": "new_remote_worker_ids"})
    if terminal.get("event") != "resume46_selected_complete" or terminal.get("workers_joined") is not True or terminal.get("submitted") != len(plan.get("missing", [])) or terminal.get("finished") != len(plan.get("missing", [])):
        errors.append({"kind": "new_remote_terminal"})
    old_chunks = [x for x in old_events if x.get("event") == "chunk" and int(x["m"]) + int(x["n"]) == TOTAL]
    new_chunks = [x for x in new_events if x.get("event") == "chunk"]
    table = expected_table(audited)
    errors += partition_errors(table, old_chunks, new_chunks, plan)
    return {"errors": errors, "plan_sha256": plan_sha, "old_events_sha256": sha(old_root / "events.jsonl"),
            "new_events_sha256": sha(new_root / "events.jsonl"), "old_frozen_sha256": sha(old_frozen),
            "old_remote_sha256": sha(old_remote), "new_frozen_sha256": hashlib.sha256(frozen).hexdigest(),
            "new_remote_resume_sha256": hashlib.sha256(remote).hexdigest(), "new_payload_sha256": hashlib.sha256(rebuilt).hexdigest(),
            "new_exit_receipt_sha256": sha(receipt_path), "audited_verifier_sha256": sha(AUDITED),
            "old_chunks": old_chunks, "new_chunks": new_chunks, "table": table,
            "receipt": receipt, "launch": launch}


def cache_key(source: str, event: dict) -> str:
    return f"{source}-m{int(event['m']):02d}n{int(event['n']):02d}-i{int(event['index']):06d}.json"


def worker(root_text: str, event: dict) -> dict:
    audited = load_audited(); errors = []
    answer = audited.verify_chunk(Path(root_text), event, CAP, errors)
    return {"job": answer["job"], "errors": errors, "verified_rows": answer["verified_rows"],
            "counts": dict(answer["counts"]), "event_gzip_sha256": event["gzip_sha256"], "event_raw_sha256": event["raw_sha256"]}


def valid_cache(cache: Path, root: Path, event: dict, binding: dict) -> dict | None:
    try:
        row = json.loads(cache.read_text())
        if row.get("binding") != binding or row.get("event") != {"gzip_sha256": event["gzip_sha256"], "raw_sha256": event["raw_sha256"]}:
            return None
        target = root / f"m{int(event['m']):02d}n{int(event['n']):02d}" / f"i{int(event['index']):06d}.jsonl.gz"
        if not target.is_file() or sha(target) != event["gzip_sha256"]:
            return None
        answer = row.get("result")
        return answer if isinstance(answer, dict) and not answer.get("errors") else None
    except (OSError, ValueError, KeyError, TypeError):
        return None


def atomic(path: Path, value: dict) -> None:
    temp = path.with_suffix(path.suffix + ".tmp")
    temp.write_text(json.dumps(value, sort_keys=True) + "\n")
    temp.replace(path)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--old-run", type=Path, required=True); ap.add_argument("--new-run", type=Path, required=True)
    ap.add_argument("--plan", type=Path, required=True); ap.add_argument("--frozen", type=Path, required=True)
    ap.add_argument("--remote-resume", type=Path, required=True); ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--cache", type=Path, required=True); ap.add_argument("--workers", type=int, default=4)
    args = ap.parse_args()
    if not Path("/Volumes/data").is_mount() or not 1 <= args.workers <= 4:
        raise SystemExit("data mount / worker contract failed")
    began = time.monotonic(); manifest = build_manifest(args.old_run, args.new_run, args.plan, args.frozen, args.remote_resume)
    report = {"format": "a30-combined-total46-independent-verifier-v1", "status": "rejected_manifest", "total": 46,
              "cap": 6, "workers": args.workers, "manifest": {k: v for k, v in manifest.items() if k not in ("old_chunks", "new_chunks", "table")},
              "chunk_results": [], "counts": {}, "cache": {}}
    if not manifest["errors"]:
        args.cache.mkdir(parents=True, exist_ok=True)
        binding = {k: manifest[k] for k in ("plan_sha256", "old_events_sha256", "new_events_sha256", "old_frozen_sha256", "new_frozen_sha256", "new_remote_resume_sha256", "new_payload_sha256", "new_exit_receipt_sha256", "audited_verifier_sha256")}
        jobs = [("old", args.old_run, x) for x in manifest["old_chunks"]] + [("new", args.new_run, x) for x in manifest["new_chunks"]]
        totals = Counter(raw=0, canonical=0, excluded=0, fair=0, unknown=0); pending = []; reused = 0
        for source, root, event in jobs:
            cached = valid_cache(args.cache / cache_key(source, event), root, event, binding)
            if cached is None: pending.append((source, root, event))
            else:
                cached["source"] = source; report["chunk_results"].append(cached); totals.update(cached["counts"]); reused += 1
        with ProcessPoolExecutor(max_workers=args.workers) as pool:
            futures = {pool.submit(worker, str(root), event): (source, root, event) for source, root, event in pending}
            for future in as_completed(futures):
                source, root, event = futures[future]; answer = future.result(); answer["source"] = source
                report["chunk_results"].append(answer); totals.update(answer["counts"])
                atomic(args.cache / cache_key(source, event), {"binding": binding, "event": {"gzip_sha256": event["gzip_sha256"], "raw_sha256": event["raw_sha256"]}, "result": {k: v for k, v in answer.items() if k != "source"}})
        report["chunk_results"].sort(key=lambda x: (x["source"], x["job"])); report["counts"] = dict(totals)
        report["cache"] = {"binding": binding, "reused": reused, "scheduled": len(pending), "directory": str(args.cache)}
        report["status"] = "passed_complete_combined_total46" if not any(x["errors"] for x in report["chunk_results"]) else "rejected_chunk_content"
    report["elapsed_seconds"] = time.monotonic() - began
    report["scope"] = "All total-degree-46 support-index chunks only; old and resumed sources are independently bound and recomputed. Unknown is retained, never excluded."
    args.out.parent.mkdir(parents=True, exist_ok=True); args.out.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"status": report["status"], "chunks": len(report["chunk_results"]), "counts": report["counts"], "elapsed_seconds": report["elapsed_seconds"]}, sort_keys=True))
    if report["status"] != "passed_complete_combined_total46": raise SystemExit(2)


if __name__ == "__main__":
    main()
