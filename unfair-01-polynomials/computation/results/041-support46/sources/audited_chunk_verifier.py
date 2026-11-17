#!/usr/bin/env python3
"""Independent verifier for round-038 streamed support ledgers.

This file intentionally reimplements enumeration, canonical normalization, and
the singleton-unit exclusion certificate.  It does not import the searcher's
``frozen_support`` module or any of its helpers.
"""
from __future__ import annotations

import argparse
import gzip
import hashlib
import itertools
import json
import os
import re
import sys
from collections import Counter, defaultdict
from pathlib import Path
from typing import Iterable


def sha256_bytes(x: bytes) -> str:
    return hashlib.sha256(x).hexdigest()


def independent_supports(degree: int, cap: int) -> Iterable[tuple[int, ...]]:
    """All endpoint-normalized 0/1 supports, in the claimed lexicographic order."""
    for size in range(2, min(degree + 1, cap) + 1):
        for inside in itertools.combinations(range(1, degree), size - 2):
            yield (0, *inside, degree)


def independently_canonical(I: tuple[int, ...], J: tuple[int, ...]) -> bool:
    if (len(I), I) > (len(J), J) and I[-1] == J[-1]:
        return False
    if math_gcd_many((*I, *J)) != 1:
        return False
    pair = (I, J)
    reverse = (tuple(I[-1] - x for x in reversed(I)), tuple(J[-1] - x for x in reversed(J)))
    if I[-1] == J[-1] and (len(reverse[0]), reverse[0]) > (len(reverse[1]), reverse[1]):
        reverse = (reverse[1], reverse[0])
    return pair_key(pair) <= pair_key(reverse)


def math_gcd_many(xs: tuple[int, ...]) -> int:
    ans = 0
    for x in xs:
        ans = __import__("math").gcd(ans, x)
    return ans


def pair_key(pair: tuple[tuple[int, ...], tuple[int, ...]]):
    return ((len(pair[0]), pair[0]), (len(pair[1]), pair[1]))


def independently_classify(I: tuple[int, ...], J: tuple[int, ...]):
    diagonals: dict[int, list[tuple[int, int]]] = {}
    for i in I:
        for j in J:
            diagonals.setdefault(i + j, []).append((i, j))
    unit_i, unit_j = {0, I[-1]}, {0, J[-1]}
    # Collect all singleton-unit certificates before testing any collision;
    # this is the stated certificate semantics, independent of dict order.
    for edges in diagonals.values():
        if len(edges) == 1:
            i, j = edges[0]
            unit_i.add(i)
            unit_j.add(j)
    collision = False
    for diagonal in sorted(diagonals):
        edges = diagonals[diagonal]
        if len(edges) > 1:
            collision = True
            for i, j in edges:
                if i in unit_i and j in unit_j:
                    return "excluded", [diagonal, i, j]
    return ("unknown", None) if collision else ("fair", None)


def record_error(errors: list[dict], kind: str, **details) -> None:
    errors.append({"kind": kind, **details})


def ranges(values: set[int]) -> list[list[int]]:
    """Compact consecutive ranges from an explicit expected index set."""
    out, begin, previous = [], None, None
    for x in sorted(values):
        if begin is None:
            begin = previous = x
        elif x == previous + 1:
            previous = x
        else:
            out.append([begin, previous]); begin = previous = x
    if begin is not None:
        out.append([begin, previous])
    return out


def expected_jobs(start: dict) -> dict[tuple[int, int], list[tuple[int, tuple[int, ...]]]]:
    cap = int(start["cap"])
    if start.get("smoke"):
        # This is an independently stated test contract, not imported from the
        # generator: smoke is exactly two known, twenty-record chunks.
        return {(2, 4): [(0, (0, 2))], (3, 5): [(1, (0, 1, 3))]}
    answer = {}
    for total in range(int(start["total_min"]), int(start["total_max"]) + 1):
        for m in range(1, total // 2 + 1):
            n = total - m
            answer[(m, n)] = list(enumerate(independent_supports(m, cap)))
    return answer


def expected_payload(source_text: str, remote_text: str, smoke: bool) -> bytes:
    """Exact deterministic bridge payload construction, reproduced independently."""
    prefix, marker, _ = source_text.partition("\ndef main():")
    if not marker:
        raise ValueError("frozen support source has no main separator")
    source_hash = sha256_bytes(source_text.encode())
    return (
        "SOURCE_SHA=" + repr(source_hash) + "\n" + prefix + "\n" + remote_text
        + f"\nrun(smoke={smoke!r},workers={2 if smoke else 15})\n"
    ).encode()


def verify_chunk(root: Path, event: dict, cap: int, errors: list[dict]) -> dict:
    m, n, index = int(event["m"]), int(event["n"]), int(event["index"])
    I = tuple(event["I"])
    label = f"m{m:02}n{n:02}/i{index:06d}"
    path = root / f"m{m:02}n{n:02}" / f"i{index:06d}.jsonl.gz"
    result = {"job": label, "path": str(path), "verified_rows": 0, "counts": Counter()}
    if not path.is_file():
        record_error(errors, "missing_chunk_file", job=label, path=str(path))
        return result
    packed = path.read_bytes()
    if sha256_bytes(packed) != event.get("gzip_sha256"):
        record_error(errors, "gzip_hash_mismatch", job=label)
        return result
    try:
        raw = gzip.decompress(packed)
    except Exception as exc:
        record_error(errors, "gzip_decode_failure", job=label, error=repr(exc))
        return result
    if sha256_bytes(raw) != event.get("raw_sha256"):
        record_error(errors, "uncompressed_hash_mismatch", job=label)
        return result
    try:
        lines = raw.decode("utf-8").splitlines()
        rows = [json.loads(line) for line in lines]
    except Exception as exc:
        record_error(errors, "ledger_json_failure", job=label, error=repr(exc))
        return result
    # Expected sequence includes every normalized/canonical J, exactly once.
    expected = []
    raw_count = 0
    expected_unknown = []
    expected_counts = Counter(raw=0, canonical=0, excluded=0, fair=0, unknown=0)
    for J in independent_supports(n, cap):
        raw_count += 1
        expected_counts["raw"] += 1
        if not independently_canonical(I, J):
            continue
        status, witness = independently_classify(I, J)
        expected.append([[*I], [*J], status, witness])
        expected_counts["canonical"] += 1
        expected_counts[status] += 1
        if status == "unknown":
            expected_unknown.append([[*I], [*J]])
    if rows != expected:
        # Locate the first semantic mismatch without suppressing hash evidence.
        at = next((q for q, (x, y) in enumerate(itertools.zip_longest(rows, expected)) if x != y), None)
        record_error(errors, "record_sequence_or_certificate_mismatch", job=label, first_difference=at,
                     observed_rows=len(rows), expected_rows=len(expected))
    claimed_counts = event.get("counts")
    if claimed_counts != dict(expected_counts):
        record_error(errors, "event_counts_mismatch", job=label, observed=claimed_counts, expected=dict(expected_counts))
    if event.get("unknown") != expected_unknown:
        record_error(errors, "unknown_listing_mismatch", job=label)
    result["verified_rows"] = len(rows)
    result["counts"] = expected_counts
    result["raw_support_count"] = raw_count
    return result


def verify(root: Path, *, manifest_only: bool = False) -> dict:
    errors: list[dict] = []
    events_file = root / "events.jsonl"
    report = {"format": "round038-independent-support-ledger-verifier-v1", "root": str(root),
              "manifest_only": manifest_only, "errors": errors, "verified_chunks": 0, "verified_rows": 0}
    if not events_file.is_file():
        record_error(errors, "missing_events_file", path=str(events_file))
        report["status"] = "rejected"
        return report
    try:
        events = [json.loads(x) for x in events_file.read_text().splitlines() if x.strip()]
    except Exception as exc:
        record_error(errors, "events_json_failure", error=repr(exc)); report["status"] = "rejected"; return report
    starts = [x for x in events if x.get("event") == "started"]
    if len(starts) != 1:
        record_error(errors, "started_event_count", count=len(starts)); report["status"] = "rejected"; return report
    start = starts[0]
    if not events or events[0] is not start:
        record_error(errors, "started_event_not_first_or_not_unique")
    required = ("source_sha256", "total_min", "total_max", "cap", "smoke")
    for key in required:
        if key not in start:
            record_error(errors, "started_field_missing", field=key)
    if start.get("total_min") != 45 or start.get("total_max") != 46 or start.get("cap") != 6 or not isinstance(start.get("smoke"), bool):
        record_error(errors, "fixed_production_scope_contract_failure",
                     observed={k: start.get(k) for k in ("total_min", "total_max", "cap", "smoke")})
    if start.get("workers") != (2 if start.get("smoke") else 15):
        record_error(errors, "started_workers_contract_failure", observed=start.get("workers"),
                     expected=2 if start.get("smoke") else 15)
    source = root.parent / "frozen_support.py"
    remote = root.parent / "remote.py"
    source_text = source.read_text() if source.is_file() else ""
    if not source.is_file() or start.get("source_sha256") != sha256_bytes(source_text.encode()):
        record_error(errors, "frozen_source_hash_binding_failure", source=str(source))
    payload, exit_file = root / "payload.py", root / "exit.json"
    if payload.is_file():
        payload_hash = sha256_bytes(payload.read_bytes())
        if not remote.is_file():
            record_error(errors, "missing_remote_source", path=str(remote))
        else:
            try:
                rebuilt = expected_payload(source_text, remote.read_text(), bool(start.get("smoke")))
                if payload.read_bytes() != rebuilt:
                    record_error(errors, "reconstructed_payload_binding_failure")
            except Exception as exc:
                record_error(errors, "payload_reconstruction_failure", error=repr(exc))
    else:
        payload_hash = None
        record_error(errors, "missing_payload")
    exit_state = None
    if exit_file.is_file():
        try:
            exit_state = json.loads(exit_file.read_text())
            if payload_hash and exit_state.get("payload_sha256") != payload_hash:
                record_error(errors, "payload_hash_binding_failure")
        except Exception as exc:
            record_error(errors, "exit_json_failure", error=repr(exc))
    jobs = expected_jobs(start)
    expected_by_pair = {pair: {i: I for i, I in values} for pair, values in jobs.items()}
    chunks = [x for x in events if x.get("event") == "chunk"]
    seen: dict[tuple[int, int], set[int]] = defaultdict(set)
    aggregate = Counter(raw=0, canonical=0, excluded=0, fair=0, unknown=0)
    for event in chunks:
        try:
            m, n, index, I = int(event["m"]), int(event["n"]), int(event["index"]), tuple(event["I"])
        except Exception as exc:
            record_error(errors, "malformed_chunk_event", event=event, error=repr(exc)); continue
        pair = (m, n)
        expected_i = expected_by_pair.get(pair, {}).get(index)
        if expected_i is None:
            record_error(errors, "chunk_index_outside_expected_coverage", m=m, n=n, index=index)
            continue
        if I != expected_i:
            record_error(errors, "chunk_I_does_not_match_support_index", m=m, n=n, index=index, observed=list(I), expected=list(expected_i))
            continue
        if index in seen[pair]:
            record_error(errors, "duplicate_chunk_index", m=m, n=n, index=index); continue
        seen[pair].add(index)
        if not manifest_only:
            checked = verify_chunk(root, event, int(start["cap"]), errors)
            report["verified_chunks"] += 1
            report["verified_rows"] += checked["verified_rows"]
            aggregate.update(checked["counts"])
    coverage = []
    all_covered = True
    for (m, n), expected in sorted(expected_by_pair.items()):
        missing = ranges(set(expected) - seen[(m, n)])
        all_covered &= not missing
        coverage.append({"m": m, "n": n, "expected_I_count": len(expected), "seen_I_count": len(seen[(m, n)]), "missing_index_ranges": missing})
    report["coverage"] = coverage
    report["recomputed_counts_from_verified_chunks"] = dict(aggregate)
    terminals = [x for x in events if x.get("event") in ("complete", "partial_timeout")]
    if len(terminals) > 1:
        record_error(errors, "terminal_event_not_unique", count=len(terminals))
    if terminals and events[-1] is not terminals[0]:
        record_error(errors, "terminal_event_not_last")
    terminal = terminals[-1] if terminals else None
    report["terminal_event"] = terminal
    claimed_complete = bool(terminal and terminal.get("event") == "complete" and terminal.get("complete") is True)
    if terminal and terminal.get("workers_joined") is not True:
        record_error(errors, "terminal_workers_not_joined")
    if claimed_complete and not all_covered:
        record_error(errors, "complete_claim_without_full_degreepair_and_I_coverage")
    if claimed_complete and (not exit_state or exit_state.get("returncode") != 0):
        record_error(errors, "complete_claim_without_successful_bridge_exit")
    if claimed_complete and terminal and (terminal.get("submitted") != sum(len(v) for v in jobs.values()) or terminal.get("finished") != sum(len(v) for v in jobs.values())):
        record_error(errors, "complete_claim_submitted_finished_mismatch")
    if errors:
        report["status"] = "rejected"
    elif claimed_complete and all_covered and manifest_only:
        report["status"] = "passed_manifest_only"
    elif claimed_complete and all_covered:
        report["status"] = "passed_complete"
    else:
        report["status"] = "partial_not_complete"
    report["claimed_complete"] = claimed_complete
    report["all_expected_coverage_present"] = all_covered
    return report


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("ledger_root", type=Path)
    ap.add_argument("--report", type=Path, required=True)
    ap.add_argument("--manifest-only", action="store_true")
    ap.add_argument("--require-complete", action="store_true")
    a = ap.parse_args()
    report = verify(a.ledger_root, manifest_only=a.manifest_only)
    a.report.parent.mkdir(parents=True, exist_ok=True)
    a.report.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"status": report["status"], "verified_chunks": report["verified_chunks"], "verified_rows": report["verified_rows"], "errors": len(report["errors"])}, indent=2))
    if report["status"] == "rejected" or (a.require_complete and report["status"] != "passed_complete"):
        raise SystemExit(2)


if __name__ == "__main__":
    main()
