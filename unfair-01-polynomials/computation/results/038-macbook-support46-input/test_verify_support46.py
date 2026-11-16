#!/usr/bin/env python3
"""Positive and adversarial tests for the independent round-038 verifier."""
from __future__ import annotations

import gzip
import hashlib
import json
import shutil
from pathlib import Path

from verify_support46 import verify


HERE = Path(__file__).resolve().parent
SOURCE = HERE
FIXTURES = Path("/Volumes/data/CodexTemp/T2/unfair-01-polynomials/round038-parent/support46-regressions")


def h(x: bytes) -> str:
    return hashlib.sha256(x).hexdigest()


def install(name: str) -> Path:
    target = FIXTURES / name
    if target.exists():
        shutil.rmtree(target)
    target.parent.mkdir(parents=True, exist_ok=True)
    shutil.copytree(SOURCE / "smoke", target)
    shutil.copy2(SOURCE / "frozen_support.py", FIXTURES / "frozen_support.py")
    shutil.copy2(SOURCE / "remote.py", FIXTURES / "remote.py")
    return target


def events(root: Path) -> list[dict]:
    return [json.loads(x) for x in (root / "events.jsonl").read_text().splitlines() if x]


def write_events(root: Path, rows: list[dict]) -> None:
    (root / "events.jsonl").write_text("".join(json.dumps(x, separators=(",", ":")) + "\n" for x in rows))


def replace_gzip_and_rebind(root: Path, event: dict, rows: list) -> None:
    path = root / f"m{event['m']:02d}n{event['n']:02d}" / f"i{event['index']:06d}.jsonl.gz"
    raw = b"".join(json.dumps(x, separators=(",", ":")).encode() + b"\n" for x in rows)
    packed = gzip.compress(raw, compresslevel=1, mtime=0)
    path.write_bytes(packed)
    event["raw_sha256"], event["gzip_sha256"] = h(raw), h(packed)


def require_rejected(report: dict, contains: str) -> None:
    assert report["status"] == "rejected", report
    assert any(contains in x["kind"] for x in report["errors"]), report["errors"]


def main() -> None:
    results = {}
    legal = install("legal")
    positive = verify(legal)
    assert positive["status"] == "passed_complete" and positive["verified_rows"] == 20, positive
    results["legal_20_rows"] = {"status": positive["status"], "rows": positive["verified_rows"]}
    manifest_positive = verify(legal, manifest_only=True)
    assert manifest_positive["status"] == "passed_manifest_only", manifest_positive
    results["manifest_only_never_content_complete"] = manifest_positive["status"]

    witness = install("bad-witness")
    e = events(witness); chunk = next(x for x in e if x["event"] == "chunk")
    path = witness / f"m{chunk['m']:02d}n{chunk['n']:02d}" / f"i{chunk['index']:06d}.jsonl.gz"
    rows = [json.loads(x) for x in gzip.decompress(path.read_bytes()).decode().splitlines()]
    rows[1][3] = [999, 999, 999]
    replace_gzip_and_rebind(witness, chunk, rows); write_events(witness, e)
    bad_witness = verify(witness); require_rejected(bad_witness, "record_sequence_or_certificate_mismatch")
    results["changed_witness"] = bad_witness["status"]

    missing = install("missing-line")
    e = events(missing); chunk = next(x for x in e if x["event"] == "chunk")
    path = missing / f"m{chunk['m']:02d}n{chunk['n']:02d}" / f"i{chunk['index']:06d}.jsonl.gz"
    rows = [json.loads(x) for x in gzip.decompress(path.read_bytes()).decode().splitlines()]
    replace_gzip_and_rebind(missing, chunk, rows[:-1]); write_events(missing, e)
    bad_missing = verify(missing); require_rejected(bad_missing, "record_sequence_or_certificate_mismatch")
    results["missing_line"] = bad_missing["status"]

    wrong = install("wrong-I-index")
    e = events(wrong); chunk = [x for x in e if x["event"] == "chunk"][1]
    chunk["I"] = [0, 2, 3]  # valid-looking normalized support, wrong for index 1
    write_events(wrong, e)
    bad_index = verify(wrong); require_rejected(bad_index, "chunk_I_does_not_match_support_index")
    results["wrong_I_index"] = bad_index["status"]

    wrong_index = install("wrong-index")
    e = events(wrong_index); chunk = [x for x in e if x["event"] == "chunk"][1]
    chunk["index"] = 999999
    write_events(wrong_index, e)
    bad_index_number = verify(wrong_index); require_rejected(bad_index_number, "chunk_index_outside_expected_coverage")
    results["wrong_index_number"] = bad_index_number["status"]

    forged = install("forged-complete")
    e = events(forged); started = next(x for x in e if x["event"] == "started")
    started["smoke"] = False  # falsely recast two test chunks as the full 45/46 production run
    write_events(forged, e)
    bad_complete = verify(forged); require_rejected(bad_complete, "complete_claim_without_full_degreepair_and_I_coverage")
    results["forged_complete"] = bad_complete["status"]

    narrowed = install("narrowed-scope")
    e = events(narrowed); started = next(x for x in e if x["event"] == "started")
    started["total_min"] = 46
    write_events(narrowed, e)
    bad_scope = verify(narrowed); require_rejected(bad_scope, "fixed_production_scope_contract_failure")
    results["narrowed_production_scope"] = bad_scope["status"]

    modified_payload = install("modified-payload")
    payload = modified_payload / "payload.py"
    payload.write_text(payload.read_text().replace("SOURCE_SHA", "SOURCE_SHB", 1))
    exit_path = modified_payload / "exit.json"
    exit_state = json.loads(exit_path.read_text())
    exit_state["payload_sha256"] = h(payload.read_bytes())
    exit_path.write_text(json.dumps(exit_state, indent=2) + "\n")
    bad_payload = verify(modified_payload); require_rejected(bad_payload, "reconstructed_payload_binding_failure")
    results["modified_payload_rehashed_exit"] = bad_payload["status"]

    (HERE / "test-results.json").write_text(json.dumps(results, indent=2) + "\n")
    print(json.dumps(results, indent=2))


if __name__ == "__main__":
    main()
