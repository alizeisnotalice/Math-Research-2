#!/usr/bin/env python3
"""Independent metadata and submitted-certificate audit for recovered round 16.

This is intentionally a small read-only collector.  It does not replay the
multi-billion-record gzip ledgers and therefore never represents itself as a
second full record replay; each shard's existing verifier record is checked
for the corresponding ledger and source hash.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
import sys
from collections import Counter
from datetime import datetime, timezone
from math import comb
from pathlib import Path


ap=argparse.ArgumentParser()
ap.add_argument('--batch-dir',type=Path,required=True)
ap.add_argument('--out',type=Path,required=True)
args=ap.parse_args()
TOPIC=Path(__file__).resolve().parents[1]
OUT=args.batch_dir

SEARCH = TOPIC / "computation" / "support_search_layer.py"
VERIFY = TOPIC / "computation" / "verify_support_layer.py"
ROUND19 = TOPIC / "computation" / "check_round19_candidates.py"
EXCEPTIONS = TOPIC / "computation" / "check_layer19_exceptions.py"
CANDIDATES = TOPIC / "computation" / "results" / "019-candidate-certificates.json"


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def support_key(row: dict) -> tuple[tuple[int, ...], tuple[int, ...]]:
    return tuple(row["I"]), tuple(row["J"])


def number_of_supports(degree: int, cap: int) -> int:
    return sum(comb(degree - 1, j) for j in range(min(degree - 1, cap - 2) + 1))


def main() -> None:
    assert Path("/Volumes/data").is_mount()
    batch = json.loads((OUT / "batch.json").read_text())
    assert batch["all_workers_exited"] is True
    assert len(batch["jobs"]) == 36
    assert all(job["status"] == "verified" and job["process_exited"] is True for job in batch["jobs"])
    search_sha, verifier_sha = sha(SEARCH), sha(VERIFY)
    group_records = []
    unknown_rows: list[dict] = []
    expected_configs = {
        (39, 40, 6): set(range(1, 21)),
        (31, 32, 9): set(range(1, 17)),
    }
    for config_key, expected_m in expected_configs.items():
        low, high, cap = config_key
        tag = "six40" if cap == 6 else "nine32"
        selected = [job for job in batch["jobs"] if job["tag"] == tag]
        assert {job["m"] for job in selected} == expected_m
        totals = Counter()
        degree_pairs = set()
        group_unknown = []
        for m in sorted(expected_m):
            folder = OUT / f"{tag}-m{m:02d}"
            state = json.loads((folder / "job.json").read_text())
            summary = json.loads((folder / "summary.json").read_text())
            verification = json.loads((folder / "verification.json").read_text())
            assert state["status"] == "verified" and state["process_exited"] is True
            assert summary["complete"] is True and verification["verified"] is True
            assert summary["counts"] == verification["counts"]
            assert summary["ledger_sha256_uncompressed"] == verification["ledger_sha256_uncompressed"]
            assert summary["script_sha256"] == search_sha
            assert verification["verifier_sha256"] == verifier_sha
            subtotal = Counter()
            for shard in summary["shards"]:
                n = shard["n"]
                assert shard["m"] == m and m <= n and low <= m + n <= high
                assert (m, n) not in degree_pairs
                degree_pairs.add((m, n))
                assert shard["raw"] == number_of_supports(m, cap) * number_of_supports(n, cap)
                assert shard["canonical"] == shard["excluded"] + shard["fair"] + shard["unknown"]
                subtotal.update({key: shard[key] for key in summary["counts"]})
            assert dict(subtotal) == summary["counts"]
            rows = [json.loads(line) for line in (folder / "unknown.jsonl").read_text().splitlines() if line]
            assert len(rows) == summary["counts"]["unknown"]
            totals.update(subtotal)
            group_unknown.extend(rows)
        expected_pairs = {(m, total - m) for total in range(low, high + 1) for m in range(1, total // 2 + 1)}
        assert degree_pairs == expected_pairs
        group_records.append({"range": [low, high], "max_terms": cap, "degree_pairs": len(degree_pairs),
                              "counts": dict(totals), "unknown_rows": len(group_unknown)})
        unknown_rows.extend(group_unknown)
    unknown_keys = {support_key(row) for row in unknown_rows}
    assert len(unknown_keys) == len(unknown_rows) == 57
    certificate_rows = json.loads(CANDIDATES.read_text())["cases"]
    assert len(certificate_rows) == 57
    certificate_keys = {support_key(row) for row in certificate_rows}
    assert len(certificate_keys) == 57 and certificate_keys == unknown_keys
    direct = [row for row in certificate_rows if row["certificate"] is not None]
    special = [row for row in certificate_rows if row["certificate"] is None]
    assert len(direct) == 52 and len(special) == 5
    candidates_run = subprocess.run([sys.executable, str(ROUND19)], cwd=TOPIC / "computation", text=True,
                                    capture_output=True, check=True)
    exceptions_run = subprocess.run([sys.executable, str(EXCEPTIONS)], cwd=TOPIC / "computation", text=True,
                                    capture_output=True, check=True)
    result = {
        "status": "metadata_and_57_candidate_mapping_verified",
        "scope": "round016 recovered support layers only",
        "audit_sha256": sha(Path(__file__)),
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "source_sha256": {"support_search_layer.py": search_sha, "verify_support_layer.py": verifier_sha,
                            "check_round19_candidates.py": sha(ROUND19), "check_layer19_exceptions.py": sha(EXCEPTIONS)},
        "batch": {"jobs": len(batch["jobs"]), "all_workers_exited": True,
                  "recovery_source_sha256": batch["recovery_source_sha256"]},
        "groups": group_records,
        "unknown_certificate_mapping": {"unknowns": len(unknown_rows), "distinct_unknowns": len(unknown_keys),
                                        "candidate_entries": len(certificate_rows), "direct_certificates": len(direct),
                                        "special_certificates": len(special), "sets_match_exactly": True},
        "certificate_reruns": {"round19_candidates": json.loads(candidates_run.stdout),
                                 "round19_exceptions": json.loads(exceptions_run.stdout)},
        "conclusion": "Both target layers have complete independently verified shard ledgers and their 57 residual supports exactly match rerun strict certificates.",
        "limit": "This collector validates completed verifier metadata and reruns the 57 certificates. It does not itself replay every compressed ledger record; that replay is the existing verify_support_layer.py execution recorded per shard.",
    }
    output = args.out
    assert not output.exists(), "refusing to overwrite audit evidence"
    output.write_text(json.dumps(result, ensure_ascii=False, indent=2) + "\n")
    print(json.dumps({"status": result["status"], "groups": group_records,
                      "mapping": result["unknown_certificate_mapping"]}, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
