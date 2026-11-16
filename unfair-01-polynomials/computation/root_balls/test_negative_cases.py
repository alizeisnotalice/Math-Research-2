#!/usr/bin/env python3
"""Negative tests for the hardened R8--L9 rational root-ball verifier.

Each case starts from the archived certificate/input pair, makes only the
stated adversarial mutation in a data-volume temporary directory, and requires
the verifier to reject it.  In particular, the wrong-target case updates both
the candidate C hash and the certificate input-file hash, so rejection cannot
be attributed merely to a stale hash.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import subprocess
import sys
import tempfile
from copy import deepcopy
from pathlib import Path


HERE = Path(__file__).resolve().parent
VERIFIER = HERE / "verify_root_balls.py"
ap = argparse.ArgumentParser()
ap.add_argument('--certificate', type=Path, required=True)
ap.add_argument('--input', type=Path, required=True)
ap.add_argument('--record', type=Path, required=True)
ap.add_argument('--work', type=Path, required=True)
args = ap.parse_args()
CERTIFICATE, INPUT, REPORT = args.certificate, args.input, args.record



def write_json(path: Path, value: object) -> None:
    path.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n")


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def coefficient_sha256(coefficients: list[int]) -> str:
    return hashlib.sha256(json.dumps(coefficients, separators=(",", ":")).encode()).hexdigest()


def verify_rejects(label: str, certificate: object, source: object, expected: str, work: Path) -> dict[str, object]:
    certificate_path = work / f"{label}-certificate.json"
    input_path = work / f"{label}-input.json"
    record_path = work / f"{label}-record.json"
    write_json(input_path, source)
    # Every fixture is serialized anew, so keep the outer file binding valid.
    # This makes every rejection about the intended adversarial field rather
    # than about a stale source-file digest.  The wrong-target fixture also
    # independently refreshes its internal C_sha256 above.
    certificate["input_sha256"] = sha(input_path)
    write_json(certificate_path, certificate)
    env = dict(os.environ)
    env["PYTHONDONTWRITEBYTECODE"] = "1"
    env["OMP_NUM_THREADS"] = "1"
    result = subprocess.run(
        [sys.executable, str(VERIFIER), "--certificate", str(certificate_path), "--input", str(input_path), "--record", str(record_path)],
        cwd=HERE, text=True, capture_output=True, env=env, check=False,
    )
    combined = result.stdout + result.stderr
    if result.returncode == 0:
        raise AssertionError(f"{label}: mutated input was accepted")
    if expected not in combined:
        raise AssertionError(f"{label}: expected rejection {expected!r}, got {combined!r}")
    return {"case": label, "rejected": True, "returncode": result.returncode, "expected_message": expected}


def main() -> None:
    if not Path("/Volumes/data").is_mount():
        raise SystemExit("/Volumes/data is not mounted; refusing to create test files elsewhere")
    original_certificate = json.loads(CERTIFICATE.read_text())
    original_source = json.loads(INPUT.read_text())
    cases = []
    with tempfile.TemporaryDirectory(prefix="root-ball-negative-", dir=args.work) as temporary:
        work = Path(temporary)
        certificate = deepcopy(original_certificate)
        certificate["groups"][0]["representative"]["radius"] = "-1/1"
        cases.append(verify_rejects("negative_radius", certificate, deepcopy(original_source), "nonpositive disk radius", work))

        certificate = deepcopy(original_certificate)
        certificate["groups"][0]["representative"]["radius"] = "0/1"
        cases.append(verify_rejects("zero_radius", certificate, deepcopy(original_source), "nonpositive disk radius", work))

        certificate = deepcopy(original_certificate)
        source = deepcopy(original_source)
        candidate = source["candidates"][0]
        candidate["C_coefficients_ascending"][0] = 0
        candidate["C_sha256"] = coefficient_sha256(candidate["C_coefficients_ascending"])
        cases.append(verify_rejects(
            "wrong_target_with_rehashed_input", certificate, source,
            "source candidate C is not the exact R8(x)R8(x^9) target", work,
        ))

        certificate = deepcopy(original_certificate)
        certificate["groups"][0]["representative"]["center"]["re"] = "0/1"
        cases.append(verify_rejects(
            "changed_center", certificate, deepcopy(original_source),
            "mother center was not copied exactly from source suggestion", work,
        ))

        certificate = deepcopy(original_certificate)
        certificate["groups"][1]["source_group_index"] = certificate["groups"][0]["source_group_index"]
        cases.append(verify_rejects(
            "duplicate_source_coverage", certificate, deepcopy(original_source),
            "duplicate source coverage in certificate groups", work,
        ))
    report = {
        "status": "all_negative_cases_rejected", "certificate": str(CERTIFICATE), "certificate_sha256": sha(CERTIFICATE),
        "input": str(INPUT), "input_sha256": sha(INPUT), "cases": cases,
        "note": "The wrong-target case updates candidate C_sha256 and certificate input_sha256 before invoking the verifier.",
    }
    write_json(REPORT, report)
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
