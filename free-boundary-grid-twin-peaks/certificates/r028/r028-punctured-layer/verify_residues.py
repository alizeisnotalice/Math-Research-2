#!/usr/bin/env python3
"""Reconstruct punctured-layer residues and compare finite exact tables.

This is intentionally separate from the C++ recurrence.  It uses Python
integers for CRT and reads only JSON checkpoints produced one-prime-at-a-time.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

PRIMES = [
    9223372036854775783,
    9223372036854775643,
    9223372036854775549,
    9223372036854775507,
]


def crt(residues: list[int]) -> int:
    modulus = 1
    for p in PRIMES:
        modulus *= p
    value = 0
    for a, p in zip(residues, PRIMES):
        other = modulus // p
        value = (value + a * other * pow(other, -1, p)) % modulus
    return value


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read_group(directory: Path, n: int, v: int) -> tuple[int, int, dict]:
    rows = []
    for pi, p in enumerate(PRIMES):
        path = directory / f"n{n}-v{v}-p{pi}.json"
        row = json.loads(path.read_text())
        assert row["n"] == n and row["v"] == v and row["prime_index"] == pi
        assert int(row["prime"]) == p
        rows.append(row)
    assert len({r["algorithm"] for r in rows}) == 1
    assert len({r["independent_sets"] for r in rows}) == 1
    a = crt([int(r["A2v_residue"]) for r in rows])
    b = crt([int(r["Bv_residue"]) for r in rows])
    return a, b, rows[0]


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--runs", type=Path, required=True)
    ap.add_argument("--table-root", type=Path, required=True)
    ap.add_argument("--table-n6", type=Path, required=True)
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()

    refs: dict[int, dict] = {}
    source_hashes = {}
    for n in (3, 4, 5):
        path = args.table_root / f"punctured-table-n{n}.json"
        refs[n] = json.loads(path.read_text())
        source_hashes[path.name] = sha256(path)
    n6_path = args.table_n6
    refs[6] = json.loads(n6_path.read_text())
    source_hashes[n6_path.name] = sha256(n6_path)

    expected_vs = {3: [0, 1, 4], 4: [0, 1, 5], 5: [0, 1, 2, 6, 7, 12], 6: [0, 1, 2, 7, 8, 14]}
    out = {"status": "passed", "method": "four-prime CRT of modular rolling IE checkpoints", "checks": [],
           "source_table_sha256": source_hashes,
           "modulus_bits": __import__('math').prod(PRIMES).bit_length()}
    for n, vertices in expected_vs.items():
        for v in vertices:
            a, b, meta = read_group(args.runs, n, v)
            record = {"n": n, "v": v, "A2_v": str(a), "B_v": str(b),
                      "independent_sets": meta["independent_sets"],
                      "max_modeled_payload_bytes": meta["max_modeled_payload_bytes"]}
            ref = refs[n]["per_vertex"][v]
            record["table_crosscheck"] = (a == int(ref["A2_v"]) and b == int(ref["B_v"]))
            if not record["table_crosscheck"]:
                raise AssertionError((n, v, record, ref))
            out["checks"].append(record)
    args.output.write_text(json.dumps(out, indent=2, sort_keys=True) + "\n")
    print(json.dumps(out, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
