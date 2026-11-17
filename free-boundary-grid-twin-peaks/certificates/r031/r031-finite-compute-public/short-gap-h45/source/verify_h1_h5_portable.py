#!/usr/bin/env python3
"""Portable finite h1..h5 replay against r029 and r030 full L histograms.

This program reads saved four-prime checkpoints; it deliberately does not run
the n=7 computation.  Its reference locations are explicit arguments so the
source has no private absolute-path dependency.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from math import prod
from pathlib import Path

PRIMES = [9223372036854775783, 9223372036854775643,
          9223372036854775549, 9223372036854775507]
SOURCE_SHA256 = "08c31508f6d3cb0ae50f74973aa372a32e3bac8439d9d6bfd51d4bf6c0312fc8"
N5_FILES = {0: "tail-full-n5-v0.json", 1: "tail-full-n5-v1.json",
            2: "tail-full-n5-v2.json", 6: "tail-full-n5-v6.json",
            7: "tail-full-n5-v7.json", 12: "tail-full-n5-center.json"}


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def crt(residues: list[int]) -> int:
    modulus = prod(PRIMES)
    return sum(x * (modulus // p) * pow(modulus // p, -1, p)
               for x, p in zip(residues, PRIMES)) % modulus


def checkpoint_values(directory: Path, n: int, vertex: int) -> list[int]:
    rows = [json.loads((directory / f"n{n}-v{vertex}-p{i}.json").read_text())
            for i in range(4)]
    for i, row in enumerate(rows):
        assert row["n"] == n and row["v"] == vertex
        assert row["prime_index"] == i and int(row["prime"]) == PRIMES[i]
        assert row["algorithm"] == "short-gap-h1-to-h5-modular-layer-rolling"
    return [crt([int(row[f"h{j}_residue"]) for row in rows]) for j in range(1, 6)]


def hole_map(path: Path) -> dict[int, dict]:
    document = json.loads(path.read_text())
    return {int(item["v"]): item for item in document["holes"]}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--checkpoints", type=Path, required=True)
    parser.add_argument("--r029-tail", type=Path, required=True,
                        help="r029/r029-tail-full-n5 directory")
    parser.add_argument("--r030-geometry", type=Path, required=True,
                        help="r030/r030-gap-geometry directory")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    assert sha256(args.source) == SOURCE_SHA256, "frozen source SHA-256 mismatch"
    r29_n3 = args.r029_tail / "tail-full-n3.json"
    r29_n4 = args.r029_tail / "tail-full-n4.json"
    r29_n5 = {v: args.r029_tail / name for v, name in N5_FILES.items()}
    r30_n3 = args.r030_geometry / "gap-n3-all.json"
    r30_n5 = {v: args.r030_geometry / f"gap-n5-v{v}.json" for v in N5_FILES}
    for path in [r29_n3, r29_n4, r30_n3, *r29_n5.values(), *r30_n5.values()]:
        assert path.is_file(), f"missing reference: {path}"

    checks: list[dict] = []
    references: dict[str, dict[str, str]] = {}
    r29_by_n = {3: hole_map(r29_n3), 4: hole_map(r29_n4)}
    # r029 is the independent complete L-histogram source for all 31 points.
    for n in (3, 4):
        references[f"r029_n{n}"] = {"sha256": sha256((r29_n3, r29_n4)[n - 3])}
        for v, item in r29_by_n[n].items():
            got = checkpoint_values(args.checkpoints, n, v)
            want = list(map(int, item["histogram_L"][:5]))
            assert got == want, (n, v, got, want)
            checks.append({"n": n, "v": v, "h": list(map(str, got)),
                           "r029_complete_histogram_match": True})
    # r030 separately carries complete n=3 and n=5 geometry tables.  Read
    # them and require equality with r029 before using the same five entries.
    r30_n3_holes = hole_map(r30_n3)
    for v, item in r29_by_n[3].items():
        assert item["histogram_L"] == r30_n3_holes[v]["histogram_L"], ("r030", 3, v)
    references["r030_n3"] = {"sha256": sha256(r30_n3)}
    for v, r29_path in r29_n5.items():
        r29 = hole_map(r29_path)[v]
        r30_path = r30_n5[v]
        r30 = hole_map(r30_path)[v]
        assert r29["histogram_L"] == r30["histogram_L"], ("r030", 5, v)
        assert r29["A2_v_from_fibres"] == r30["A2_v_from_fibres"], ("r030 A2", v)
        assert r29["B_v_from_D"] == r30["B_v_from_D"], ("r030 B", v)
        got = checkpoint_values(args.checkpoints, 5, v)
        want = list(map(int, r29["histogram_L"][:5]))
        assert got == want, (5, v, got, want)
        checks.append({"n": 5, "v": v, "h": list(map(str, got)),
                       "r029_r030_complete_histogram_match": True})
        references[f"r029_n5_v{v}"] = {"sha256": sha256(r29_path)}
        references[f"r030_n5_v{v}"] = {"sha256": sha256(r30_path)}
    assert len(checks) == 31
    result = {
        "status": "passed",
        "scope": "saved n=3,4,5 checkpoint CRT replay against complete independent r029/r030 L histograms; no n=7 computation or asymptotic inference",
        "source_sha256": SOURCE_SHA256,
        "checkpoint_count": 124,
        "point_count": len(checks),
        "crt_modulus_bits": prod(PRIMES).bit_length(),
        "references": references,
        "checks": checks,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps({k: result[k] for k in ("status", "checkpoint_count", "point_count", "crt_modulus_bits", "scope")}, indent=2))


if __name__ == "__main__":
    main()
