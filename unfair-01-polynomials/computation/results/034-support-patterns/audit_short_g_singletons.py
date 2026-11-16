#!/usr/bin/env python3
"""Independent one-core audit for the q <= 20 short-G support ledger.

This program deliberately does not import the producer.  It enumerates every
odd q-bit integer mask, so bit 0 of T is present by construction.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import time
from collections import Counter
from math import gcd
from pathlib import Path


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def support_from_mask(mask: int, q: int) -> tuple[int, ...]:
    """Return T union (T+q), with T represented by q low bits of mask."""
    return tuple(i for i in range(2 * q) if (mask >> (i % q)) & 1)


def c_injective(mask: int, p: int, q: int) -> bool:
    """Test the exact no-p-difference criterion on S=T union (T+q)."""
    gmask = mask | (mask << q)
    return (gmask & (gmask >> p)) == 0


def no_mod_p_singleton(mask: int, p: int, q: int) -> bool:
    counts = [0] * p
    gmask = mask | (mask << q)
    for i in range(2 * q):
        if (gmask >> i) & 1:
            counts[i % p] += 1
    return all(c != 1 for c in counts)


def normalized_record(p: int, q: int, mask: int) -> bytes:
    t = [i for i in range(q) if (mask >> i) & 1]
    s = list(support_from_mask(mask, q))
    return json.dumps(
        {"G_support": s, "T": t, "p": p, "q": q},
        sort_keys=True, separators=(",", ":"),
    ).encode() + b"\n"


def normalized_key(key: tuple[int, int, tuple[int, ...], tuple[int, ...]]) -> bytes:
    p, q, t, s = key
    return json.dumps({"G_support": list(s), "T": list(t), "p": p, "q": q},
                      sort_keys=True, separators=(",", ":")).encode() + b"\n"


def record_key(row: dict, p: int, q: int) -> tuple[int, int, tuple[int, ...], tuple[int, ...]]:
    if set(row) != {"T", "G_support"}:
        raise ValueError(f"unexpected survivor fields for {(p, q)}: {sorted(row)}")
    t, s = row["T"], row["G_support"]
    if (any(type(x) is not int for x in t) or
            any(type(x) is not int for x in s)):
        raise ValueError(f"non-integer support datum for {(p, q)}")
    if t != sorted(set(t)) or s != sorted(set(s)):
        raise ValueError(f"non-canonical survivor datum for {(p, q)}")
    if not t or t[0] != 0 or any(x < 0 or x >= q for x in t):
        raise ValueError(f"invalid T for {(p, q)}")
    expected_s = sorted(set(t) | {x + q for x in t})
    if s != expected_s:
        raise ValueError(f"G_support is not T union (T+q) for {(p, q)}")
    return p, q, tuple(t), tuple(s)


def archive_mask(row: dict, p: int, q: int) -> int:
    mask = 0
    for i in row["T"]:
        mask |= 1 << i
    return mask


def expected_pairs(max_q: int) -> list[tuple[int, int]]:
    return [(p, q) for q in range(3, max_q + 1) for p in range(2, q)
            if gcd(p, q) == 1]


def direct_c_map_injective(mask: int, p: int, q: int) -> bool:
    """Slow literal check of (a,b)->pa+qb, used only for boundary samples."""
    s = support_from_mask(mask, q)
    values = [p * a + q * b for a in s for b in s]
    return len(values) == len(set(values))


def main() -> None:
    parser=argparse.ArgumentParser()
    parser.add_argument("--input",type=Path,required=True)
    parser.add_argument("--record",type=Path,required=True)
    args=parser.parse_args()
    input_path=args.input
    started = time.monotonic()
    data = json.loads(input_path.read_text())
    if data.get("schema") != "short-G-singleton-support-v1":
        raise ValueError("unexpected input schema")
    if data.get("range") != {"q<=": 20, "coprime_1<p<q": True,
                             "q_ge_2p_only": False}:
        raise ValueError("unexpected claimed coverage range")
    if type(data.get("rows")) is not list:
        raise ValueError("rows must be a list")
    pairs = expected_pairs(20)
    archive_by_pair: dict[tuple[int, int], dict] = {}
    for row in data["rows"]:
        if not isinstance(row, dict) or type(row.get("p")) is not int or type(row.get("q")) is not int:
            raise ValueError("bad p/q row")
        pair = row["p"], row["q"]
        if pair in archive_by_pair:
            raise ValueError(f"duplicated p/q row: {pair}")
        archive_by_pair[pair] = row
    if set(archive_by_pair) != set(pairs) or len(archive_by_pair) != len(pairs):
        raise ValueError("p/q coverage is not exactly the coprime grid")

    archive_stream = hashlib.sha256()
    generated_stream = hashlib.sha256()
    totals = Counter()
    per_pair: list[dict] = []
    q_less_2p_zero: list[tuple[int, int]] = []
    direct_samples = 0
    for p, q in pairs:
        row = archive_by_pair[p, q]
        survivors = row.get("no_mod_p_singleton_survivors")
        if type(survivors) is not list or type(row.get("T_patterns_binary_C_compatible")) is not int or type(row.get("no_mod_p_singleton_count")) is not int:
            raise ValueError(f"bad counts/list for {(p,q)}")
        archived: set[tuple[int, int, tuple[int, ...], tuple[int, ...]]] = set()
        for survivor in survivors:
            key = record_key(survivor, p, q)
            if key in archived:
                raise ValueError(f"duplicated archived survivor {key}")
            archived.add(key)
        if len(archived) != row["no_mod_p_singleton_count"]:
            raise ValueError(f"archived survivor count mismatch for {(p,q)}")

        generated: set[tuple[int, int, tuple[int, ...], tuple[int, ...]]] = set()
        binary_count = 0
        for mask in range(1, 1 << q, 2):
            ci = c_injective(mask, p, q)
            # Literal product-map checks at both ends and the first surviving
            # pattern independently validate the bit criterion during this run.
            if mask in (1, (1 << q) - 1):
                if ci != direct_c_map_injective(mask, p, q):
                    raise AssertionError(f"criterion mismatch at boundary {(p,q,mask)}")
                direct_samples += 1
            if not ci:
                continue
            binary_count += 1
            if no_mod_p_singleton(mask, p, q):
                t = tuple(i for i in range(q) if (mask >> i) & 1)
                s = support_from_mask(mask, q)
                key = p, q, t, s
                if key in generated:
                    raise AssertionError("duplicate generated mask")
                generated.add(key)
                if direct_samples < 3 * len(pairs):
                    if not direct_c_map_injective(mask, p, q):
                        raise AssertionError(f"criterion mismatch at survivor {(p,q,mask)}")
                    direct_samples += 1
        if binary_count != row["T_patterns_binary_C_compatible"]:
            raise ValueError(f"binary-compatible count mismatch for {(p,q)}")
        if generated != archived:
            raise ValueError(f"concrete survivor omission/addition for {(p,q)}")
        # The digest order is defined here, rather than inherited from either
        # program's DFS/mask traversal order.
        for key in sorted(archived):
            archive_stream.update(normalized_key(key))
        for key in sorted(generated):
            generated_stream.update(normalized_key(key))
        if q < 2 * p:
            if generated:
                raise ValueError(f"q<2p has a survivor for {(p,q)}")
            q_less_2p_zero.append((p, q))
        totals["binary"] += binary_count
        totals["survivors"] += len(generated)
        totals["nonempty_pairs"] += bool(generated)
        per_pair.append({"p": p, "q": q, "binary": binary_count,
                         "survivors": len(generated),
                         "q_less_2p": q < 2 * p})

    if data.get("summary") != {
        "pairs": len(pairs), "binary_C_compatible_T_patterns": totals["binary"],
        "no_singleton_survivors": totals["survivors"],
        "pairs_with_survivor": totals["nonempty_pairs"],
    }:
        raise ValueError("top-level summary mismatch")
    receipt = {
        "audit": "independent-short-G-mask-enumeration-v1",
        "scope": "Necessary support-level conditions only; not a factorization or P33 counterexample.",
        "input": str(input_path), "input_sha256": sha256_file(input_path),
        "audit_source_sha256": sha256_file(Path(__file__)),
        "method": "sequential odd q-bit masks; no producer import; gmask & (gmask>>p)==0; residue counts modulo p",
        "one_cpu_environment": {k: os.environ.get(k) for k in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS")},
        "pairs": len(pairs), "totals": dict(totals), "q_less_2p_zero_pairs": q_less_2p_zero,
        "normalized_archive_survivor_stream_sha256": archive_stream.hexdigest(),
        "normalized_generated_survivor_stream_sha256": generated_stream.hexdigest(),
        "direct_product_map_samples": direct_samples,
        "per_pair": per_pair,
        "elapsed_seconds": time.monotonic() - started,
        "status": "complete",
    }
    args.record.write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"status": receipt["status"], "totals": receipt["totals"],
                      "elapsed_seconds": receipt["elapsed_seconds"],
                      "archive_digest": receipt["normalized_archive_survivor_stream_sha256"],
                      "generated_digest": receipt["normalized_generated_survivor_stream_sha256"]}, sort_keys=True))


if __name__ == "__main__":
    main()
