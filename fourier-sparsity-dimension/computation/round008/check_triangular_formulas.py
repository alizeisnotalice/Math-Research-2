#!/usr/bin/env python3
"""Exact truth-table audit for the Round-8 triangular-family formulas.

Enumerates k=1,2,3 only. All arithmetic is integer/F2; L values are exact
Fractions. Writes a self-describing JSON receipt instead of relying on stdout.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from fractions import Fraction
from pathlib import Path
from sys import version

HERE = Path(__file__).resolve().parent


def rank(rows: list[int]) -> int:
    piv: dict[int, int] = {}
    for x in rows:
        while x:
            b = x.bit_length() - 1
            if b in piv:
                x ^= piv[b]
            else:
                piv[b] = x
                break
    return len(piv)


def floor_log2_fraction(x: Fraction) -> int:
    d = x.numerator.bit_length() - x.denominator.bit_length()
    while Fraction(1 << d, 1) > x:
        d -= 1
    while Fraction(1 << (d + 1), 1) <= x:
        d += 1
    return d


def triangular_basis(k: int) -> tuple[list[int], list[list[int]], list[int]]:
    q = 1 << k
    vertical = [sum(1 << (u*q + v) for u in range(q)) for v in range(q)]
    levels: list[list[int]] = []
    for j in range(k):
        level = []
        for w in range(1 << j):
            mask = 0
            for u in range(q):
                for v in range(q):
                    if ((u >> j) & 1) and (v & ((1 << j) - 1)) == w:
                        mask |= 1 << (u*q + v)
            level.append(mask)
        levels.append(level)
    affine = [(1 << (q*q)) - 1]
    for j in range(k):
        affine.append(sum(1 << (u*q + v) for u in range(q) for v in range(q) if (u >> j) & 1))
        affine.append(sum(1 << (u*q + v) for u in range(q) for v in range(q) if (v >> j) & 1))
    return vertical + [x for lev in levels for x in lev], levels, affine


def all_words(basis: list[int]) -> list[int]:
    words = [0]
    for b in basis:
        words += [x ^ b for x in words]
    return words


def frac_text(x: Fraction) -> str:
    return f"{x.numerator}/{x.denominator}" if x.denominator != 1 else str(x.numerator)


def audit_k(k: int) -> dict[str, object]:
    q = 1 << k
    basis, _levels, affine = triangular_basis(k)
    words = all_words(basis)
    assert len(words) == 1 << (2*q - 1)
    minimum = [x for x in words if x.bit_count() == q]
    min_rank = rank(minimum)
    assert len(minimum) == 3*q
    assert min_rank == 3*q//2
    c_rank = rank(affine)
    cm_rank = rank(affine + minimum)
    assert c_rank == 2*k + 1
    assert c_rank + min_rank - cm_rank == k + 2
    assert len(basis) - cm_rank == q//2 - k

    thresholds = sorted(L for L in {Fraction(1), Fraction(3, 2), Fraction(2), Fraction(q, 4), Fraction(q, 2)} if L >= 1)
    threshold_rows = []
    for L in thresholds:
        d = floor_log2_fraction(L)
        J = max(0, k - 1 - d)
        low = [x for x in words if Fraction(x.bit_count(), q) <= L]
        low_rank = rank(low)
        c_low_rank = rank(affine + low)
        expected_low = 2*q - (1 << J)
        expected_v_residual = (1 << J) - 1
        expected_c_residual = (1 << J) - 1 - J
        assert low_rank == expected_low
        assert len(basis) - low_rank == expected_v_residual
        assert len(basis) - c_low_rank == expected_c_residual
        threshold_rows.append({
            "L": frac_text(L), "J": J, "low_span_rank": low_rank,
            "V_mod_low": len(basis) - low_rank,
            "V_mod_C_plus_low": len(basis) - c_low_rank,
        })
    return {
        "k": k, "q": q, "v_dim": len(basis), "minimum_count": len(minimum),
        "minimum_span_rank": min_rank, "c_dim": c_rank,
        "c_intersection_minimum": c_rank + min_rank - cm_rank,
        "V_mod_C_plus_minimum": len(basis) - cm_rank,
        "thresholds": threshold_rows,
    }


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", type=Path, default=HERE / "triangular-formulas-results.json")
    ns = ap.parse_args()
    source = Path(__file__)
    result = {
        "status": "PASS", "method": "complete truth-table F2 enumeration",
        "python": version, "source": str(source),
        "source_sha256": hashlib.sha256(source.read_bytes()).hexdigest(),
        "k_values": [1, 2, 3], "results": [audit_k(k) for k in (1, 2, 3)],
    }
    ns.out.parent.mkdir(parents=True, exist_ok=True)
    ns.out.write_text(json.dumps(result, ensure_ascii=False, indent=2) + "\n")
    print(json.dumps({"status": "PASS", "out": str(ns.out)}, ensure_ascii=False))


if __name__ == "__main__":
    main()
