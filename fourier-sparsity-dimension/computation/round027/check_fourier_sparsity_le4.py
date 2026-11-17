#!/usr/bin/env python3
"""Exact exhaustive check of Boolean functions on F_2^4 using integer FWHT."""

from __future__ import annotations

import hashlib
import itertools
import json
import platform
import sys
import time
from pathlib import Path

N = 4
Q = 1 << N
OUT = Path(__file__).with_name("fourier_sparsity_le4_result.json")


def fwht(values: list[int]) -> list[int]:
    values = values[:]
    block = 1
    while block < len(values):
        for start in range(0, len(values), 2 * block):
            for offset in range(block):
                a = values[start + offset]
                b = values[start + offset + block]
                values[start + offset] = a + b
                values[start + offset + block] = a - b
        block *= 2
    return values


def gf2_rank(rows: list[int]) -> int:
    rows = rows[:]
    rank = 0
    for bit in range(N - 1, -1, -1):
        pivot = next((i for i in range(rank, len(rows)) if (rows[i] >> bit) & 1), None)
        if pivot is None:
            continue
        rows[rank], rows[pivot] = rows[pivot], rows[rank]
        for i in range(len(rows)):
            if i != rank and ((rows[i] >> bit) & 1):
                rows[i] ^= rows[rank]
        rank += 1
    return rank


def anf_coefficients(table_mask: int) -> list[int]:
    coeff = [(table_mask >> x) & 1 for x in range(Q)]
    for coordinate in range(N):
        step = 1 << coordinate
        for monomial in range(Q):
            if monomial & step:
                coeff[monomial] ^= coeff[monomial ^ step]
    return coeff


def quadratic_rank(coeff: list[int]) -> int:
    # Rows of the polar alternating matrix of the degree-two ANF part.
    rows = [0] * N
    for i in range(N):
        for j in range(i + 1, N):
            if coeff[(1 << i) | (1 << j)]:
                rows[i] |= 1 << j
                rows[j] |= 1 << i
    return gf2_rank(rows)


def is_affine_two_flat(support: list[int]) -> bool:
    if len(support) != 4:
        return False
    origin = support[0]
    return gf2_rank([point ^ origin for point in support[1:]]) == 2


def main() -> None:
    started = time.perf_counter()
    support_histogram: dict[int, int] = {}
    sparse_count = 0
    support_one_count = 0
    support_four_count = 0
    failures: list[dict[str, object]] = []
    for truth_table in range(1 << Q):
        signs = [1 if ((truth_table >> x) & 1) == 0 else -1 for x in range(Q)]
        spectrum = fwht(signs)
        support = [frequency for frequency, value in enumerate(spectrum) if value != 0]
        sparsity = len(support)
        support_histogram[sparsity] = support_histogram.get(sparsity, 0) + 1
        if sparsity <= 4:
            sparse_count += 1
            if sparsity == 1:
                support_one_count += 1
                coeff = anf_coefficients(truth_table)
                degree = max((monomial.bit_count() for monomial, bit in enumerate(coeff) if bit), default=0)
                if abs(spectrum[support[0]]) != 16 or degree > 1:
                    failures.append(
                        {
                            "truth_table": truth_table,
                            "failure": "support-one structure",
                            "amplitude": spectrum[support[0]],
                            "anf_degree": degree,
                        }
                    )
            elif sparsity == 4:
                support_four_count += 1
                coeff = anf_coefficients(truth_table)
                degree = max((monomial.bit_count() for monomial, bit in enumerate(coeff) if bit), default=0)
                qrank = quadratic_rank(coeff)
                if not (
                    is_affine_two_flat(support)
                    and all(abs(spectrum[frequency]) == 8 for frequency in support)
                    and degree == 2
                    and qrank == 2
                ):
                    failures.append(
                        {
                            "truth_table": truth_table,
                            "failure": "support-four structure",
                            "support": support,
                            "amplitudes": [spectrum[frequency] for frequency in support],
                            "anf_degree": degree,
                            "quadratic_rank": qrank,
                        }
                    )
            else:
                failures.append({"truth_table": truth_table, "failure": "unexpected sparsity", "sparsity": sparsity})
    result = {
        "claim_scope": "All 2^16 Boolean truth tables on F_2^4; exact integer FWHT.",
        "python": sys.version,
        "platform": platform.platform(),
        "enumerated_functions": 1 << Q,
        "support_histogram": {str(key): support_histogram[key] for key in sorted(support_histogram)},
        "sparsity_le_4_count": sparse_count,
        "sparsity_1_count": support_one_count,
        "sparsity_4_count": support_four_count,
        "failure_count": len(failures),
        "failures": failures[:10],
        "elapsed_seconds": time.perf_counter() - started,
        "script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    }
    OUT.write_text(json.dumps(result, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result, ensure_ascii=False, indent=2))
    if failures:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
