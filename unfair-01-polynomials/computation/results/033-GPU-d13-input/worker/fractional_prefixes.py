#!/usr/bin/env python3
"""Exact low-prefix enumeration for a first fractional coefficient at d.

Only Python's standard library is used.  A code digit 0,1,2 denotes
(a_i,b_i)=(0,0),(1,0),(0,1), respectively, for i=1,...,d-1.
"""
from __future__ import annotations

import hashlib
import itertools
import json
from pathlib import Path


PAIR_OPTIONS = ((0, 0), (1, 0), (0, 1))


def convolution_prefix(a: list[int], b: list[int], end: int) -> list[int]:
    return [sum(a[i] * b[k-i] for i in range(max(0, k-len(b)+1), min(k, len(a)-1) + 1))
            for k in range(end + 1)]


def decode(code: str) -> tuple[list[int], list[int]]:
    a, b = [1], [1]
    for char in code:
        ai, bi = PAIR_OPTIONS[int(char)]
        a.append(ai)
        b.append(bi)
    return a, b


def legal_code(code: str, d: int) -> bool:
    if len(code) != d - 1 or any(char not in '012' for char in code):
        return False
    a, b = decode(code)
    low = convolution_prefix(a, b, d - 1)
    if any(value not in (0, 1) for value in low):
        return False
    cross = sum(a[i] * b[d-i] for i in range(1, d))
    return cross == 0


def enumerate_valid_codes(d: int) -> list[str]:
    if d < 1:
        raise ValueError('d must be positive')
    return [''.join(chars) for chars in itertools.product('012', repeat=d-1)
            if legal_code(''.join(chars), d)]


def strict_extension_filter_reasons(code: str, d: int) -> list[str]:
    """Necessary (not sufficient) nonextension tests for every 0<t<1."""
    a, b = decode(code)
    prefix_product = convolution_prefix(a, b, 2*d-2)
    reasons = []
    if any(value > 1 for value in prefix_product):
        reasons.append('prefix_product_coefficient_gt_1')
    for i in range(1, d):
        # At degree d+i, the fixed fractional pair contributes
        # t*b_i+(1-t)*a_i, strictly positive whenever a_i+b_i>0.
        baseline = prefix_product[d+i] if d+i < len(prefix_product) else 0
        if baseline >= 1 and a[i] + b[i] > 0:
            reasons.append('forced_positive_fractional_cross_at_d_plus_i')
            break
    return reasons


def report(max_d: int) -> dict[str, object]:
    rows = []
    for d in range(1, max_d + 1):
        codes = enumerate_valid_codes(d)
        survivors = [code for code in codes if not strict_extension_filter_reasons(code, d)]
        reason_counts = {
            'prefix_product_coefficient_gt_1': sum('prefix_product_coefficient_gt_1' in strict_extension_filter_reasons(code, d) for code in codes),
            'forced_positive_fractional_cross_at_d_plus_i': sum('forced_positive_fractional_cross_at_d_plus_i' in strict_extension_filter_reasons(code, d) for code in codes),
        }
        digest = hashlib.sha256(json.dumps(codes, separators=(',', ':')).encode()).hexdigest()
        survivor_digest = hashlib.sha256(json.dumps(survivors, separators=(',', ':')).encode()).hexdigest()
        rows.append({'d': d, 'raw_prefix_codes': 3**(d-1), 'valid_count': len(codes),
                     'code_sha256': digest, 'valid_codes': codes,
                     'strict_necessary_filter': {
                         'survivor_count': len(survivors), 'survivor_code_sha256': survivor_digest,
                         'survivor_codes': survivors, 'rejected_by_reason_counts': reason_counts,
                         'meaning': 'Survival is necessary only; it does not prove a global extension.'}})
    return {
        'model': 'A0=B0=1; i<d pairs in {(0,0),(1,0),(0,1)}; low convolution binary; d-cross zero',
        'encoding': {'0': [0, 0], '1': [1, 0], '2': [0, 1]},
        'max_d': max_d, 'rows': rows,
        'claim_boundary': 'Exact low-prefix feasibility is separate from the strict necessary nonextension filter. Neither survival nor rejection by unlisted conditions claims a complete nonnegative factorization.',
    }


def main() -> None:
    output = Path(__file__).resolve().parent / 'fractional-prefixes-d1-d9.json'
    value = report(9)
    output.write_text(json.dumps(value, indent=2) + '\n')
    print(json.dumps({'output': str(output), 'rows': [(row['d'], row['valid_count']) for row in value['rows']],
                      'output_sha256': hashlib.sha256(output.read_bytes()).hexdigest()}, indent=2))


if __name__ == '__main__':
    main()
