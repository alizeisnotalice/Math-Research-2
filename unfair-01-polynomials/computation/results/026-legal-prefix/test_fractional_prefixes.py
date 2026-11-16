#!/usr/bin/env python3
"""Independent exhaustive test for fractional_prefixes.py, standard library only."""
from __future__ import annotations

import hashlib
import itertools
import json
from pathlib import Path

import fractional_prefixes as subject


def independently_legal(code: str, d: int) -> bool:
    pairs = {'0': (0, 0), '1': (1, 0), '2': (0, 1)}
    a = [1] + [pairs[c][0] for c in code]
    b = [1] + [pairs[c][1] for c in code]
    for k in range(d):
        coefficient = 0
        for i in range(k + 1):
            coefficient += a[i] * b[k-i]
        if coefficient not in (0, 1):
            return False
    return sum(a[i] * b[d-i] for i in range(1, d)) == 0


def independently_filter_reasons(code: str, d: int) -> list[str]:
    pairs = {'0': (0, 0), '1': (1, 0), '2': (0, 1)}
    a = [1] + [pairs[c][0] for c in code]
    b = [1] + [pairs[c][1] for c in code]
    product = []
    for k in range(2*d-1):
        product.append(sum(a[i]*b[k-i] for i in range(max(0, k-len(b)+1), min(k, len(a)-1)+1)))
    reasons = []
    if any(value > 1 for value in product):
        reasons.append('prefix_product_coefficient_gt_1')
    if any((product[d+i] if d+i < len(product) else 0) >= 1 and a[i]+b[i] > 0 for i in range(1, d)):
        reasons.append('forced_positive_fractional_cross_at_d_plus_i')
    return reasons


def main() -> None:
    target = Path(__file__).resolve().parent / 'fractional-prefixes-d1-d9.json'
    saved = json.loads(target.read_text())
    assert saved['max_d'] == 9
    for row in saved['rows']:
        d = row['d']
        expected = [''.join(chars) for chars in itertools.product('012', repeat=d-1)
                    if independently_legal(''.join(chars), d)]
        assert row['valid_codes'] == expected == subject.enumerate_valid_codes(d)
        assert row['valid_count'] == len(expected)
        assert row['code_sha256'] == hashlib.sha256(json.dumps(expected, separators=(',', ':')).encode()).hexdigest()
        assert all(subject.legal_code(code, d) for code in expected)
        assert all((code in expected) == subject.legal_code(code, d)
                   for code in (''.join(chars) for chars in itertools.product('012', repeat=d-1)))
        survivors = [code for code in expected if not independently_filter_reasons(code, d)]
        strict = row['strict_necessary_filter']
        assert strict['survivor_codes'] == survivors
        assert strict['survivor_count'] == len(survivors)
        assert strict['survivor_code_sha256'] == hashlib.sha256(json.dumps(survivors, separators=(',', ':')).encode()).hexdigest()
        assert strict['rejected_by_reason_counts'] == {
            reason: sum(reason in independently_filter_reasons(code, d) for code in expected)
            for reason in ('prefix_product_coefficient_gt_1', 'forced_positive_fractional_cross_at_d_plus_i')}
    receipt = {'status': 'independent_exhaustive_test_passed', 'd_range': [1, 9],
               'report_sha256': hashlib.sha256(target.read_bytes()).hexdigest(),
               'test_source_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
               'float_arithmetic_used': False}
    out = Path(__file__).resolve().parent / 'fractional-prefixes-independent-test.json'
    out.write_text(json.dumps(receipt, indent=2) + '\n')
    print(json.dumps(receipt, indent=2))


if __name__ == '__main__':
    main()
