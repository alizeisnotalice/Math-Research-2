#!/usr/bin/env python3
"""Exact necessary-balance layout count for H7(x^3)H7(x^4).

H7 has one negative real root and three nonreal conjugate pairs.  The output
is a feasibility/cardinality audit, not a numerical root cache or a proof that
every listed selection has rational/nonnegative coefficients.
"""
from __future__ import annotations
import json
from math import comb
from pathlib import Path


def main() -> None:
    a, b, complex_blocks = 3, 4, 3
    terminal = a*b
    # At the negative real mother root: x^3-r gives one linear plus one
    # quadratic factor, while x^4-r gives two quadratic factors.
    real_a = {0: 1, 1: 1, 2: 1, 3: 1}
    real_b = {0: 1, 2: 2, 4: 1}
    by_E = {}
    for u in range(a+1):
        for v in range(b+1):
            E = b*u + a*v
            if 0 < E < terminal and v in real_b:
                assert E not in by_E
                by_E[E] = (u, v)
    rows = []
    for E, (u, v) in sorted(by_E.items()):
        complex_options = comb(a, u)*comb(b, v)
        real_options = real_a[u]*real_b[v]
        rows.append({'E': E, 'u': u, 'v': v, 'complex_block_options': complex_options,
                     'real_block_options': real_options,
                     'canonical_count': complex_options**complex_blocks * real_options})
    terminal_count = 2**complex_blocks  # Fix first complex block to all_a.
    total = sum(row['canonical_count'] for row in rows) + terminal_count
    assert rows == [
        {'E': 4, 'u': 1, 'v': 0, 'complex_block_options': 3, 'real_block_options': 1, 'canonical_count': 27},
        {'E': 6, 'u': 0, 'v': 2, 'complex_block_options': 6, 'real_block_options': 2, 'canonical_count': 432},
        {'E': 8, 'u': 2, 'v': 0, 'complex_block_options': 3, 'real_block_options': 1, 'canonical_count': 27},
        {'E': 10, 'u': 1, 'v': 2, 'complex_block_options': 18, 'real_block_options': 2, 'canonical_count': 11664},
    ]
    assert total == 12158
    result = {
        'candidate': 'H7(x^3)H7(x^4)', 'a': a, 'b': b, 'mother_degree': 7,
        'complex_parent_blocks': complex_blocks, 'real_parent_roots': 1,
        'root_group_layout': {'quadratic_groups': 24, 'real_linear_groups': 1, 'total_degree': 49},
        'balance': 'E=4*u+3*v equal across all seven mother roots (necessary only)',
        'real_block_selection_counts': {'a_lifts': real_a, 'b_lifts': real_b},
        'internal_rows': rows, 'terminal_E': terminal,
        'terminal_canonical_count': terminal_count, 'canonical_necessary_superset_count': total,
        'engine_feasibility': 'requires a real-linear-group extension; the current all-quadratic n-block cache is not directly valid',
        'claim_boundary': 'This is a necessary-balance search layout only; it is not a proof of rationality, nonnegativity, or a counterexample.'
    }
    Path('H7-a3b4-real-block-layout.json').write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps(result, indent=2))


if __name__ == '__main__': main()
