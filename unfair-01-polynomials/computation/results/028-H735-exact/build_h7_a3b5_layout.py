#!/usr/bin/env python3
"""Exact necessary root-selection ledger for H7(x^3) H7(x^5).

The bit masks enumerate only the Galois/balance necessary superset; they do
not establish that a listed selection gives real, nonnegative coefficients.
No roots are approximated here.
"""
from __future__ import annotations

import hashlib
import itertools
import json
from collections import Counter
from pathlib import Path


ROOT = Path(__file__).resolve().parent
H = [1, 1, 1, 0, 0, 0, 1, 1]  # ascending coefficients
A, B = 3, 5
TERMINAL = A * B
MOTHER_CERT = Path('/Volumes/data/CodexProjects/T2/unfair-01-polynomials/computation/results/025-Hd-frontier/rev3/H7-mother-certificate.json')
MOTHER_MODULAR = Path('/Volumes/data/CodexProjects/T2/unfair-01-polynomials/computation/results/025-Hd-frontier/rev3/H7-mother-modular-verification.json')


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def substituted_support(power: int) -> list[int]:
    return [i * power for i, c in enumerate(H) if c]


def exact_binary_target() -> dict[str, object]:
    left, right = substituted_support(A), substituted_support(B)
    coefficients = Counter(i + j for i in left for j in right)
    assert len(coefficients) == 25 and set(coefficients.values()) == {1}
    assert min(coefficients) == 0 and max(coefficients) == 56
    return {'H_ascending': H, 'a': A, 'b': B, 'degree': 56,
            'left_support': left, 'right_support': right,
            'support': sorted(coefficients), 'term_count': len(coefficients),
            'all_coefficients_are_exactly_one': True}


def mask(indices: tuple[int, ...], size: int) -> str:
    return ''.join('1' if i in indices else '0' for i in range(size))


def complex_options() -> dict[int, list[str]]:
    """One nonreal parent pair: 3 a-quadratics and 5 b-quadratics."""
    out: dict[int, list[str]] = {}
    for u in range(A + 1):
        for v in range(B + 1):
            e = B * u + A * v
            for ia in itertools.combinations(range(A), u):
                for ib in itertools.combinations(range(A, A + B), v):
                    out.setdefault(e, []).append(mask(ia + ib, A + B))
    for values in out.values():
        values.sort()
    return out


def real_options() -> dict[int, list[str]]:
    """Real parent: a linear+quadratic and b linear+two quadratics."""
    out: dict[int, list[str]] = {}
    for bits in itertools.product((0, 1), repeat=5):
        eps_a, eta_a, eps_b, w0, w1 = bits
        e = B * eps_a + 2 * B * eta_a + A * eps_b + 2 * A * (w0 + w1)
        out.setdefault(e, []).append(''.join(map(str, bits)))
    for values in out.values():
        values.sort()
    return out


def complement(bitmask: str) -> str:
    return ''.join('0' if bit == '1' else '1' for bit in bitmask)


def selected_indices(bitmask: str) -> list[int]:
    return [i for i, bit in enumerate(bitmask) if bit == '1']


def block_index_complex(bitmask: str) -> int:
    return B * bitmask[:A].count('1') + A * bitmask[A:].count('1')


def block_index_real(bitmask: str) -> int:
    eps_a, eta_a, eps_b, w0, w1 = map(int, bitmask)
    return B * eps_a + 2 * B * eta_a + A * eps_b + 2 * A * (w0 + w1)


def full_mask(c0: str, c1: str, c2: str, r: str) -> str:
    return c0 + c1 + c2 + r


def iter_canonical(copts: dict[int, list[str]], ropts: dict[int, list[str]]):
    """All actual 29-bit masks modulo simultaneous factor complementation."""
    for e in range(1, TERMINAL):
        if e not in copts or e not in ropts:
            continue
        for c0, c1, c2, r in itertools.product(copts[e], copts[e], copts[e], ropts[e]):
            yield e, full_mask(c0, c1, c2, r)
    # At E=15 complement stays in this E layer. Fix the first complex block
    # to its all-a mask 11100000; complement has its first block all-b.
    all_a = '1' * A + '0' * B
    assert block_index_complex(all_a) == TERMINAL
    for c1, c2, r in itertools.product(copts[TERMINAL], copts[TERMINAL], ropts[TERMINAL]):
        yield TERMINAL, full_mask(all_a, c1, c2, r)


def verify_entry(e: int, bits: str, terminal: bool) -> None:
    assert len(bits) == 29
    blocks = (bits[:8], bits[8:16], bits[16:24], bits[24:])
    assert [block_index_complex(item) for item in blocks[:3]] == [e, e, e]
    assert block_index_real(blocks[3]) == e
    comp = complement(bits)
    cblocks = (comp[:8], comp[8:16], comp[16:24], comp[24:])
    assert [block_index_complex(item) for item in cblocks[:3]] == [30-e, 30-e, 30-e]
    assert block_index_real(cblocks[3]) == 30-e
    if terminal:
        assert bits[:8] == '11100000' and comp[:8] == '00011111'


def main() -> None:
    target = exact_binary_target()
    copts, ropts = complex_options(), real_options()
    assert copts[15] == ['00011111', '11100000']
    expected_rows = {
        3: (5, 1, 125), 5: (3, 1, 27), 6: (10, 2, 2000),
        8: (15, 1, 3375), 9: (10, 2, 2000), 10: (3, 1, 27),
        11: (30, 2, 54000), 12: (5, 1, 125), 13: (15, 1, 3375),
        14: (30, 2, 54000)}
    rows = []
    for e, (c, r, count) in expected_rows.items():
        assert len(copts[e]) == c and len(ropts[e]) == r and c**3*r == count
        rows.append({'E': e, 'complex_block_option_count': c, 'real_block_option_count': r,
                     'actual_canonical_mask_count': count})
    assert len(copts[15]) == 2 and len(ropts[15]) == 2

    ledger = ROOT/'H7-a3b5-canonical-masks.jsonl'
    count_by_e: Counter[int] = Counter()
    all_masks: set[str] = set()
    with ledger.open('w', encoding='utf-8') as handle:
        for e, bits in iter_canonical(copts, ropts):
            verify_entry(e, bits, e == TERMINAL)
            assert bits not in all_masks
            all_masks.add(bits); count_by_e[e] += 1
            handle.write(json.dumps({'E': e, 'mask': bits, 'selected_group_indices': selected_indices(bits)}, separators=(',', ':'))+'\n')
    assert count_by_e == Counter({row['E']: row['actual_canonical_mask_count'] for row in rows} | {15: 8})
    assert len(all_masks) == 119062

    group_order = []
    for c in range(3):
        group_order += [f'complex_{c}:a_quadratic_lift_{j}' for j in range(3)]
        group_order += [f'complex_{c}:b_quadratic_lift_{j}' for j in range(5)]
    group_order += ['real:a_linear_lift', 'real:a_quadratic_lift', 'real:b_linear_lift',
                    'real:b_quadratic_lift_0', 'real:b_quadratic_lift_1']
    assert len(group_order) == 29
    summary = {
        'format': 'H7-a3b5-exact-necessary-root-layout-v1',
        'target': target,
        'mother_root_schema': {
            'H7': '1+x+x^2+x^6+x^7', 'nonreal_conjugate_parent_pairs': 3,
            'negative_real_parent_roots': 1,
            'certificate_source_sha256': digest(MOTHER_CERT),
            'modular_certificate_source_sha256': digest(MOTHER_MODULAR),
            'source_claim_used': 'For a monic real factor A with constant coefficient 1, the selected-root product P=(-1)^deg(A) is in {+1,-1}, hence P^15 is in Q*. The existing S7 certificate is used only for the resulting common 15th-power rationality necessary balance; A itself is not assumed rational.'},
        'groups': {'bit_order': group_order, 'quadratic_group_count': 27, 'linear_group_count': 2,
                   'total_bit_count': 29,
                   'root_group_assignment_convention': 'bit 1 selects that root group into factor A and bit 0 selects it into factor B. The bitwise complement exchanges A and B; it does not label individual groups by a sign of the overall root product P.',
                   'complex_block': '3 a quadratic lifts + 5 b quadratic lifts',
                   'real_block': 'a linear + a quadratic + b linear + 2 b quadratics'},
        'balance': {'complex_parent_block_E': '5*u + 3*v',
                    'real_parent_block_E': '5*eps_a + 10*eta_a + 3*eps_b + 6*w0 + 6*w1',
                    'common_E_range': [0, 30], 'trivial_E_excluded': [0, 30],
                    'complement_rule': 'simultaneous A/B exchange (bitwise complement) sends E to 30-E; it is not a per-group root-product-sign operation.',
                    'terminal_rule': 'At E=15 retain only masks whose first complex block is all-a (11100000).'},
        'internal_E_rows': rows,
        'terminal_E': {'E': 15, 'raw_mask_count_before_complement_quotient': 16,
                       'canonical_mask_count_after_first_complex_all_a': 8},
        'canonical_actual_mask_count': len(all_masks),
        'ledger': {'path': str(ledger), 'sha256': digest(ledger), 'record_count': len(all_masks),
                   'every_record_verified_blockwise_and_under_complement': True},
        'reproduction': {'generator': 'build_h7_a3b5_layout.py', 'compact_compressor': 'compress_h7_a3b5_ledger.py',
                           'public_archival_policy': 'Do not archive the 13 MiB uncompressed ledger; regenerate it or use the deterministic gzip representation under 10 MiB.'},
        'claim_boundary': 'The 119062 records are an exact necessary root-layout superset only. They are not an exact factorization enumeration, a nonnegativity proof, or a counterexample.'}
    out = ROOT/'H7-a3b5-layout-summary.json'
    out.write_text(json.dumps(summary, indent=2)+'\n', encoding='utf-8')
    print(json.dumps({'summary': str(out), 'summary_sha256': digest(out),
                      'ledger_sha256': digest(ledger), 'canonical_actual_mask_count': len(all_masks),
                      'rows': rows, 'terminal': summary['terminal_E']}, indent=2))


if __name__ == '__main__':
    main()
