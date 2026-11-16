#!/usr/bin/env python3
"""Exact elementary checks used in the G(x)G(x^L) row proof.

No searcher or floating point arithmetic is used.  The script checks the
support gaps and the coefficient identities for GQ quoted in the proof.
"""
from __future__ import annotations

I = {0, 2, 3, 5, 6}

def conv_support(a, b):
    return {u + v for u in a for v in b}

def coeff_gq(d):
    # Q has possible coefficients q_0,q_2,q_3,q_5,q_6 only.
    return tuple(sorted(d - i for i in I if d - i in I))

def main():
    # These are sets of Q-subscripts occurring in [x^d]GQ.
    expected = {4: (2,), 7: (2, 5), 8: (2, 3, 5, 6), 9: (3, 6)}
    actual = {d: coeff_gq(d) for d in expected}
    assert actual == expected, (actual, expected)

    for L in range(7, 101):
        S = {i + L*j for i in I for j in I}
        # The base calculation needs these three globally zero coefficients.
        assert all(k not in S for k in (7, 8, 9))
        # Every row residue outside I is zero; the carry locations selected in
        # the extraction lemma have residue 1 when L=7 or L=8.
        for r in range(7):
            assert r * L + 4 not in S
            if L == 7:
                assert r * L + 8 not in S
            elif L == 8:
                assert r * L + 7 not in S and r * L + 9 not in S
            else:
                assert r * L + 7 not in S and r * L + 8 not in S

    print('exact coefficient identities:', actual)
    print('support/carry zero checks: L=7..100 passed')

if __name__ == '__main__':
    main()
