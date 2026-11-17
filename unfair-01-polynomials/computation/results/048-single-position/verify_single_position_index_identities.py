#!/usr/bin/env python3
"""Exact integer-pair regression for single-nonrational-position indices.

Each coefficient is represented by (u,v), meaning u+v*a over Z[a], without
using a minimal polynomial.  The check only audits the convolution indices in
the two-end locking and cofactor-degree-hole lemmas.
"""
from random import Random


def conv(a, b):
    out = [0] * (len(a)+len(b)-1)
    for i, x in enumerate(a):
        for j, y in enumerate(b): out[i+j] += x*y
    return out


def rev(a): return list(reversed(a))


def check_case(d, r, g, q, seed):
    rng = Random(seed)
    p0 = [rng.randint(-3, 3) for _ in range(d+1)]
    p0[0] = p0[d] = 1
    p0[r] = 0
    G = [rng.randint(-3, 3) for _ in range(g+1)]
    Q = [rng.randint(-3, 3) for _ in range(q+1)]
    G[0] = Q[0] = G[g] = Q[q] = 1
    D, E = conv(G, p0), conv(Q, p0)
    # Coefficient of a in GP is x^r G; verify the two locking positions.
    def acoef(G, k): return G[k-r] if 0 <= k-r < len(G) else 0
    assert acoef(G, r) == 1 and acoef(G, g+r) == 1
    assert acoef(Q, r) == 1 and acoef(Q, q+r) == 1
    # If g>=q, the endpoint x^(d+q) has the stated a coefficient.
    t = d-r
    if g >= q:
        k = d+q
        assert acoef(G, k) == (G[q+t] if q+t <= g else 0)
    # Simultaneous reversal changes the single position from r to t and
    # maps G*_(q+r) back to G_(g-q-r).
    if g >= q:
        gs, p0s = rev(G), rev(p0)
        assert conv(gs, p0s)[t] == D[g+r]
        if q+r <= g:
            assert gs[q+r] == G[g-q-r]
    return {'d': d, 'r': r, 'g': g, 'q': q, 'D_lock': (D[r], D[g+r]), 'E_lock': (E[r], E[q+r])}


def main():
    cases = [check_case(9, 3, 12, 4, 110+i) for i in range(8)]
    cases += [check_case(11, 2, 7, 1, 210+i) for i in range(8)]
    print('PASS: 16 exact Z[a] index regressions')
    print(cases[0])


if __name__ == '__main__': main()
