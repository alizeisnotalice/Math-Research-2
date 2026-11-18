"""Exact checks of affine sign alignment and maximal common quadratic factors."""
from collections import Counter
from itertools import combinations, product
import json
import platform


def bit(mask, x):
    return (mask >> x) & 1


def walsh(mask, n):
    N = 1 << n
    return [sum((-1)**(bit(mask, x) ^ ((x & y).bit_count() & 1))
                for x in range(N)) for y in range(N)]


def verify_odd(n):
    N = 1 << n
    masks = [s for s in range(1 << N) if s.bit_count() & 1]
    w = {s: walsh(s, n) for s in masks}
    hist = Counter()
    for f, g in product(masks, repeat=2):
        b = 0
        for x in range(N):
            if bit(f ^ g, x):
                b ^= x
        k = (g.bit_count() - f.bit_count())//2
        bad = 0
        for xi, (z, t) in enumerate(zip(w[f], w[g])):
            aligned = t * (-1)**((k + (b & xi).bit_count()) & 1)
            assert (z-aligned) % 8 == 0
            assert (z != aligned) == (abs(z) != abs(t))
            bad += z != aligned
        hist[bad] += 1
    return {"ordered_pairs": sum(hist.values()), "numbers_of_unmatched_magnitudes": dict(sorted(hist.items()))}


def span(generator):
    H = {0}
    for x in generator:
        if x not in H:
            H |= {v ^ x for v in H}
    return H


def subspaces(n):
    N = 1 << n
    groups = {frozenset([0])}
    for d in range(1, n+1):
        groups.update(frozenset(span(c)) for c in combinations(range(1, N), d))
    return sorted(groups, key=lambda h: (len(h), sorted(h)))


def quadratic_masks(n):
    N = 1 << n
    monomials = [1 << j for j in range(n)]
    monomials.extend((1 << j) | (1 << k) for j, k in combinations(range(n), 2))
    return sorted(span(sum(1 << x for x in range(N) if (x & m) == m) for m in monomials))


def normalized_shift(f, b, n):
    return sum((bit(f, x ^ b) ^ bit(f, b)) << x for x in range(1 << n))


def finite_quadratic(n):
    N = 1 << n
    X = list(range(N))
    laws = quadratic_masks(n)
    spaces = subspaces(n)
    splits = [(u, v) for u in spaces for v in spaces if len(u)*len(v) == N and u & v == {0}]
    candidates = {}
    for f in laws:
        dct = {}
        for j, (U, V) in enumerate(splits):
            if all(bit(f, u ^ v) == bit(f, u) ^ bit(f, v) for u in U for v in V):
                dct[j] = frozenset(u for u in U if bit(f, u))
        candidates[f] = dct
    exact = {}
    for f, g in product(laws, repeat=2):
        exact[f, g] = max(len(splits[j][0]).bit_length()-1
                          for j, sig in candidates[f].items() if candidates[g].get(j) == sig)
    shifts = {g: {normalized_shift(g, b, n) for b in X} for g in laws}
    def dimension(H):
        return len(H).bit_length()-1
    def pairing(f, x, y):
        return bit(f, x ^ y) ^ bit(f, x) ^ bit(f, y)
    def radical(f, H):
        return {x for x in H if all(not pairing(f, x, y) for y in H)}
    roots = {f: radical(f, X) for f in laws}
    corr = Counter()
    for f, g in product(laws, repeat=2):
        R = f ^ g
        K = roots[R]
        K0 = {x for x in K if not bit(R, x)}
        Z0 = roots[f] & K0
        initial = dimension(K0)-dimension(radical(f, K0))+dimension(Z0)
        sheared = dimension(K)-dimension(radical(f, K))+dimension(Z0)
        assert initial == exact[f, g]
        assert sheared == max(exact[f, u] for u in shifts[g])
        assert n-sheared <= 2*(n-dimension(K))+1
        corr[sheared-initial] += 1
    return {"ordered_pairs": len(exact), "complementary_direct_sum_splits": len(splits),
            "increase_in_common_dimension_after_shear": dict(sorted(corr.items()))}


def bent_correlation(n):
    N = 1 << n
    fixed_quadratics = [m for m in quadratic_masks(n) if all(not bit(m, 1 << x) for x in range(n))]
    bent = {m: walsh(m, n) for m in fixed_quadratics}
    bent = {m: w for m, w in bent.items() if all(abs(k) == 1 << (n//2) for k in w)}
    dist = Counter()
    for f, g in product(bent, repeat=2):
        dual_sum = sum(1 << i for i in range(N) if bent[f][i]*bent[g][i] < 0)
        best = max(map(abs, walsh(dual_sum,n)))
        primal = max(map(abs, walsh(f ^ g,n)))
        assert best == primal
        dist[best] += 1
    return {"ordered_pairs_of_nonsingular_polar_forms": len(bent)**2,
            "max_affine_signed_correlation_numerator": dict(sorted(dist.items())), "denominator": N}


print(json.dumps({"python": platform.python_version(),
                  "all_odd_child_pairs": {n: verify_odd(n) for n in (2,3)},
                  "all_quadratic_child_pairs_modulo_global_signs": {n: finite_quadratic(n) for n in (2,3)},
                  "four_variable_bent_pairs": bent_correlation(4)}, indent=2))
