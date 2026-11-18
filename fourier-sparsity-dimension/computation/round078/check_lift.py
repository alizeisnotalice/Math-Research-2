"""Integer replay of the shared two-label tree and its forced relation module."""
from collections import Counter, defaultdict
from fractions import Fraction
from itertools import combinations
from math import ceil
import json
import platform


def rank(vectors):
    basis = {}
    for x in vectors:
        while x and x.bit_length() in basis:
            x ^= basis[x.bit_length()]
        if x:
            basis[x.bit_length()] = x
    return len(basis)


def walsh(values):
    v = list(values)
    t = 1
    while t < len(v):
        for i in range(0, len(v), 2*t):
            for a in range(i, i+t):
                lo, hi = v[a], v[a+t]
                v[a], v[a+t] = lo+hi, lo-hi
        t *= 2
    return v


def orbit_density(orbits):
    best = Fraction(0)
    for k in range(1, 1+len(orbits)):
        for subset in combinations(orbits, k):
            union = set.union(*subset)
            best = max(best, Fraction(len(union), rank(union)))
    return best


def direct_checks(spec, denominator, j):
    L = 1 << j
    d = 3*L+1
    truth = []
    for u in range(1 << d):
        node = 1
        h = 0
        for _ in range(j):
            side = (u >> (node-1)) & 1
            h ^= side
            node = 2*node + side
        leaf = node-L
        v = (u >> (L-1)) & 1
        w = (u >> L) & 1
        z = (u >> (L+1+2*leaf)) & 3
        truth.append(w ^ (v*h) ^ int(v == 0 and z == 0))
    ww = walsh([1-2*x for x in truth])
    assert {a: c*denominator for a, c in enumerate(ww) if c} == {
        a: c*(1 << d) for a, c in spec.items()}
    anf = truth[:]
    for i in range(d):
        bit = 1 << i
        for a in range(1 << d):
            if a & bit:
                anf[a] ^= anf[a ^ bit]
    degree = max(a.bit_count() for a, c in enumerate(anf) if c)
    assert degree == j+3
    return len(truth)


def rec(j):
    L = 1 << j
    denominator = 4*L
    v, w = 1 << (L-1), 1 << L
    addrmask = v-1
    spec = defaultdict(int)
    unique_labels = {}
    relations = []
    forced = 0
    for leaf in range(L):
        node = 1
        signs = [(0, 1)]
        for k in range(j-1, -1, -1):
            side = (leaf >> k) & 1
            bit = 1 << (node-1)
            signs += [(a | bit, s*(1-2*side)) for a, s in signs]
            node = 2*node + side
        h = leaf.bit_count() % 2
        row = {}
        for t in range(8):
            label = w | ((t & 1)*v) | ((t >> 1) << (L+1+2*leaf))
            if label not in unique_labels:
                unique_labels[label] = len(unique_labels)
            c = 3 if t == h else -1
            row[label] = c
            for a, s in signs:
                spec[a | label] += c*s
        groups = defaultdict(list)
        for a, b in combinations(row, 2):
            groups[a^b].append((row[a]*row[b], (1 << unique_labels[a]) |
                                 (1 << unique_labels[b])))
        for pairs in groups.values():
            assert sorted(c for c, _ in pairs) == [-3, 1, 1, 1]
            forced += 1
            relations.extend(q ^ pairs[0][1] for _, q in pairs[1:])
    spec = {f: c for f, c in spec.items() if c}
    bands = defaultdict(set)
    for f in spec:
        bands[f & ~addrmask].add(f & addrmask)
    K = rank(x ^ min(b) for b in bands.values() for x in b)
    R = rank(relations)
    d = rank(spec)
    s = len(spec)
    assert K == L-1 and R == 4*L and d == K+len(bands)-R == 3*L+1
    assert s == (8 if j == 0 else 2+13*L*L//2)
    assert len(bands) == 2+6*L
    assert sum(c*c for c in spec.values()) == denominator**2
    convolution = defaultdict(int)
    for (a, b), (c, e) in combinations(spec.items(), 2):
        convolution[a ^ c] += b * e
    assert all(t == 0 for t in convolution.values())
    delta = (1-max(Fraction(abs(c), denominator) for c in spec.values()))/2
    assert delta == (Fraction(1, 8) if not j else
                     (1-max(Fraction(1, 4), Fraction(1, L)))/2)
    private = {f for f in spec if f >> (L+1)}
    orbits = [set(spec)]
    if j:
        c0 = {v | w, w}
        c1 = set(spec) - c0 - private
        orbits = [c0, c1, private]
        assert list(map(rank, orbits)) == [2, L, d]
        assert rank(c0 | c1) == L+1
    density = orbit_density(orbits)
    assert density == Fraction(s, d)
    return dict(j=j, leaves=L, s=s, d=d, c=ceil(density), density=str(density),
                delta_star=str(delta), degree=j+3, K=K, relation_rank=R,
                intersection_rank=R, data_labels=len(bands),
                forced_matchings=forced, denominator=denominator,
                signed_numerators=dict(sorted(Counter(spec.values()).items())),
                nonzero_differences=len(convolution),
                truth_entries_checked=direct_checks(spec, denominator, j) if j<=2 else 0)


def fake(m):
    # A necessary condition: each v-slice on the independent e_i is a
    # Boolean linear form, so at most one coefficient can be nonzero.
    # Both f_i^+=a_i+b_i and f_i^-=a_i-b_i cannot vanish for even one
    # nonzero pair. Test every pair of allowed slice supports.
    all_bits = (1 << m)-1
    singleton_or_empty = [0] + [1 << i for i in range(m)]
    return any((a|b) == all_bits for a in singleton_or_empty for b in singleton_or_empty)


if __name__ == '__main__':
    rows = [rec(j) for j in range(5)]
    rejects = [m for m in range(3, 65) if not fake(m)]
    assert len(rejects) == 62
    print(json.dumps(dict(python=platform.python_version(), tolerance=0,
                          generations=rows, rejected_fake_m=rejects), indent=2))
