"""Exact magnitude classes and genuinely split zero-sum modules on four bits."""
from collections import Counter, defaultdict
from fractions import Fraction
from itertools import combinations
from math import ceil
import json
import platform

from check_lift import rank, walsh


def audit_odd4():
    groups = defaultdict(list)
    for word in range(1 << 16):
        if not word.bit_count() % 2:
            continue
        coeffs = walsh([1-2*((word >> x) & 1) for x in range(16)])
        mag = tuple(abs(a) for a in coeffs)
        phase = sum((a < 0) << xi for xi, a in enumerate(coeffs))
        groups[mag].append(phase)
    affine = {sum((((a & x).bit_count() % 2) ^ t) << x for x in range(16))
              for a in range(16) for t in range(2)}
    assert len(affine) == 32
    split_spectra = split_matchings = 0
    for mags, phases in groups.items():
        base = phases[0]
        assert {base ^ p for p in phases} == affine and len(phases) == 32
        coeffs = [a*(1-2*((base >> i) & 1)) for i, a in enumerate(mags)]
        relations = []
        saw_split = False
        for d in range(1, 16):
            pairs = [(i, i ^ d) for i in range(16) if i < i ^ d]
            sums = [coeffs[i]*coeffs[h] for i, h in pairs]
            assert sum(sums) == 0
            masks = [tuple(range(len(pairs)))]
            for size in range(2, 5):
                first = next((ix for ix in combinations(range(len(pairs)), size)
                              if sum(sums[k] for k in ix) == 0), None)
                if first:
                    masks = [first, sorted(set(range(8)) - set(first))]
                    saw_split = True
                    split_matchings += 1
                    break
            for part in masks:
                mask = [(1 << pairs[q][0]) | (1 << pairs[q][1]) for q in part]
                relations.extend(mask[0] ^ z for z in mask[1:])
        assert rank(relations) == 11
        split_spectra += saw_split
    return dict(weight_odd_truths=sum(map(len, groups.values())),
                magnitude_classes=len(groups), affine_signs_per_class=len(affine),
                classes_with_split=split_spectra,
                genuinely_split_matchings=split_matchings,
                relation_rank_after_greedy_split=11)


def sharp(j):
    L = 1 << j
    w = 1 << (L-1)
    coeff = defaultdict(int)
    denom = 4*L
    orbits = []
    leaves = []
    for i in range(L):
        path = [(0, 1)]
        node = 1
        for level in reversed(range(j)):
            side = (i >> level) & 1
            bit = 1 << (node-1)
            path += [(a ^ bit, b*(1-2*side)) for a, b in path]
            node = 2*node + side
        private = set()
        for u in range(8):
            heavy = int(i == L-1)
            c = 3 if u == heavy else -1
            data = w | (u << (L+3*i))
            for a, b in path:
                coeff[a | data] += c*b
                if u:
                    private.add(a | data)
        leaves.append(private)
    coeffs = {x: c for x, c in coeff.items() if c}
    assert sum(c*c for c in coeffs.values()) == denom**2
    cor = defaultdict(int)
    for (x, c), (y, b) in combinations(coeffs.items(), 2):
        cor[x ^ y] += c*b
    assert not any(cor.values())
    common = set(coeffs)-set.union(*leaves)
    s = len(coeffs)
    d = rank(coeffs)
    assert s == 7*L*L+L and d == 4*L and len(common) == L
    delta = (1-max(Fraction(abs(c),denom) for c in coeffs.values()))/2
    assert delta == (Fraction(1,8) if L == 1 else
                     (1-max(abs(Fraction(3,4)-Fraction(1,L)),Fraction(1,L)))/2)
    maxden = Fraction(s,d)
    # Sweep every subset of leaves to replay the density majorant.
    checked = 0
    if L <= 8:
        for m in range(1,L):
            for subset in combinations(leaves, m):
                pts = set.union(*subset)
                assert Fraction(len(pts)+L,rank(pts)) < maxden
                checked += 1
    truth_checked = 0
    if j <= 2:
        truth = []
        for x in range(1 << d):
            node = 1
            for step in range(j):
                node = 2*node + ((x >> (node-1)) & 1)
            i = node-L
            z = (x >> (L+3*i)) & 7
            truth.append(int(bool(x & w)) ^ int(not z) ^ ((z & 1) if i == L-1 else 0))
        spectra = walsh([1-2*y for y in truth])
        assert {xi:c*denom for xi,c in enumerate(spectra) if c} == {
            xi:c*(1 << d) for xi,c in coeffs.items()}
        for k in range(d):
            for x in range(1 << d):
                if x & (1 << k):
                    truth[x] ^= truth[x ^ (1 << k)]
        assert max(x.bit_count() for x, c in enumerate(truth) if c) == j+3
        truth_checked = len(truth)
    return dict(j=j,s=s,d=d,c=ceil(maxden),delta_star=str(delta),degree=j+3,
                K=L-1,intersection_rank=4*L,data_labels=1+7*L,denominator=denom,
                signed_numerators=dict(sorted(Counter(coeffs.values()).items())),
                nonzero_differences=len(cor),subset_density_checks=checked,
                truth_entries_checked=truth_checked)


if __name__ == '__main__':
    print(json.dumps(dict(python=platform.python_version(),tolerance=0,
                          odd_four_bit=audit_odd4(),
                          sharp_generations=[sharp(j) for j in range(5)]),indent=2))
