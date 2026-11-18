"""Integer Walsh verification of the contaminated matching construction."""
from collections import Counter, defaultdict
from fractions import Fraction
from itertools import product
from random import Random
import json
import platform


def basis(rows):
    result = {}
    for v, bit in rows:
        while v:
            i = v.bit_length()
            if i not in result:
                result[i] = (v, bit)
                break
            a, b = result[i]
            v ^= a
            bit ^= b
    return tuple(result[i] for i in sorted(result, reverse=True))


def reduce_vector(v, pivots):
    bit = 0
    for a, b in pivots:
        if (v >> (a.bit_length() - 1)) & 1:
            v ^= a
            bit ^= b
    return v, bit


def rank(vectors):
    return len(basis((v, 0) for v in vectors))


def fold(source, pivots):
    sums = defaultdict(int)
    for v, c in source.items():
        p, b = reduce_vector(v, pivots)
        sums[p] += (1 - 2 * b) * c
    return {v: c for v, c in sums.items() if c}


def signed_span(rows):
    result = [(0, 0)]
    for v, a in rows:
        result += [(v ^ w, a ^ b) for w, b in result]
    return result


def fwht(values):
    data = values.copy()
    step = 1
    while step < len(data):
        for i in range(0, len(data), 2 * step):
            for j in range(i, i + step):
                x, y = data[j], data[j + step]
                data[j], data[j + step] = x + y, x - y
        step *= 2
    return data


def main():
    seed = 85029
    rng = Random(seed)
    edges = [(0, 1), (2, 4), (3, 7), (5, 6)]
    offsets = [3, 6, 9, 12]
    widths = [2, 2, 2, 3]
    masks = [((1 << (w + 1)) - 1) << o for o, w in zip(offsets, widths)]
    guards = [1 << (o + w) for o, w in zip(offsets, widths)]
    seed_masks = [m ^ g for m, g in zip(masks, guards)]
    vertices = {v: i for i, edge in enumerate(edges) for v in edge}
    augment = 1 << 16
    truth = []
    for x in range(augment):
        i = vertices[x & 7]
        datum = (x & seed_masks[i]) >> offsets[i]
        negative = (datum == (0 if i == 3 else 3)) ^ bool(x & guards[i])
        truth.append(1 - 2 * negative)
    full = fwht(truth)
    S = {v: c for v, c in enumerate(full) if c}
    n, N, r, s, b, a = 3, 4, 3, 4, 8, 4
    D = rank(augment | v for v in S) + 1
    assert D == 3*r+a+n+1 and len(S) == N*(3*s+b)
    off = {v: val for v, val in S.items() if v & seed_masks[-1]}
    assert len(off) == N*(b-1) and set(map(abs, off.values())) == {2*augment//(N*b)}
    assert rank(augment | v for v in S if v not in off) + 3 == D-1
    l1 = sum(Fraction(abs(c), augment) for c in off.values())
    energy = sum(Fraction(c, augment)**2 for c in off.values())
    assert l1 == 2*(1-Fraction(1,b)) and l1*l1 == len(off)*energy
    parts = [{v: c for v, c in S.items() if v & g} for g in guards]
    labels = sorted({v & ~7 for v in S})
    elementary = sorted({(v ^ w) | j for v in labels for w in labels for j in range(8)}-{0})
    tally = Counter()
    tested = 0
    for h in range(1, 5):
        for _ in range(150):
            H = rng.sample(elementary, h)
            if rank(H) != h:
                continue
            for bits in product((0,1), repeat=h):
                rows = list(zip(H,bits))
                allrows = signed_span(rows)
                final = basis(rows)
                fused = fold(S, final)
                newD = rank(augment | v for v in fused) + 1
                pure = [(v,z) for v,z in allrows if v and not v & ~7]
                if pure:
                    C = {x for x in range(8) if all((v & x).bit_count()%2 == z for v,z in pure)}
                    f = sum(x in C and y in C for x,y in edges)
                    t = len(C)-2*f
                    assert 2*len(fused) >= len(C)*s*(f+2*t-4*h-2)
                    assert newD <= (f+t)*r+int(any(x in C for x in edges[-1]))*(a-r)+n+1
                    tally['pure'] += 1
                    continue
                local = [[(v,z) for v,z in allrows if v & mask and not v & ~7 & ~mask] for mask in masks]
                # With data bits above address, elimination leaves private projections canonical.
                local_basis = basis(y for group in local for y in group)
                local_parts = [fold(part,local_basis) for part in parts]
                local_ranks = [rank(v & ~7 for v in part) for part in local_parts]
                post_labels = [{reduce_vector(v,final)[0]&~7 for v in part} for part in local_parts]
                c = int(any(post_labels[i] & post_labels[j] for i in range(4) for j in range(i)))
                losses = [r-ri for ri in local_ranks[:3]]
                X, Q = sum(losses), sum(mi*mi for mi in losses)
                u = int(bool(local[-1]))
                assert local_ranks[-1] <= a-u
                assert newD <= D-X-u-c
                assert (len(S)-len(fused))*r*r <= N*(s*(2*X*r-Q)+(b*u+4*(1<<h)*c)*r*r)
                tally[f'no_pure_rogue{u}_collision{c}'] += 1
                tested += 1
    return dict(python=platform.python_version(),seed=seed,tolerance=0,truth_values=augment,
                parent_rank_plus_one=D,sparsity=len(S),cut_codimension=3,outside_count=len(off),
                outside_l1=str(l1),outside_energy=str(energy),query_answers=dict(tally),no_pure=tested)


def two_special_leaves():
    edges = [(0,1), (2,4), (3,7), (5,6)]
    offsets = [3,6,9,13]
    widths = [2,2,3,3]
    size = 1 << 17
    parent = []
    for x in range(size):
        i = next(i for i, edge in enumerate(edges) if x & 7 in edge)
        o, w = offsets[i], widths[i]
        datum = (x >> o) & ((1 << w) - 1)
        sign = (datum == (3 if i < 2 else 0)) ^ bool(x & (1 << (o+w)))
        parent.append(1-2*sign)
    spectrum = {v: c for v,c in enumerate(fwht(parent)) if c}
    full_rank = rank(size | v for v in spectrum)
    masks = [7 << offsets[2], 7 << offsets[3]]
    supports = [{v for v in spectrum if v & mask} for mask in masks]
    assert not (supports[0] & supports[1])
    rows = []
    for mask, outside in [(masks[0],supports[0]),(masks[1],supports[1]),
                          (masks[0]|masks[1],supports[0]|supports[1])]:
        loss = full_rank-rank(size | v for v in spectrum if not (v & mask))
        energy = sum(Fraction(spectrum[v],size)**2 for v in outside)
        l1 = sum(Fraction(abs(spectrum[v]),size) for v in outside)
        assert loss == 3*(len(outside)//28)
        assert l1 == Fraction(7,4)*(len(outside)//28)
        assert len(outside)*energy == l1*l1
        rows.append(dict(codimension=loss,points=len(outside),l1=str(l1),energy=str(energy)))
    return dict(truth_values=size,parent_rank_plus_one=full_rank+1,sparsity=len(spectrum),cuts=rows)


if __name__ == '__main__':
    result = main()
    result['two_special_leaves'] = two_special_leaves()
    print(json.dumps(result, indent=2))
