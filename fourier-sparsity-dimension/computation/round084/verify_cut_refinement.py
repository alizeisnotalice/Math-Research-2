"""Exact checks of the cost cases and two genuinely different rank-two parents."""
import itertools
import json
import sys
from fractions import Fraction

from verify_nonfull import walsh


def cut_cost(rank, projected, outside):
    scale = 1 << (rank // 2)
    return Fraction(scale) if outside else scale * (1 - Fraction(1, 1 << projected))


def cost_cases():
    total = 0
    special_four = 0
    options = [
        (rank, projected, outside, cut_cost(rank, projected, outside))
        for rank in range(0, 9, 2)
        for projected in range(min(4, rank) + 1)
        for outside in (0, 1)
    ]
    for n in range(1, 5):
        for pattern in itertools.combinations_with_replacement(options, n):
            paid = sum(x[3] for x in pattern) + max(x[3] for x in pattern)
            outside = sum(x[2] for x in pattern)
            for t in range(1 + min(4, sum(x[1] for x in pattern))):
                if outside + t and paid < outside + t + 1:
                    assert t == 4
                    assert [(r, p) for r, p, _, _ in pattern if p] == [(2, 2), (2, 2)]
                    special_four += 1
            total += 1
    for bits in range(3):
        for rank in range(2 * bits, 13, 2):
            for projected in range(2 * bits, rank + 1):
                assert cut_cost(rank, projected, False) >= Fraction(3 * bits, 2)
                assert cut_cost(rank, projected, True) >= 1 + bits
    return total, special_four


def points(mask):
    while mask:
        first = mask & -mask
        yield first.bit_length() - 1
        mask ^= first


def subspaces(dimension):
    full = (1 << (1 << dimension)) - 1
    current = {1}
    for rank in range(dimension + 1):
        following = set()
        for space in current:
            yield space, rank
            support = tuple(points(space))
            outside = full ^ space
            while outside:
                sample = (outside & -outside).bit_length() - 1
                coset = sum(1 << (sample ^ x) for x in support)
                following.add(space | coset)
                outside ^= coset
        current = following


def check_l(values, address_bits):
    support = {i | len(values): abs(a) for i, a in enumerate(walsh(values)) if a}
    total = sum(support.values())
    rank = len(values).bit_length()
    author_cut_values = {}
    maximum = Fraction(0)
    tested = 0
    within_address = (1 << (1 << address_bits)) - 1
    for subspace, cut_rank in subspaces(rank):
        if cut_rank == rank:
            continue
        inside = sum(weight for label, weight in support.items() if (subspace >> label) & 1)
        value = Fraction((rank - cut_rank + 1) * len(values), total - inside)
        maximum = max(maximum, value)
        j = (subspace & within_address).bit_count().bit_length() - 1
        t = address_bits - j
        author_cut_values[t] = max(author_cut_values.get(t, Fraction(0)), value)
        tested += 1
    return {"s":len(support), "R":rank, "L":str(maximum), "tested_cuts":tested,
            "max_L_by_t": {t:str(x) for t,x in sorted(author_cut_values.items())}}


def quadratic(point, pairs):
    return sum(((point >> a) & 1) * ((point >> b) & 1) for a,b in pairs) & 1


def parent(constant, shift0, shift1):
    result = []
    for point in range(64):
        u = point & 15
        x = point >> 4
        x ^= quadratic(u, shift0) | (quadratic(u, shift1) << 1)
        result.append((-1) ** (int(x == 3) ^ quadratic(u, constant)))
    return tuple(result)


if __name__ == "__main__":
    base = check_l((1, 1, 1, -1), 0)
    common = check_l(parent(((0, 1),), ((0, 2),), ((0, 3),)), 4)
    disjoint = check_l(parent(((0, 1),), ((0, 1), (2, 3)), ((0, 2),)), 4)
    high = check_l(tuple((-1)**quadratic(x, ((0,1),(2,3),(4,5))) for x in range(64)), 6)
    mixed = check_l(parent((), ((0, 1),), ((2, 3),)), 4)
    assert base["L"] == common["L"] == disjoint["L"] == "2"
    assert mixed["s"] == 25 and mixed["L"] == "16/9"
    patterns, special_four = cost_cases()
    assert high["L"] == "1"
    print(json.dumps({"python":sys.version.split()[0], "cost_patterns":patterns,
                      "special_t4_cost_patterns":special_four,
                      "mixed_rank_parent":mixed,
                      "seed":base, "common_rank_two":common, "noncommon_rank_two":disjoint,
                      "rank_six_character_twist":high}, indent=2))
