"""Exact finite checks for the rank-at-most-four extension."""
import itertools
import json
import random
import sys
from fractions import Fraction

from verify_phases import affstats, alternating_rows, rank, single_flip

SEED = 8204


def evaluate_pencil(rows, data_dimension):
    ranks = [rank(matrix) for matrix in rows]
    assert max(ranks) <= 4
    total = sum(1 << r for r in ranks)
    address_rank = rank(vector for matrix in rows for vector in matrix)
    assert (1 << data_dimension) * (data_dimension + 2 + address_rank) ** 2 <= (
        (data_dimension + 2) ** 2 * total
    )
    return max(ranks)


def exhaustive_four_address():
    rows = [alternating_rows(4, mask) for mask in range(1 << 6)]
    checked = rank_four = 0
    for a, b, c in itertools.product(range(1 << 6), repeat=3):
        rank_max = evaluate_pencil((rows[a], rows[b], rows[c], rows[a ^ b ^ c]), 2)
        checked += 1
        rank_four += rank_max == 4
    return checked, rank_four


def wedge(size, x, y):
    return tuple(
        (y if (x >> bit) & 1 else 0) ^ (x if (y >> bit) & 1 else 0)
        for bit in range(size)
    )


def combination(linear_terms, beta):
    value = linear_terms[0]
    for bit, term in enumerate(linear_terms[1:]):
        if (beta >> bit) & 1:
            value ^= term
    return value


def random_structured(rng):
    total = rank_four = 0
    for size in range(4, 13):
        for data_dimension in range(2, 8):
            for _ in range(40):
                fixed = [rng.getrandbits(size), rng.getrandbits(size)]
                terms = [[rng.getrandbits(size) for _ in range(data_dimension + 1)] for _ in fixed]
                matrices = []
                for beta in range(1 << data_dimension):
                    rows0, rows1 = (
                        wedge(size, base, combination(coefficients, beta))
                        for base, coefficients in zip(fixed, terms)
                    )
                    matrices.append(tuple(x ^ y for x, y in zip(rows0, rows1)))
                rank_four += evaluate_pencil(matrices, data_dimension) == 4
                total += 1
    return total, rank_four


def verify_vertices():
    proofs = [
        (Fraction(1, 4), Fraction(1, 4), Fraction(4, 3)),
        (Fraction(1, 4), Fraction(3, 4), Fraction(23, 14)),
        (Fraction(1, 2), Fraction(0), Fraction(19, 10)),
        (Fraction(1), Fraction(0), Fraction(3)),
    ]
    checks = 0
    for p, q, lower in proofs:
        assert (1 + lower) ** 2 <= 1 + 3 * q + 15 * p
        w = q + 3 * p
        for v in range(2, 49):
            assert (v + 2) * lower - v * w - 4 + v * (3 - w) / ((1 << v) - 1) >= 0
            checks += 1
    return checks


def mixed_walsh():
    g = single_flip(2)
    values = []
    for point in range(64):
        address = point & 15
        data = point >> 4
        shift = int(address & 3 == 3) | (int(address & 12 == 12) << 1)
        values.append(g[data ^ shift])
    stats = affstats(values)
    assert stats == (25, 7)
    return {"s_parent": stats[0], "R_parent": stats[1], "ranks": [0, 2, 2, 4]}


if __name__ == "__main__":
    rng = random.Random(SEED)
    exhaustive, rank_four = exhaustive_four_address()
    sampled, sampled_four = random_structured(rng)
    print(json.dumps({
        "python": sys.version.split()[0], "seed": SEED,
        "all_four_address_two_data_pencils": exhaustive,
        "of_which_max_rank_four": rank_four,
        "sampled_structured_pencils": sampled,
        "sampled_of_which_max_rank_four": sampled_four,
        "exact_vertex_inequalities": verify_vertices(),
        "mixed_rank_exact_Walsh": mixed_walsh(),
    }, indent=2))
