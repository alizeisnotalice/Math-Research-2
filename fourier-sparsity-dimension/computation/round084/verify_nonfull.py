"""Exact Walsh and affine-basis checks for the nonfull-spectrum bounds."""
from collections import Counter
from fractions import Fraction
from functools import cache
import itertools
import json
import math
import random
import sys

SEED = 8402


def walsh(values):
    values = list(values)
    width = 1
    while width < len(values):
        for base in range(0, len(values), 2 * width):
            for index in range(base, base + width):
                left, right = values[index], values[index + width]
                values[index], values[index + width] = left + right, left - right
        width *= 2
    return values


def rank(vectors):
    pivots = {}
    for vector in vectors:
        while vector:
            pivot = vector.bit_length()
            if pivot not in pivots:
                pivots[pivot] = vector
                break
            vector ^= pivots[pivot]
    return len(pivots)


def spectrum(values):
    return tuple(index for index, value in enumerate(walsh(values)) if value)


def stats(values):
    support = spectrum(values)
    return len(support), rank(label | len(values) for label in support)


def eigenphase(address, even_rank):
    return sum(
        ((address >> bit) & 3) == 3 for bit in range(0, even_rank, 2)
    ) & 1


def four_leaf_truth(point, even_rank):
    root = (point & 1) ^ eigenphase(point >> 7, even_rank)
    second = (point >> (1 + root)) & 1
    terminal = 3 + 2 * root + second
    return (-1) ** ((point >> terminal) & 1)


def gamma_uniform_fixed_base(quotient, homogeneous_marker):
    columns = tuple(x | homogeneous_marker for x in quotient)
    base_size = rank(columns)
    bases = [base for base in itertools.combinations(columns, base_size) if rank(base) == base_size]
    incidences = Counter(x for base in bases for x in base)
    return max(Fraction(incidences[x] * len(columns), len(bases) * base_size) for x in columns), len(bases)


def check_true_cylinder():
    seed_truth = tuple(four_leaf_truth(x, 0) for x in range(128))
    seed_support = spectrum(seed_truth)
    assert stats(seed_truth) == (16, 7)
    assert all(beta ^ 1 in seed_support for beta in seed_support)
    quotient = tuple(beta for beta in seed_support if beta & 1 == 0)
    cap, basis_count = gamma_uniform_fixed_base(quotient, 1 << 7)
    assert cap == 1
    checked = []
    for even_rank in (2, 4, 6, 8):
        values = tuple(four_leaf_truth(point, even_rank) for point in range(1 << (7 + even_rank)))
        sparsity, homogeneous_rank = stats(values)
        assert (sparsity, homogeneous_rank) == (8 * (1 + (1 << even_rank)), 7 + even_rank)
        assert (homogeneous_rank + 1) ** 2 <= 4 * sparsity
        checked.append(dict(polar_rank=even_rank, s=sparsity, R=homogeneous_rank))
    return {"quotient_bases": basis_count, "marginal_cap_multiplier": str(cap), "parents": checked}


def selector_truth(point, n, k, r):
    leaf_count = 1 << k
    address = point & (leaf_count - 1)
    seed_x = (point >> k) & ((1 << n) - 1)
    guard = (point >> (k + n + address)) & 1
    offset = k + n + leaf_count
    if address == 0:
        return (-1) ** (guard ^ int(seed_x == 0))
    phase = eigenphase(point >> (offset + (address - 1) * r), r)
    return (-1) ** (guard ^ phase)


def check_log_examples():
    checked = []
    for n, k, r in ((2, 1, 2), (4, 1, 4), (4, 2, 2), (6, 2, 2)):
        slots = 1 << k
        base_size = n + k + slots
        parent_size = base_size + (slots - 1) * r
        base = stats(tuple(selector_truth(x, n, k, 0) for x in range(1 << base_size)))
        parent = stats(tuple(selector_truth(x, n, k, r) for x in range(1 << parent_size)))
        assert base == (slots * ((1 << n) + slots - 1), base_size)
        assert parent == (slots * ((1 << n) + (slots - 1) * (1 << r)), parent_size)
        checked.append(dict(n=n, k=k, polar_rank=r, seed_s_R=base, parent_s_R=parent))
    scaling = []
    for n in (64, 128, 512, 2048, 8192):
        k = (n**2).bit_length() - 1
        slots = 1 << k
        r = 2 * ((n - 2 * k) // 2)
        numerator = n + k + slots + (slots - 1) * r + 1
        denominator = n + k + slots + 1
        density_factor = Fraction((1 << n) + (slots - 1) * (1 << r), (1 << n) + slots - 1)
        log_seed = n + k + math.log2(1 + float(Fraction(slots - 1, 1 << n)))
        normalized = numerator / denominator / math.sqrt(float(density_factor)) / log_seed
        scaling.append(dict(n=n, k=k, polar_rank=r, gain_over_log_seed=normalized))
    return {"exact_Walsh": checked, "asymptotic_formula_evaluations": scaling}


@cache
def positions(size):
    return tuple(itertools.combinations(range(size), 2))


def alternating_rows(size, bits):
    vectors = [0] * size
    for index, (i, j) in enumerate(positions(size)):
        if (bits >> index) & 1:
            vectors[i] ^= 1 << j
            vectors[j] ^= 1 << i
    return tuple(vectors)


def masks_affine_pencil(generators, data_labels):
    out = []
    for label in data_labels:
        mask = generators[0]
        for index, generator in enumerate(generators[1:]):
            if (label >> index) & 1:
                mask ^= generator
        out.append(mask)
    return out


def check_matched_examples():
    rng = random.Random(SEED)
    checked = rank_counts = 0
    expected_caps = []
    geometries = (((0, 1, 2), (3,), 3), ((0, 1, 2, 3, 4), (5, 6), 4))
    for quotient, unmatched, v in geometries:
        cap, count = gamma_uniform_fixed_base(quotient, 1 << v)
        assert cap == (1 if v == 3 else Fraction(5, 4))
        expected_caps.append(dict(v=v, cap=str(cap), checked_bases=count))
        eta = 1 << (v - 1)
        paired = tuple(x ^ t for x in quotient for t in (0, eta))
        all_labels = paired + unmatched
        assert rank(x | (1 << v) for x in paired) == v + 1
        tau = Fraction(len(paired), len(all_labels))
        for size in range(2, 13):
            for _ in range(20):
                generators = [rng.getrandbits(len(positions(size))) for _ in range(v + 1)]
                d = rank(alternating_rows(size, generators[-1]))
                if not d:
                    continue
                masks = masks_affine_pencil(generators, all_labels)
                matrices = [alternating_rows(size, mask) for mask in masks]
                ranks = [rank(matrix) for matrix in matrices]
                rank_counts = max(rank_counts, max(ranks))
                a = rank(vector for matrix in matrices for vector in matrix)
                B = Fraction(sum(ranks[:len(paired)]), len(paired))
                assert a <= d + cap * v * (B - Fraction(d, 2))
                s = len(all_labels)
                s_new = sum(1 << r for r in ranks)
                theta2 = Fraction((v + 2 + a) ** 2 * s, (v + 2) ** 2 * s_new)
                constant2 = (1 + (cap - 1) * Fraction(v, v + 2)) ** 2 / tau
                assert theta2 <= constant2
                R = v + 1
                assert (1 << a) * R**R <= (R + s_new - s) ** R
                checked += 1
    return {"examples": checked, "max_polar_rank": rank_counts, "caps": expected_caps}


if __name__ == "__main__":
    print(json.dumps({
        "python": sys.version.split()[0], "seed": SEED,
        "cylinder": check_true_cylinder(),
        "log_growth": check_log_examples(),
        "matched_abstract_spectra": check_matched_examples(),
    }, indent=2))
