#!/usr/bin/env python3
"""Direct 3x3 verification of the stable second-birth survival identity."""
from collections import Counter, defaultdict
from fractions import Fraction
from itertools import permutations
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
N = 9


def xy(v):
    return v % 3, v // 3


ADJ = [set() for _ in range(N)]
for v in range(N):
    x, y = xy(v)
    for dx, dy in ((1, 0), (-1, 0), (0, 1), (0, -1)):
        xx, yy = x + dx, y + dy
        if 0 <= xx < 3 and 0 <= yy < 3:
            ADJ[v].add(yy * 3 + xx)


def births(word):
    exposed, result = set(), []
    for rank, v in enumerate(word, start=1):
        if not (ADJ[v] & exposed):
            result.append((rank, v))
        exposed.add(v)
    return result


def distance(u, v):
    x, y = xy(u)
    xx, yy = xy(v)
    return abs(x - xx) + abs(y - yy)


def euclidean_square(u, v):
    x, y = xy(u)
    xx, yy = xy(v)
    return (x - xx) ** 2 + (y - yy) ** 2


def doubled_center(v):
    x, y = xy(v)
    return 2 * x - 2, 2 * y - 2


def root_moments(x, y):
    """The F,G integrands used in the root-moment exact account."""
    tx, ty = doubled_center(x)
    ux, uy = doubled_center(y)
    return tx * tx + ty * ty + ux * ux + uy * uy, tx * ux + ty * uy


def sigma_and_L(word, tau, y):
    sigma = (word[0], y) + word[1:tau - 1] + word[tau:]
    assert births(sigma)[0:2] == [(1, word[0]), (2, y)]
    assert len(births(sigma)) == 2
    j = next(rank for rank, v in enumerate(sigma, start=1)
             if rank >= 3 and v in ADJ[y])
    L = j - 2
    assert L >= tau - 1
    return sigma, L


def main():
    direct = Counter()       # (distance, original second-birth rank)
    sigma_records = {}       # each tau=2 word -> (distance,L)
    direct_tau_sum = 0
    direct_euclidean_square = direct_F = direct_G = 0
    for word in permutations(range(N)):
        b = births(word)
        if len(b) != 2:
            continue
        (one, x), (tau, y) = b
        assert one == 1
        d = distance(x, y)
        direct[(d, tau)] += 1
        direct_tau_sum += tau
        direct_euclidean_square += euclidean_square(x, y)
        f, g = root_moments(x, y)
        direct_F += f
        direct_G += g
        sigma, L = sigma_and_L(word, tau, y)
        previous = sigma_records.get(sigma)
        if previous is None:
            sigma_records[sigma] = (d, L)
        else:
            assert previous == (d, L)

    # Independent enumeration of every tau=2 word and all of its legal slots.
    independent_sigma = {}
    clean_sigma = Counter()
    clean_weighted = Counter()
    generated = Counter()
    generated_words = set()
    for sigma in permutations(range(N)):
        b = births(sigma)
        if len(b) != 2 or b[1][0] != 2:
            continue
        x, y = b[0][1], b[1][1]
        d = distance(x, y)
        j = next(rank for rank, v in enumerate(sigma, start=1)
                 if rank >= 3 and v in ADJ[y])
        L = j - 2
        independent_sigma[sigma] = (d, L)
        eta = sigma[:1] + sigma[2:]
        # "Clean deletion" means that removing the second peak leaves exactly
        # one birth in the induced order.  This is independently useful when
        # comparing the tau=2 slice with one-peak words on G-y.
        if len(births(eta)) == 1:
            clean_sigma[d] += 1
            clean_weighted[d] += L
        for slot_after in range(1, L + 1):
            word = eta[:slot_after] + (y,) + eta[slot_after:]
            bword = births(word)
            assert len(bword) == 2 and bword[1] == (slot_after + 1, y)
            generated[(d, slot_after + 1)] += 1
            generated_words.add(word)

    # Independently enumerate a deleted vertex y and a one-birth order on
    # G-y.  Insertion of y at rank 2 must recover exactly the clean subset.
    clean_from_onepeak = Counter()
    clean_weighted_from_onepeak = Counter()
    for y in range(N):
        remaining = tuple(v for v in range(N) if v != y)
        for eta in permutations(remaining):
            if len(births(eta)) != 1 or eta[0] in ADJ[y]:
                continue
            sigma = (eta[0], y) + eta[1:]
            assert len(births(sigma)) == 2 and births(sigma)[1] == (2, y)
            j = next(rank for rank, v in enumerate(sigma, start=1)
                     if rank >= 3 and v in ADJ[y])
            d, L = distance(eta[0], y), j - 2
            clean_from_onepeak[d] += 1
            clean_weighted_from_onepeak[d] += L

    assert independent_sigma == sigma_records
    assert generated == direct
    assert sum(generated.values()) == len(generated_words)
    assert clean_from_onepeak == clean_sigma
    assert clean_weighted_from_onepeak == clean_weighted

    source = json.loads((Path("/Volumes/data/CodexTemp/T2/a022-distance-birth-exact-r020") /
                         "n3-direct.json").read_text())
    source_joint = {(int(d), int(t)): int(c)
                    for d, row in source["joint_distance_birth"].items()
                    for t, c in row.items()}
    assert dict(direct) == source_joint

    weighted_euclidean_square = weighted_F = weighted_G = 0
    for sigma, (_, L) in independent_sigma.items():
        x, y = sigma[:2]
        weighted_euclidean_square += L * euclidean_square(x, y)
        f, g = root_moments(x, y)
        weighted_F += L * f
        weighted_G += L * g
    assert weighted_euclidean_square == direct_euclidean_square
    assert weighted_F == direct_F
    assert weighted_G == direct_G
    root_source = json.loads((Path("/Volumes/data/CodexTemp/T2/a022-archive-prepare-r021") /
                              "organized/exact-root-n7/n3.json").read_text())
    assert direct_euclidean_square == int(root_source["squared_euclidean_total"])
    assert direct_F == int(root_source["firstsum"])
    assert direct_G == int(root_source["secondsum"])

    survival = {}
    for d in sorted({d for d, _ in direct}):
        for tau in range(2, N + 1):
            survival[(d, tau)] = sum(1 for dd, L in independent_sigma.values()
                                     if dd == d and L >= tau - 1)
            assert survival[(d, tau)] == direct[(d, tau)]
    sum_L = sum(L for _, L in independent_sigma.values())
    sum_L2 = sum(L * L for _, L in independent_sigma.values())
    A_tau2 = len(independent_sigma)
    A2 = sum(direct.values())
    assert A2 == sum_L
    assert A_tau2 == sum(direct[(d, 2)] for d in range(2, 5))
    mean_tau_formula = Fraction(sum_L2 + 3 * sum_L, 2 * sum_L)
    mean_tau_direct = Fraction(direct_tau_sum, A2)
    assert mean_tau_formula == mean_tau_direct

    result = {
        "n": 3, "method": "independent all-9!-word direct enumeration plus all tau=2 words and slots",
        "A2": A2, "A_tau2": A_tau2, "sum_L": sum_L, "sum_L2": sum_L2,
        "E_tau2_L": str(Fraction(sum_L, A_tau2)),
        "E2_tau_direct": str(mean_tau_direct),
        "E2_tau_from_L": str(mean_tau_formula),
        "weighted_euclidean_square_total": weighted_euclidean_square,
        "full_euclidean_square_total": direct_euclidean_square,
        "weighted_root_F_total": weighted_F,
        "full_root_F_total": direct_F,
        "weighted_root_G_total": weighted_G,
        "full_root_G_total": direct_G,
        "clean_deletion_tau2_counts_by_distance": {str(d): clean_sigma[d]
                                                    for d in sorted(clean_sigma)},
        "clean_deletion_L_weight_by_distance": {str(d): clean_weighted[d]
                                                   for d in sorted(clean_weighted)},
        "direct_joint_distance_birth": {str(d): {str(t): direct[(d, t)] for t in range(2, N + 1)
                                                     if direct[(d, t)]}
                                        for d in range(2, 5)},
        "survival_from_tau2_L": {str(d): {str(t): survival[(d, t)] for t in range(2, N + 1)
                                             if survival[(d, t)]}
                                  for d in range(2, 5)},
        "checks": {"stable_move_bijection": True, "all_slot_words_distinct": True,
                   "matches_r020_distance_birth_table": True,
                   "A2_equals_sum_L": True, "mean_tau_identity": True,
                   "L_weighted_euclidean_square_identity": True,
                   "L_weighted_root_F_identity": True,
                   "L_weighted_root_G_identity": True,
                   "clean_deletion_onepeak_bridge": True},
    }
    (HERE / "second-birth-survival-n3.json").write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
