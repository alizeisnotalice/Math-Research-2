from itertools import product
from math import comb, isqrt
from fractions import Fraction
import platform


def walsh(values):
    values = list(values)
    width = 1
    while width < len(values):
        for start in range(0, len(values), 2 * width):
            for index in range(start, start + width):
                a, b = values[index], values[index + width]
                values[index], values[index + width] = a + b, a - b
        width *= 2
    return values


def character(frequency, point):
    return 1 - 2 * ((frequency & point).bit_count() & 1)


def rank(vectors):
    pivots = {}
    for x in vectors:
        while x:
            i = x.bit_length()
            if i not in pivots:
                pivots[i] = x
                break
            x ^= pivots[i]
    return len(pivots)


def support(values):
    return tuple(i for i, value in enumerate(values) if value)


def affine_rank(points, ambient):
    return rank((1 << ambient) | point for point in points)


def valuation(value):
    return (value & -value).bit_length() - 1 if value else float("inf")


def verify_columns():
    n, N = 3, 8
    checks = 0
    for z in product((-1, 0, 1), repeat=N):
        if not any(z):
            continue
        T = support(z)
        k = len(T)
        spectrum = walsh(z)
        p = len(support(spectrum))
        assert k * p >= N
        assert (N - p) * (1 << affine_rank(T, n)) <= N * comb(k, k // 2)
        if k % 2:
            assert p == N
        moment = [sum(value for u, value in enumerate(z) if u & I == I) for I in range(N)]
        weight = [I.bit_count() + valuation(value) for I, value in enumerate(moment)]
        j = min(weight)
        assert j == min(map(valuation, spectrum))
        active = [I for I in range(N) if weight[I] == j]
        degree = max(I.bit_count() for I in active)
        for alpha in range(N):
            forecast = sum(1 for I in active if alpha & I == I) & 1
            assert forecast == ((spectrum[alpha] >> j) & 1)
        assert p >= (N >> degree)
        checks += 1
    z = [1, 1, 1, 0, 1, 0, 0, 0]
    assert len(support(walsh(z))) == 5
    return checks


def merge(left, right):
    return [a + b for a, b in zip(left, right)] + [a - b for a, b in zip(left, right)]


def one_direction(left, right):
    L, R = support(left), support(right)
    left_rank = affine_rank(L, 4)
    right_rank = affine_rank(R, 4)
    union_rank = affine_rank(set(L) | set(R), 4)
    j = left_rank + right_rank - union_rank
    k = len(set(L) & set(R))
    e = sum(abs(left[i]) == abs(right[i]) for i in set(L) & set(R))
    U = left_rank + right_rank + 2
    assert 2 * (len(L) + len(R)) * (2 * U * j - j * j) < (2 * k + e) * U * U
    return j, k, e


def verify_witness(masks, expected):
    rows = [walsh(1 - 2 * ((mask >> i) & 1) for i in range(8)) for mask in masks]
    table = list(zip(*rows))
    occurrences = [(i, beta) for i, row in enumerate(rows) for beta in support(row)]
    L = len(occurrences)
    U = sum(affine_rank(support(row), 3) + 1 for row in rows)
    labels = [i for i, values in enumerate(table) if any(values)]
    R = affine_rank(labels, 3)
    J = U - 4 - R
    address_supports = {beta: support(walsh(table[beta])) for beta in labels}
    spectrum = [beta + 8 * alpha for beta in labels for alpha in address_supports[beta]]
    parent_rank = affine_rank(spectrum, 5)
    t_eff = parent_rank - R
    Delta = J + 4 - t_eff - 1
    k = tuple(sum(value != 0 for value in table[beta]) for beta in labels)
    n = tuple(len(address_supports[beta]) for beta in labels)
    assert all(isqrt(a * b)**2 == a * b for a, b in zip(k, n))
    Dstar = sum(Fraction(a) - Fraction(isqrt(a * b), 2) for a, b in zip(k, n))
    h = L - rank(((8 | beta) << (4 * u)) | (1 << (16 + beta))
                 for u, beta in occurrences)
    Q = L - len(labels)
    assert J >= Q - h
    jke = [one_direction(merge(rows[i], rows[j]), merge(rows[k1], rows[l]))
           for i, j, k1, l in ((0, 1, 2, 3), (0, 2, 1, 3), (0, 3, 1, 2))]
    assert Fraction(Delta, U) >= Dstar / L
    if masks[0] == 204:
        assert Fraction(Delta, U) < Fraction(Q, L)
    assert (len(spectrum), parent_rank, L, U, J, t_eff, Delta, Dstar) == expected
    assert (parent_rank + 1)**2 < 4 * len(spectrum)
    return {"k": k, "n": n, "shared_nullity": h, "one_bit_profiles": jke}


def check_odd_profiles():
    t = 4
    address_count = 1 << t
    pairs = [(1 << i) | (1 << j) for i in range(t) for j in range(i + 1, t)]
    checks = 0
    for a in range(2, 9):
        data_count = 1 << a
        for m in range(min(a, 5) + 1):
            for no_cheap in (False, True):
                b = [sum((u & pair == pair) << index for index, pair in enumerate(pairs[:m]))
                     for u in range(address_count)]
                eps = [1 - 2 * (no_cheap and u & pairs[m] == pairs[m])
                       for u in range(address_count)]
                phases = [support(walsh(sign * character(beta, value)
                                        for sign, value in zip(eps, b)))
                          for beta in range(data_count)]
                for beta in range(data_count):
                    assert phases[beta] == phases[beta & ((1 << m) - 1)]
                cheap = [beta for beta, band in enumerate(phases) if len(band) == 1]
                assert len(cheap) == (0 if no_cheap else 1 << (a - m))
                diffs = [root ^ band[0] for band in phases for root in band]
                width_rank = rank(diffs)
                values = [eps[u] * (1 - 2 * (x == b[u])) for u in range(address_count)
                          for x in range(data_count)]
                full_support = support(walsh(values))
                expected = tuple(sorted(beta + (alpha << a) for beta, band in enumerate(phases)
                                        for alpha in band))
                assert full_support == expected
                parent_rank = affine_rank(full_support, a + t)
                assert parent_rank == a + 1 + width_rank
                if m or no_cheap:
                    assert 2 * len(full_support) >= 5 * data_count
                if a >= 6 or (a >= 5 and not no_cheap):
                    optimum2 = max(Fraction(4), *(Fraction((affine_rank(band,t) + 1)**2, len(band))
                                                   for band in phases))
                    assert Fraction((parent_rank + 1)**2, len(full_support)) <= optimum2
                checks += 1
    return checks


if __name__ == "__main__":
    print("Python", platform.python_version())
    print("exact integer column profiles:", verify_columns())
    print("first true two-bit witness:", verify_witness((178,17,71,40), (13,6,16,16,8,2,9,9)))
    print("refined true two-bit witness:",
          verify_witness((204,123,156,219), (13,6,13,14,6,2,7,Fraction(13,2))))
    print("exact odd-weight isospectral Boolean profiles:", check_odd_profiles())
