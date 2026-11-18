"""Exact integer checks for odd data weights and address Fourier circuits."""
from itertools import combinations, product
import json
import platform


def walsh(values):
    result = list(values)
    step = 1
    while step < len(values):
        for start in range(0, len(values), 2 * step):
            for i in range(start, start + step):
                a, b = result[i], result[i + step]
                result[i], result[i + step] = a + b, a - b
        step *= 2
    return result


def parity(x):
    return x.bit_count() % 2


def sign(f, u):
    return 1 - 2 * parity(f & u)


def xor_all(array):
    acc = 0
    for x in array:
        acc ^= x
    return acc


def rank(array):
    rows = {}
    for value in array:
        while value:
            bit = value.bit_length()
            if bit not in rows:
                rows[bit] = value
                break
            value ^= rows[bit]
    return len(rows)


def circuit(space, support):
    for size in range(4, len(support) + 1, 2):
        for sub in combinations(support, size):
            if not xor_all(sub):
                return sub
    return None


def check_circuit(space, support, sub):
    alpha = sub[0]
    extras = sorted(set(support) - set(sub))
    values = []
    for u in range(space):
        weight = sum(parity((a ^ alpha) & u) for a in sub[1:])
        assert weight % 2 == 0
        values.append(sign(alpha, u) * (1 - weight) +
                      4 * sum(sign(a, u) for a in extras))
    coeffs = walsh(values)
    assert [i for i, c in enumerate(coeffs) if c] == list(support)
    phase = [((values[0] - v) // 2) % 2 for v in values]
    bits = range(space.bit_length() - 1)
    polar = [sum((phase[(1 << i) ^ (1 << j)] ^ phase[1 << i] ^
                  phase[1 << j]) << j for j in bits) for i in bits]
    assert rank(polar) == len(sub) - 2
    n_data = max(2, max(map(abs, values)).bit_length() + 1)
    n = 1 << n_data
    for a in values:
        weight = n // 2 - a
        assert 0 <= weight <= n and weight % 2 == 1
        assert sum([-1] * weight + [1] * (n - weight)) == 2 * a


def check_independent(space, support):
    tests = 0
    for coeffs in product((-2, -1, 1, 2), repeat=len(support)):
        if sum(coeffs) % 2 == 0:
            continue
        label = xor_all(a for a, c in zip(support, coeffs) if c % 2)
        values = [sum(c * sign(a, u) for a, c in zip(support, coeffs))
                  for u in range(space)]
        assert all((v - values[0] * sign(label, u)) % 4 == 0
                   for u, v in enumerate(values))
        tests += 1
    return tests


def check_two_lifts():
    data = range(8)
    slices = [tuple(1 - 2 * ((word >> x) & 1) for x in data)
              for word in range(256) if parity(word)]
    seen = {}
    for pair in product(slices, repeat=2):
        # Embed an arbitrary pair of genuine odd slices in a four-address
        # selector using a fixed controlled translation and modulation.
        rows = [walsh([sign(2, u) * pair[parity(u & 3)][x ^ (3 if u & 2 else 0)]
                       for x in data]) for u in range(4)]
        labels = []
        widths = []
        for xi in data:
            values = [r[xi] // 2 for r in rows]
            scaled = walsh(values)
            assert all(c % 4 == 0 for c in scaled)
            coeffs = [c // 4 for c in scaled]
            width = sum(c != 0 for c in coeffs)
            assert 1 <= width <= 2
            labels.append(xor_all(i for i, c in enumerate(coeffs) if c % 2))
            assert coeffs[labels[-1]] % 2
            widths.append(width)
        for xi in data:
            for eta in data:
                assert labels[xi ^ eta] == labels[xi] ^ labels[eta] ^ labels[0]
        key = str(sum(w == 2 for w in widths))
        seen[key] = seen.get(key, 0) + 1
    return seen


def check_heavy_partitions():
    counts = {}
    total = 0

    def visit(word):
        nonlocal total
        if len(word) < 8:
            for label in range(max(word) + 2):
                visit(word + [label])
            return
        total += 1
        normals = []
        for label in sorted(set(word)):
            subset = [u for u, v in enumerate(word) if v == label]
            scaled = walsh([4 * int(v == label) - 1 for v in word])
            support = [i for i, c in enumerate(scaled) if c]
            if len(support) > 3:
                return
            direction_rank = rank([u ^ subset[0] for u in subset])
            assert 2 ** direction_rank == len(subset) and direction_rank >= 1
            normals.extend(a ^ support[0] for a in support)
        r = rank(normals)
        assert r <= min(3, len(set(word)) - 1)
        counts[r] = counts.get(r, 0) + 1

    visit([0])
    assert total == 4140
    return {"tested_partitions": total, "accepted_by_address_rank": counts}


def main():
    independent = dependent = polynomials = 0
    space = 8
    for size in range(1, space + 1):
        for support in combinations(range(space), size):
            sub = circuit(space, support)
            if sub is None:
                assert rank([space | a for a in support]) == len(support)
                independent += 1
                polynomials += check_independent(space, support)
            else:
                dependent += 1
                check_circuit(space, support, sub)
    lengths = []
    for length in (4, 6, 8, 10, 12):
        dimension = length - 2
        support = tuple(sorted([0, (1 << dimension) - 1] +
                               [1 << i for i in range(dimension)]))
        check_circuit(1 << dimension, support, support)
        lengths.append(length)
    counts = check_two_lifts()
    out = dict(python=platform.python_version(), tolerance=0,
               address_three_bit_independent_supports=independent,
               independent_integer_polynomials=polynomials,
               address_three_bit_dependent_supports=dependent,
               circuit_lengths=lengths, three_bit_heavy_partitions=check_heavy_partitions(),
               constructed_selectors_by_number_of_two_lift_bands=counts)
    print(json.dumps(out, indent=2, sort_keys=True))


if __name__ == '__main__':
    main()
