from collections import Counter
from functools import lru_cache
from fractions import Fraction
from itertools import combinations
import platform


def character(frequency, point):
    return 1 - 2 * ((frequency & point).bit_count() & 1)


def fwht(values):
    values = list(values)
    width = 1
    while width < len(values):
        for start in range(0, len(values), 2 * width):
            for index in range(start, start + width):
                left, right = values[index], values[index + width]
                values[index], values[index + width] = left + right, left - right
        width *= 2
    return values


@lru_cache(None)
def affine_nullity(support, n):
    pivots = {}
    relations = []
    for index, frequency in enumerate(support):
        column = (1 << n) | frequency
        relation = 1 << index
        while column:
            top = column.bit_length()
            if top not in pivots:
                pivots[top] = column, relation
                break
            other_column, other_relation = pivots[top]
            column ^= other_column
            relation ^= other_relation
        else:
            relations.append(relation)
    return len(relations), relations[0] if len(relations) == 1 else 0


def profile(abs_four_times, negative_count, coloop_four_times):
    key = tuple(abs_four_times), negative_count % 2, tuple(coloop_four_times)
    return {
        ((1, 1, 1, 1), 0, ()): "A",
        ((1, 1, 1, 1, 1, 1), 1, ()): "B",
        ((1, 1, 1, 3), 1, ()): "C",
        ((1, 1, 1, 1), 1, (2,)): "D",
        ((2, 2, 2, 2), 1, ()): "E",
    }[key]


def profiles_four_dimensional():
    n, N = 4, 16
    characters = [
        sum(1 << point for point in range(N) if character(frequency, point) < 0)
        for frequency in range(N)
    ]
    hyperplanes = [
        tuple(frequency for frequency in range(N) if character(normal, frequency) == sign)
        for normal in range(1, N)
        for sign in (-1, 1)
    ]
    counts = Counter()
    for truth in range(1 << N):
        coeff = [N - 2 * (truth ^ mask).bit_count() for mask in characters]
        for plane in hyperplanes:
            support = tuple(freq for freq in plane if coeff[freq])
            if not support:
                continue
            nullity, rel_mask = affine_nullity(support, n)
            if nullity != 1:
                continue
            assert all(coeff[freq] % 4 == 0 for freq in support)
            circuit = [
                coeff[freq] // 4 for i, freq in enumerate(support) if (rel_mask >> i) & 1
            ]
            coloops = [
                coeff[freq] // 4 for i, freq in enumerate(support) if not ((rel_mask >> i) & 1)
            ]
            name = profile(sorted(map(abs, circuit)), sum(x < 0 for x in circuit), sorted(map(abs, coloops)))
            counts[name] += 1
    return counts


def realize_profile(name, circuit_weights, coloop_weights):
    length = len(circuit_weights)
    dim = length - 2 + len(coloop_weights)
    circuit = (0, *(1 << i for i in range(length - 2)), (1 << (length - 2)) - 1)
    coloops = tuple(1 << i for i in range(length - 2, dim))
    frequencies = circuit + coloops
    weights = circuit_weights + coloop_weights
    signed = [sum(weight * character(freq, y) for weight, freq in zip(weights, frequencies))
              for y in range(1 << dim)]
    assert all(value in (-4, 0, 4) for value in signed)
    zero = tuple(y for y, value in enumerate(signed) if value == 0)
    values = [1 if value == 0 else value // 4 for value in signed]
    values += [1 if value == 0 else -value // 4 for value in signed]
    raw = fwht(values)
    original_band = raw[1 << dim:]
    assert set(i for i, v in enumerate(original_band) if v) == set(frequencies)
    for freq, weight in zip(frequencies, weights):
        assert original_band[freq] * 4 == weight * len(values)
    assert affine_nullity(tuple(sorted(frequencies)), dim)[0] == 1
    assert len(zero) == {"A": 3, "B": 10, "C": 1, "D": 4, "E": 0}[name]
    if name in ("B", "D"):
        for normal in range(1, 1 << dim):
            for sign in (-1, 1):
                flat = {p for p in range(1 << dim) if character(normal, p) == sign}
                assert not flat <= set(zero)
    if name == "B":
        assert all(set(zero) != {point ^ shift for point in zero} for shift in range(1, 1 << dim))
    return len(zero), len(signed)


def verify_sharp_codim(t):
    length = (1 << (t + 1)) + 2
    d = length - 2
    if t <= 2:
        basis_mask = (1 << d) - 1
        spectrum = (0, *(1 << i for i in range(d)), basis_mask)
        assert affine_nullity(spectrum, d)[0] == 1
        h_num = [
            -1 + sum(character(freq, y) for freq in spectrum[1:])
            for y in range(1 << d)
        ]
        all_bits = (1 << t)
        values = []
        for aux in range(all_bits):
            for num in h_num:
                assert num % 4 == 0 and abs(num) <= 2 * all_bits
                values.append(1 if aux < (2 * all_bits + num) // 4 else -1)
        transform = fwht(values)
        band = transform[:1 << d]
        assert tuple(i for i, val in enumerate(band) if val) == tuple(sorted(spectrum))
        assert all(abs(band[i]) == (len(values) >> (t + 1)) for i in spectrum)
    else:
        for negatives in range(1, length, 2):
            numerator = length - 2 * negatives
            assert numerator % 4 == 0 and abs(numerator) <= 1 << (t + 1)
    return length


def verify_sharp_denominators():
    cases = 0
    for k in range(1, 8):
        d = k + 1
        raw = fwht([-1] + [1] * ((1 << d) - 1))
        assert all(value % 4 == 2 for value in raw)
        for t in range(k):
            for prefix in range(1 << t):
                support = tuple((prefix << (d - t)) + i for i in range(1 << (d - t)))
                K = k - t
                assert affine_nullity(support, d)[0] == (1 << (K + 1)) - K - 2
                cases += 1
    return cases


def verify_record_algebra():
    checks = 0
    for i in range(21):
        t = Fraction(i, 20)
        assert (2 - t) ** 2 - (3 - 3 * t * t) == (2 * t - 1) ** 2
        for j in range(21):
            delta = Fraction(j, 20)
            difference = (2 - t + delta) ** 2 - (Fraction(5, 2) - 3 * t * t + Fraction(35, 2) * delta)
            completed_square = Fraction(1, 2) - Fraction(29, 2) * delta + Fraction(3, 4) * delta * delta
            completed_square += 4 * (t - Fraction(1, 2) - delta / 4) ** 2
            assert difference == completed_square
            checks += 1
    return checks


def main():
    print("Python", platform.python_version())
    print("four-variable complete truth tables:", 1 << 16)
    print("nullity-one hyperplane counts:", dict(sorted(profiles_four_dimensional().items())))
    for name, circ, coloops in (
        ("A", (1, 1, 1, 1), ()),
        ("B", (-1, 1, 1, 1, 1, 1), ()),
        ("C", (-3, 1, 1, 1), ()),
        ("D", (-1, 1, 1, 1), (2,)),
        ("E", (-2, 2, 2, 2), ()),
    ):
        print("realized profile", name, realize_profile(name, circ, coloops))
    print("sharp codimension t=1..8, complete transform for t=1,2:",
          [verify_sharp_codim(t) for t in range(1, 9)])
    print("sharp denominator/nullity flat cuts:", verify_sharp_denominators())
    print("rational normalizations checked:", verify_record_algebra())


if __name__ == "__main__":
    main()
