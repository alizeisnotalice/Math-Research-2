"""Finite exact checks of the sharp projected rank theorem and the embedding."""

from collections import Counter
from fractions import Fraction
import json
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


def gf2_rank(vectors):
    pivots = {}
    for vector in vectors:
        while vector:
            pivot = vector.bit_length() - 1
            if pivot not in pivots:
                pivots[pivot] = vector
                break
            vector ^= pivots[pivot]
    return len(pivots)


def bit(value, index):
    return (value >> index) & 1


def sharp_exact(m):
    r = m.bit_length() - 1
    assert m == 1 << r
    selector_bits = r + m
    bands = {}
    for address in range(m):
        # e, f, t in the address's independent data block.
        local_labels = [(1 << (3 * address + 2)) | value
                        for value in (0, 1 << (3 * address), 1 << (3 * address + 1), 3 << (3 * address))]
        for label in local_labels:
            # Labeled by the private sign bit u_i and any address frequency.
            bands[label] = [frequency | (1 << (r + address)) for frequency in range(m)]
    full_points = [private_frequency | (label << selector_bits)
                   for label, frequencies in bands.items() for private_frequency in frequencies]
    projected_augmented_rank = gf2_rank(label | (1 << (3 * m)) for label in bands)
    assert len(full_points) == 4 * m * m
    assert gf2_rank(full_points) == r + 3 * m
    assert projected_augmented_rank == 3 * m
    assert 4 * projected_augmented_rank**2 == 9 * len(full_points)

    if m <= 4:
        # Recompute all physical data bands on the complete selector cube.
        for label, expected in bands.items():
            i = next(i for i in range(m) if label & (1 << (3 * i + 2)))
            values = [((-1) ** bit(x, r + i)) if (x & (m - 1)) == i else 0
                      for x in range(1 << selector_bits)]
            actual = {frequency for frequency, coefficient in enumerate(walsh(values)) if coefficient}
            assert actual == set(expected)

    return {"addresses": m, "projected_augmented_rank": projected_augmented_rank,
            "full_sparsity": len(full_points), "full_linear_rank": r + 3 * m,
            "per_label_actual_selector_assignments": 1 << m}


def repeated_types(m, t):
    assert m & (m - 1) == 0 and t & (t - 1) == 0 and 1 < t < m
    r = m.bit_length() - 1
    base = 1 << (3 * m)
    signals = {}
    for i in range(m):
        common_or_private = base if i < t else 1 << (3 * i + 2)
        for extra in (0, 1 << (3 * i), 1 << (3 * i + 1), 3 << (3 * i)):
            label = common_or_private | extra
            # Root label coefficient is +1/2; other signs are irrelevant
            # for disjoint private labels. Integer scale omits this half.
            signals.setdefault(label, [0] * m)[i] = -1 if extra == 3 << (3 * i) else 1
    supports = {label: [frequency for frequency, value in enumerate(walsh(signal)) if value]
                for label, signal in signals.items()}
    assert len(supports[base]) == m // t
    s = sum(map(len, supports.values()))
    projected = gf2_rank(label | (base << 1) for label in supports)
    full_rank = gf2_rank(frequency | (label << r)
                        for label, frequencies in supports.items() for frequency in frequencies)
    assert s == (4 * m - t) * m + m // t
    assert projected == 3 * m - t + 1 and full_rank == projected + r
    return {"addresses": m, "distinct_slice_types_with_common_label": t,
            "common_label_sparsity": len(supports[base]), "s": s,
            "projected_augmented_rank": projected,
            "projected_ratio_squared": str(Fraction(projected**2, s))}


def embedding():
    result = []
    for G in ([1, -1, 1, -1], [1, 1, 1, 1, 1, -1, -1, -1]):
        size = len(G)
        spectrum = [frequency for frequency, value in enumerate(walsh(G)) if value]
        a_G = gf2_rank(frequency | size for frequency in spectrum) - 1
        for a in (1, 2, 3, 4):
            values = []
            for data in range(1 << a):
                u, z = data & ((1 << (a - 1)) - 1), bit(data, a - 1)
                for g in G:
                    values.append(1 - 2 * int(u == 0 and (-1) ** z == g))
            support = [frequency for frequency, value in enumerate(walsh(values)) if value]
            d_F = gf2_rank(support)
            assert d_F == a_G + a
            assert len(support) == (len(spectrum) if a == 1 else (1 << (a - 1)) * (len(spectrum) + 1))
            result.append({"dimension_X": size.bit_length() - 1, "a": a, "s_G": len(spectrum),
                           "a_G": a_G, "s_F": len(support), "d_F": d_F})
    return result


def main():
    output = {"python": platform.python_version(),
              "exact_ratio_family": [sharp_exact(m) for m in (1, 2, 4, 16, 64)],
              "repeated_label_in_distinct_types": [repeated_types(m, t) for m, t in ((16, 4), (64, 8), (256, 16))],
              "embedding": embedding()}
    print(json.dumps(output, indent=2))


if __name__ == "__main__":
    main()
