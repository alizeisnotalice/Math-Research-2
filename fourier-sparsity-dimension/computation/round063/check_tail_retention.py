"""Check dense highest layers and exact tail contraction witnesses over F_2."""
import json
import platform
from collections import Counter
from math import gcd


def sign(bits):
    return -1 if bits.bit_count() & 1 else 1


def walsh(values):
    width = 1
    while width < len(values):
        for offset in range(0, len(values), width * 2):
            for index in range(offset, offset + width):
                first, second = values[index], values[index + width]
                values[index] = first + second
                values[index + width] = first - second
        width *= 2
    return values


def sparse_walsh(values):
    return {frequency: count for frequency, count in enumerate(walsh(values)) if count}


def bent_value(a, z):
    half = a // 2
    return sign((z >> half) & (z & ((1 << half) - 1)))


def bent_case(a):
    m = 1 << a
    n_inner = a + 1
    n_outer = a + 3
    values = []
    for index in range(1 << (n_inner + n_outer)):
        y, v = index & (m - 1), (index >> a) & 1
        z = (index >> n_inner) & (m - 1)
        t, u, w = ((index >> (n_inner + a + j)) & 1 for j in range(3))
        if not v:
            value = sign(w) if y & 1 else (-1 if y == 0 else 1)
        else:
            value = sign(u) * bent_value(a, z ^ (y if t else 0))
        values.append(value)
    return values, n_inner


def address_case(a):
    m = 1 << a
    n_inner = a + 3
    n_outer = m + 1
    values = []
    for index in range(1 << (n_inner + n_outer)):
        y, b = index & (m - 1), (index >> a) & 3
        v = (index >> (a + 2)) & 1
        w = (index >> (n_inner + m)) & 1
        if not v:
            value = sign(w) if y & 1 else (-1 if y == b == 0 else 1)
        else:
            value = sign((index >> (n_inner + y)) & 1)
        values.append(value)
    return values, n_inner


def verify(a, name):
    m = 1 << a
    values, n_inner = (bent_case if name == "bent" else address_case)(a)
    size_inner = 1 << n_inner
    size_outer = len(values) // size_inner
    spectrum = sparse_walsh(values.copy())
    divisor = gcd(*spectrum.values())
    odd = {frequency for frequency, count in spectrum.items() if (count // divisor) & 1}
    assert odd == set(range(size_inner))
    divisor_power = (len(values) // divisor).bit_length() - 1
    assert divisor_power == (a if name == "bent" else a + 2)
    histogram = Counter(frequency >> n_inner for frequency in spectrum if frequency not in odd)
    supports = {}
    if name == "bent":
        assert sum(histogram.values()) == 8 * (m - 1) + 6
        assert len(histogram) == 2 * m
    else:
        assert sum(histogram.values()) == 2 * m * m + 4
        assert len(histogram) == m + 1

    for x in range(size_inner):
        sliced = walsh(values[x::size_inner])
        if name == "bent":
            y, v = x & (m - 1), (x >> a) & 1
            expected = (
                {beta | (((beta & y).bit_count() & 1) << a) | (2 * m):
                 (size_outer >> (a // 2)) * bent_value(a, beta) for beta in range(m)}
                if v else ({4 * m: size_outer} if y & 1 else {})
            )
        else:
            y, v = x & (m - 1), (x >> (a + 2)) & 1
            expected = (
                {1 << y: size_outer} if v
                else ({1 << m: size_outer} if y & 1 else {})
            )
        actual = {frequency: count for frequency, count in enumerate(sliced) if count and frequency}
        assert actual == expected
        for beta in actual:
            supports[beta] = supports.get(beta, 0) | (1 << x)

    witnesses = 0
    if name == "bent":
        for normal in range(1, size_inner):
            first_half = sum(1 << x for x in range(size_inner) if sign(normal & x) == 1)
            second_half = ((1 << size_inner) - 1) ^ first_half
            assert any(first_half & support == 0 for support in supports.values())
            assert any(second_half & support == 0 for support in supports.values())
            witnesses += 2
    else:
        # b fixed to zero, all y and both v survive; y must be surjective
        # in the v=1 section, and the extra label forces meeting v=0.
        b_zero = sum(1 << x for x in range(size_inner) if ((x >> a) & 3) == 0)
        assert all(b_zero & support for support in supports.values())
        assert all(support.bit_count() == 4 for support in supports.values() if support != supports[1 << m])
        assert supports[1 << m] & ((1 << (4 * m)) - 1) == supports[1 << m]

    return {
        "case": name,
        "a": a,
        "ambient_variables": (len(values)).bit_length() - 1,
        "highest_denominator_exponent": divisor_power,
        "odd_layer_points": len(odd),
        "tail_points": sum(histogram.values()),
        "outer_labels": len(histogram),
        "local_sparsity_histogram": dict(sorted(Counter(histogram.values()).items())),
        "affine_hyperplanes_checked": witnesses,
        "minimum_retention_dimension": a + 1,
    }



def verify_affine_basis(m, padding):
    marked = {0: 0, **{1 << j: j + 1 for j in range(m)}}
    n_inner = m + padding
    n_outer = m + 1
    size_inner, size_outer = 1 << n_inner, 1 << n_outer
    values = []
    for index in range(size_inner * size_outer):
        x = index & ((1 << m) - 1)
        if x in marked:
            value = sign((index >> (n_inner + marked[x])) & 1)
        else:
            value = -1 if (index & (size_inner - 1)) == 3 else 1
        values.append(value)

    spectrum = sparse_walsh(values.copy())
    divisor = gcd(*spectrum.values())
    odd = {frequency for frequency, count in spectrum.items() if (count // divisor) & 1}
    assert odd == set(range(size_inner))
    divisor_power = (len(values) // divisor).bit_length() - 1
    assert divisor_power == m + padding - 1
    histogram = Counter(frequency >> n_inner for frequency in spectrum if frequency not in odd)
    assert histogram == {1 << j: 1 << m for j in range(m + 1)}
    for inner in range(size_inner):
        sliced = walsh(values[inner::size_inner])
        actual = {frequency: count for frequency, count in enumerate(sliced) if count and frequency}
        x = inner & ((1 << m) - 1)
        expected = {1 << marked[x]: size_outer} if x in marked else {}
        assert actual == expected
    return {
        "case": "affine_basis",
        "m": m,
        "padding": padding,
        "ambient_variables": n_inner + n_outer,
        "highest_denominator_exponent": divisor_power,
        "odd_layer_points": len(odd),
        "tail_points": sum(histogram.values()),
        "outer_labels_and_rank": len(histogram),
        "minimum_retention_dimension": m,
    }


if __name__ == "__main__":
    rows = [verify(a, "address") for a in (1, 2, 3)]
    rows.extend(verify(a, "bent") for a in (6, 8))
    rows.extend(verify_affine_basis(m, p) for m, p in ((2, 2), (4, 2), (6, 3)))
    print(json.dumps({"python": platform.python_version(), "arithmetic": "integer", "checks": rows}, indent=2))
