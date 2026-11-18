"""Integer Walsh checks for the ternary low-norm theorem. Uses only Python's standard library."""
from collections import Counter
from fractions import Fraction
from itertools import combinations, product
import json
import platform


def walsh(values):
    spectrum = list(values)
    step = 1
    while step < len(spectrum):
        for base in range(0, len(spectrum), 2 * step):
            for i in range(base, base + step):
                a, b = spectrum[i], spectrum[i + step]
                spectrum[i], spectrum[i + step] = a + b, a - b
        step *= 2
    return spectrum


def small_census(n):
    size = 1 << n
    counts = Counter(sum(map(abs, walsh(h))) for h in product((-1, 0, 1), repeat=size))
    assert all(4 * s >= 7 * size or 2 * s in (0, 2 * size, 3 * size) for s in counts)
    return {str(Fraction(s, size)): counts[s] for s in sorted(counts) if s <= 2 * size}


def normal_form(codim):
    size = 4 << codim
    h = [int((x & 3) == 0) - int((x >> 2) == 0) for x in range(size)]
    spectrum = walsh(h)
    magnitudes = sorted({Fraction(abs(s), size) for s in spectrum if s})
    norm = Fraction(sum(map(abs, spectrum)), size)
    assert norm == (Fraction(3, 2) if codim <= 2 else 2 - Fraction(2, 1 << codim))
    if codim == 2:
        support = {x for x, s in enumerate(spectrum) if s}
        assert support == {1, 2, 3, 4, 8, 12}
        assert {x ^ y for x, y in combinations(support, 2)} == set(range(1, 16))
    return {"norm": str(norm), "nonzero_magnitudes": list(map(str, magnitudes))}


def boolean_example():
    def linear_c(v):
        return (v >> 1) | ((((v >> 1) ^ v) & 1) << 1)

    def sign_phase(u, v):
        return (-1) ** ((u & v).bit_count() % 2)

    branches = []
    for apply_c in (False, True):
        s = [sign_phase(x & 3, linear_c(x >> 2) if apply_c else x >> 2)
             for x in range(16)]
        inverse = walsh(s)
        assert set(inverse) == {-4, 4}
        branches.extend(z // 4 for z in inverse)
    spectrum = walsh(branches)
    assert set(spectrum) == {-8, 0, 8}
    halves = [len([s for s in spectrum[i:i + 16] if s]) for i in (0, 16)]
    assert halves == [10, 6]
    points = [i for i, s in enumerate(spectrum[16:]) if s]
    assert len({x ^ y for x, y in combinations(points, 2)}) == 15
    return {"total_spectrum": sum(halves), "half_weights": halves, "six_point_cut": points}


print(json.dumps({
    "python": platform.python_version(),
    "normalization": "Walsh sums divided by table length",
    "all_ternary_tables": {str(n): {"count": 3 ** (1 << n), "norm_through_two": small_census(n)}
                           for n in (2, 3)},
    "normal_forms": {str(d): normal_form(d) for d in range(5)},
    "boolean_example": boolean_example(),
}, indent=2))
