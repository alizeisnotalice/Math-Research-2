"""Exact finite checks for the six-quarter circuit's complementary filling."""
import itertools
import json
import sys
from collections import Counter, defaultdict

Z = tuple(x for x in range(16) if ((x & 3) == 0) == (x >> 2 == 0))


def character(label, point):
    return (-1) ** ((label & point).bit_count() % 2)


def walsh(values):
    result = list(values)
    stride = 1
    while stride < len(result):
        for offset in range(0, len(result), 2 * stride):
            for index in range(offset, offset + stride):
                left, right = result[index], result[index + stride]
                result[index], result[index + stride] = left + right, left - right
        stride *= 2
    return result


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


def arm(point):
    return int((point & 3) == 0) - int(point >> 2 == 0)


def filling(negative):
    return tuple(0 if x not in Z else -1 if x in negative else 1 for x in range(16))


def check_local_fillings():
    distribution = Counter()
    quarter_states = defaultdict(list)
    half_states = defaultdict(list)
    signs = lambda a: tuple(-value for value in a)
    for mask in range(1 << len(Z)):
        negative = [x for index, x in enumerate(Z) if (mask >> index) & 1]
        values = filling(negative)
        coefficients = walsh(values)
        description = tuple(sorted(Counter(abs(x) for x in coefficients if x).items()))
        distribution[description] += 1
        support = [label for label, value in enumerate(coefficients) if value]
        assert rank(label ^ support[0] for label in support) == 4
        if description == ((4, 10),):
            quarter_states[tuple(support)].append(values)
        if description == ((4, 6), (8, 1)):
            odds = tuple(label for label in support if abs(coefficients[label]) == 4)
            peak = next(label for label in support if abs(coefficients[label]) == 8)
            half_states[odds].append((peak, values))
    assert distribution == {
        ((2, 15), (10, 1)): 32, ((4, 6), (8, 1)): 320,
        ((4, 10),): 192, ((2, 13), (6, 3)): 480,
    }
    planes = [
        (0, u, v, u ^ v)
        for u, v in itertools.combinations(Z[1:], 2)
        if u ^ v in Z
    ]
    unique_planes = {tuple(sorted(plane)) for plane in planes}
    assert len(unique_planes) == 6
    for shift in range(16):
        support = tuple(sorted(x ^ shift for x in Z))
        expected = set()
        for plane in unique_planes:
            values = filling(plane[1:])
            modulated = tuple(character(shift, x) * value for x, value in enumerate(values))
            expected.update((modulated, signs(modulated)))
        assert set(quarter_states[support]) == expected
    for shift in range(16):
        support = tuple(sorted(x ^ shift for x in range(16) if x not in Z))
        states = half_states[support]
        assert len(states) == 20 and len({peak for peak, _ in states}) == 10
    return distribution, quarter_states, half_states


def parent(first, second):
    """The length-seven parent has x bits first, then t, then two w bits."""
    values = []
    for point in range(128):
        x = point & 15
        t = character(16, point)
        state = second if point & 32 else first
        outside = t * arm(x)
        on_core = character(64, point) * state[x]
        values.append(outside + on_core)
    assert set(values) == {-1, 1}
    coefficients = walsh(values)
    support = [label for label, coefficient in enumerate(coefficients) if coefficient]
    assert len(support) == 16
    assert {abs(coefficients[label]) for label in support} == {32}
    assert sum(bool(label & 16) for label in support) == 6
    return rank(support)


def check_pairs(quarter_states, half_states):
    quarter = Counter()
    half = Counter()
    for states in quarter_states.values():
        for first, second in itertools.combinations(states, 2):
            dimension = parent(first, second)
            assert dimension == (6 if all(a == -b for a, b in zip(first, second)) else 7)
            quarter[dimension] += 1
    for states in half_states.values():
        for (peak0, first), (peak1, second) in itertools.combinations(states, 2):
            if peak0 == peak1:
                continue
            dimension = parent(first, second)
            assert dimension == 7
            half[dimension] += 1
    assert parent(filling((6, 9, 15)), filling((5, 10, 15))) == 7
    return quarter, half


def check_unbounded_fillings():
    checked = []
    for n in range(1, 9):
        seed = tuple(-1 if w == 0 else 1 for w in range(1 << n))
        seed_spectrum = [label for label, value in enumerate(walsh(seed)) if value]
        values = []
        for point in range(1 << (5 + n)):
            x = point & 15
            values.append(character(16, point) * arm(x) if x not in Z else seed[point >> 5])
        spectrum = [label for label, value in enumerate(walsh(values)) if value]
        assert len(spectrum) == 6 + 16 * len(seed_spectrum)
        assert rank(spectrum) == 5 + rank(seed_spectrum)
        checked.append({"seed_variables": n, "s": len(spectrum), "d": rank(spectrum)})
    return checked


def check_independent_children():
    checked = []
    for occupied in range(11):
        values = []
        for point in range(1 << (5 + occupied)):
            core = point & 15
            if core not in Z:
                values.append(character(16, point) * arm(core))
                continue
            child = Z.index(core)
            if child >= occupied:
                values.append(1)
            else:
                values.append(character(1 << (5 + child), point))
        spectrum = [label for label, value in enumerate(walsh(values)) if value]
        dimension = rank(spectrum)
        nonzero = len(spectrum) - int(0 in spectrum)
        assert dimension == 5 + occupied
        assert nonzero >= 6 + 16 * occupied
        assert (dimension + 1) ** 2 <= 4 * nonzero
        checked.append({"occupied_children": occupied, "q": nonzero, "d": dimension})
    return checked


if __name__ == "__main__":
    distribution, quarter_states, half_states = check_local_fillings()
    quarter, half = check_pairs(quarter_states, half_states)
    print(json.dumps({
        "python": sys.version.split()[0],
        "local_count": sum(distribution.values()),
        "local_profiles": [{"magnitude_histogram": dict(k), "count": v} for k, v in sorted(distribution.items())],
        "quarter_state_pair_ranks": quarter,
        "half_state_pair_ranks": half,
        "unrestricted_fillings": check_unbounded_fillings(),
        "independent_children": check_independent_children(),
    }, indent=2))
