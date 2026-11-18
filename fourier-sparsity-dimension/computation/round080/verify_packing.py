"""Exact Walsh, component and density checks for the two-block constructions."""
import json
import sys
from fractions import Fraction
from functools import cache


def walsh(values):
    result = list(values)
    stride = 1
    while stride < len(result):
        for start in range(0, len(result), 2 * stride):
            for index in range(start, start + stride):
                left, right = result[index], result[index + stride]
                result[index], result[index + stride] = left + right, left - right
        stride *= 2
    return result


def spectrum(values):
    return tuple(label for label, value in enumerate(walsh(values)) if value)


def rank(labels):
    pivots = {}
    for label in labels:
        while label:
            pivot = label.bit_length()
            if pivot not in pivots:
                pivots[pivot] = label
                break
            label ^= pivots[pivot]
    return len(pivots)


def bit_points(mask):
    while mask:
        least = mask & -mask
        yield least
        mask ^= least


def affine_components(labels, dimension):
    labels = tuple(sorted(labels))
    representatives = list(range(len(labels)))
    pivots = {}

    def root(index):
        while representatives[index] != index:
            representatives[index] = representatives[representatives[index]]
            index = representatives[index]
        return index

    for index, label in enumerate(labels):
        augmented, relation = label | (1 << dimension), 1 << index
        while augmented:
            pivot = augmented.bit_length()
            if pivot not in pivots:
                pivots[pivot] = augmented, relation
                break
            vector, combination = pivots[pivot]
            augmented ^= vector
            relation ^= combination
        if not augmented:
            for bit in bit_points(relation):
                representatives[root(bit.bit_length() - 1)] = root(index)
    groups = {}
    for index, label in enumerate(labels):
        groups.setdefault(root(index), set()).add(label)
    return {frozenset(group) for group in groups.values()}


@cache
def subspaces(dimension):
    everything = (1 << (1 << dimension)) - 1
    groups = {1}
    result = []
    for dim in range(dimension + 1):
        result.extend((space, dim) for space in groups)
        expanded = set()
        for space in groups:
            outside = everything ^ space
            current_points = [bit.bit_length() - 1 for bit in bit_points(space)]
            while outside:
                representative = (outside & -outside).bit_length() - 1
                coset = sum(1 << (representative ^ point) for point in current_points)
                expanded.add(space | coset)
                outside ^= coset
        groups = expanded
    return result


def brute_density(support, dimension):
    points = sum(1 << label for label in support if label)
    return max(
        Fraction((points & space).bit_count(), dim)
        for space, dim in subspaces(dimension)
        if dim > 0
    )


@cache
def profile(labels, lift):
    result = {}
    for mask in range(1 << len(labels)):
        subset = [
            label | lift
            for index, label in enumerate(labels)
            if mask & (1 << index)
        ]
        dim = rank(subset)
        result[dim] = max(result.get(dim, 0), len(subset))
    return result


def formula_density(g_spectrum, j_spectrum, g_dimension):
    if g_spectrum == tuple(range(1 << g_dimension)):
        g_profile = {0: 0, **{dim: 1 << (dim - 1) for dim in range(1, g_dimension + 2)}}
    else:
        g_profile = profile(g_spectrum, 1 << g_dimension)
    j_dimension = j_spectrum[-1].bit_length()
    if j_spectrum == tuple(range(1 << j_dimension)):
        j_profile = {dim: (1 << dim) - 1 for dim in range(j_dimension + 1)}
    else:
        j_profile = profile(tuple(x for x in j_spectrum if x), 0)
    has_zero = int(0 in j_spectrum)
    return max(
        Fraction(2 * g_count + 2 * j_count + has_zero, 1 + g_dim + j_dim)
        for g_dim, g_count in g_profile.items()
        for j_dim, j_count in j_profile.items()
    )


def quadratic(basis_dimension):
    return tuple(
        (-1) ** sum(((x >> index) & 3) == 3 for index in range(0, basis_dimension, 2))
        for x in range(1 << basis_dimension)
    )


def check_seed(g, j, brute=False):
    g_dim = len(g).bit_length() - 1
    j_dim = len(j).bit_length() - 1
    dimension = 2 + g_dim + j_dim
    t, w = 1, 2
    tag_offset = 2 + g_dim
    core_mask = ((1 << j_dim) - 1) << tag_offset
    parent = []
    for x in range(1 << dimension):
        if x & t:
            g_val = g[(x >> 2) & ((1 << g_dim) - 1)]
            parent.append(-g_val if x & w else g_val)
        else:
            parent.append(j[x >> tag_offset])
    actual = spectrum(parent)
    g_spectrum, j_spectrum = spectrum(g), spectrum(j)
    big = frozenset(w | (label << 2) | offset for label in g_spectrum for offset in (0, t))
    small = frozenset(t | (label << tag_offset) for label in j_spectrum)
    assert set(actual) == set(big) | set(small) | {label << tag_offset for label in j_spectrum}
    assert len(actual) == 2 * len(g_spectrum) + 2 * len(j_spectrum)
    assert rank(actual) == 1 + rank(j_spectrum) + rank(label | (1 << g_dim) for label in g_spectrum)
    if len(g_spectrum) >= 2 and (len(j_spectrum) == 1 or len(j_spectrum) >= 4):
        assert affine_components(big | small, dimension) == {big, small}
    assert rank(label & ~core_mask for label in big) == rank(label & ~core_mask for label in actual)
    density = formula_density(g_spectrum, j_spectrum, g_dim)
    if brute:
        assert density == brute_density(actual, dimension)
    return {"s": len(actual), "d": rank(actual), "density": str(density)}


def check_imbalance(y_dim):
    g = j = quadratic(2)
    e = quadratic(y_dim)
    dimension = 6 + y_dim
    parent = []
    for point in range(1 << dimension):
        t_sign = (-1) ** (point & 1)
        w_sign = (-1) ** ((point >> 1) & 1)
        g_value = g[(point >> 2) & 3]
        j_value = j[(point >> 4) & 3]
        e_value = e[point >> 6]
        parent.append(w_sign * g_value if t_sign == e_value else e_value * j_value)
    coefficients = walsh(parent)
    actual = {label for label, value in enumerate(coefficients) if value}
    big = {label for label in actual if label & 2}
    small = {label for label in actual if not label & 2 and label & 1}
    assert affine_components(big | small, dimension) == {frozenset(big), frozenset(small)}
    lifts = {
        base: sum((label & 15) == base for label in big)
        for base in sorted({label & 15 for label in big})
    }
    for base in range(2, 16, 4):
        assert lifts[base] == 1 and lifts[base | 1] == len(spectrum(e))
    assert all(abs(coefficients[label]) * 4 == 1 << dimension for label in small)
    return {"switch_variables": y_dim, "pair_lifts": [1, len(spectrum(e))], "s": len(actual)}


if __name__ == "__main__":
    truth_tables = [
        tuple((-1) ** ((mask >> x) & 1) for x in range(8))
        for mask in range(256)
    ]
    for g in truth_tables:
        check_seed(g, (1,), brute=True)
        check_seed(g, quadratic(2))
    equal = []
    for variables in (2, 4, 6):
        actual = check_seed(quadratic(variables), quadratic(variables), brute=variables == 2)
        expected = Fraction(4 * (1 << variables) - 1, 2 * variables + 2)
        assert actual["density"] == str(expected)
        equal.append({"variables_per_seed": variables, **actual})
    print(json.dumps({
        "python": sys.version.split()[0],
        "all_three_variable_seed_tests": 2 * len(truth_tables),
        "three_variable_seed_density_bruteforce_tests": len(truth_tables),
        "equal_full_seeds": equal,
        "imbalanced_examples": [check_imbalance(dim) for dim in (2, 4, 6, 8)],
    }, indent=2))
