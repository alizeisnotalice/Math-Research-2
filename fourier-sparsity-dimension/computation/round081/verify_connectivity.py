"""Exact Walsh and GF(2) checks of the exterior connectivity constructions."""
import itertools
import json
import random
import sys
from fractions import Fraction
from functools import cache

SEED = 8104


def walsh(values):
    result = list(values)
    stride = 1
    while stride < len(result):
        for base in range(0, len(result), 2 * stride):
            for index in range(base, base + stride):
                left, right = result[index], result[index + stride]
                result[index], result[index + stride] = left + right, left - right
        stride *= 2
    return tuple(result)


def support(values):
    return tuple(index for index, value in enumerate(walsh(values)) if value)


def points(mask):
    while mask:
        value = mask & -mask
        yield value.bit_length() - 1
        mask ^= value


def reduce_vector(vector, pivots):
    for pivot, row in pivots.items():
        if vector & (1 << pivot):
            vector ^= row
    return vector


def gf_basis(vectors):
    pivots = {}
    for vector in vectors:
        value = reduce_vector(vector, pivots)
        if value:
            pivots[value.bit_length() - 1] = value
    return pivots


def rank(vectors):
    return len(gf_basis(vectors))


def components(labels):
    groups = list(range(len(labels)))
    pivots = {}

    def parent(i):
        while groups[i] != i:
            groups[i] = groups[groups[i]]
            i = groups[i]
        return i

    for index, label in enumerate(labels):
        relation = 1 << index
        while label:
            pivot = label.bit_length()
            if pivot not in pivots:
                pivots[pivot] = label, relation
                break
            row, expression = pivots[pivot]
            label ^= row
            relation ^= expression
        if not label:
            for key in points(relation):
                groups[parent(key)] = parent(index)
    partitions = {}
    for index, value in enumerate(labels):
        partitions.setdefault(parent(index), []).append(value)
    return tuple(tuple(group) for group in partitions.values())


@cache
def subspaces(dimension):
    all_points = (1 << (1 << dimension)) - 1
    current = {1}
    result = []
    for dim in range(dimension + 1):
        result.extend((mask, dim) for mask in sorted(current))
        following = set()
        for mask in current:
            outside = all_points ^ mask
            mask_points = tuple(points(mask))
            while outside:
                new = (outside & -outside).bit_length() - 1
                coset = sum(1 << (new ^ x) for x in mask_points)
                following.add(mask | coset)
                outside ^= coset
        current = following
    return tuple(result)


def check_products():
    audits = density_cases = 0
    for set_mask in range(1, 16):
        seed_set = tuple(points(set_mask))
        base_components = components(tuple(x | 4 for x in seed_set))
        for width in (1, 2, 3):
            dimension = 3 + width
            lifted = tuple(x | 4 | (w << 3) for x in seed_set for w in range(1 << width))
            first_disconnect = None
            for kernel, dim in subspaces(dimension):
                exterior = tuple(x for x in lifted if not (kernel >> x) & 1)
                count = len(components(exterior))
                if count > 1 and first_disconnect is None:
                    first_disconnect = dim
                audits += 1
            if len(seed_set) == 1:
                expected_disconnect = {1: 0, 2: 1, 3: None}[width]
            elif width == 1 and len(base_components) > 1:
                expected_disconnect = min(rank(group) for group in base_components)
            else:
                expected_disconnect = None
            assert first_disconnect == expected_disconnect

            ordinary = tuple(x | (w << 2) for x in seed_set for w in range(1 << width))
            ordinary_mask = sum(1 << x for x in ordinary if x)
            actual_density = max(
                Fraction((ordinary_mask & space).bit_count(), dim)
                for space, dim in subspaces(2 + width) if dim > 0
            )
            expected_density = max(
                Fraction((1 << width) * (set_mask & space).bit_count() - (0 in seed_set), width + dim)
                for space, dim in subspaces(2)
            )
            assert actual_density == expected_density
            density_cases += 1
    return audits, density_cases


def balanced(leaves):
    if len(leaves) == 1:
        return leaves[0]
    left, right = leaves[:len(leaves) // 2], leaves[len(leaves) // 2:]
    return balanced(left), balanced(right)


def build_parent(seed_truth, leaves_per_child, shared):
    seed_dimension = len(seed_truth).bit_length() - 1
    next_index = 1 if shared else 0
    first_seed = None

    def allocate():
        nonlocal next_index
        result = next_index
        next_index += 1
        return result

    def fill(shape, query_mask, answer_mask):
        nonlocal first_seed
        if shape == "C":
            return ("char", 0)
        if shape == "G":
            guard = allocate()
            start = next_index
            for _ in range(seed_dimension):
                allocate()
            if first_seed is None:
                first_seed = query_mask | (1 << guard), answer_mask, start
            return ("seed", guard, start)
        bit = allocate()
        left = fill(shape[0], query_mask | (1 << bit), answer_mask)
        right = fill(shape[1], query_mask | (1 << bit), answer_mask | (1 << bit))
        return ("query", bit, left, right)

    if shared:
        left = balanced(["C"] + ["G"] * (leaves_per_child - 1))
        right = balanced([("C", "G")] + ["G"] * (leaves_per_child - 1))
    else:
        left = right = ("G", "G")
    tree = fill((left, right), 0, 0)

    def evaluate(x):
        node = tree
        while node[0] == "query":
            node = node[2 + ((x >> node[1]) & 1)]
        sign = (-1) ** ((x >> node[1]) & 1)
        if node[0] == "char":
            return sign
        return sign * seed_truth[(x >> node[2]) & (len(seed_truth) - 1)]

    values = tuple(evaluate(x) for x in range(1 << next_index))
    restriction_mask, restriction_answers, start = first_seed
    assert all(
        value == seed_truth[(x >> start) & (len(seed_truth) - 1)]
        for x, value in enumerate(values) if x & restriction_mask == restriction_answers
    )
    return values, tree[1], next_index, restriction_mask.bit_count()


def random_subspace(basis, size, rng):
    pivots = {}
    while len(pivots) < size:
        vector = 0
        for key in points(rng.getrandbits(len(basis))):
            vector ^= basis[key]
        reduced = reduce_vector(vector, pivots)
        if reduced:
            pivots[reduced.bit_length() - 1] = reduced
    return pivots


def verify_parents():
    rng = random.Random(SEED)
    outcomes = []
    exhaustive_thresholds = 0
    checked_random_kernels = 0
    for shared, width, seed_truth in itertools.chain(
        ((shared, 2, tuple((-1) ** ((mask >> x) & 1) for x in range(4)))
         for shared in (False, True) for mask in range(16)),
        ((True, 4, (1,)),),
        ((False, 2, (1,)),),
    ):
        values, root_bit, dimension, restriction_codim = build_parent(seed_truth, width, shared)
        root = 1 << root_bit
        labels = support(values)
        source = support(seed_truth)
        affine_dim = rank(vector ^ source[0] for vector in source)
        base = tuple(x for x in labels if not x & root)
        assert set(labels) == set(base) | {x ^ root for x in base}
        component_list = components(base)
        basis = tuple(gf_basis(labels).values())
        if shared:
            expected_s = 4 * width**2 * len(source) + 6 * width - 2
            expected_d = (2 * width - 1) * (affine_dim + 2) + 2
            expected_threshold = None
            assert len(component_list) == 1
        else:
            expected_s = 16 * len(source)
            expected_d = 4 * affine_dim + 7
            expected_threshold = 2 * affine_dim + 3
            assert len(component_list) == 2
            assert [rank(group) for group in component_list] == [expected_threshold] * 2
            witness = gf_basis(component_list[0])
            exterior = tuple(label for label in labels if reduce_vector(label, witness))
            assert len(components(exterior)) == 2
        assert len(labels) == expected_s and len(basis) == expected_d
        for _ in range(128):
            limit = expected_d + 1 if shared else expected_threshold
            kernel = random_subspace(basis, rng.randrange(limit), rng)
            exterior = tuple(label for label in labels if reduce_vector(label, kernel))
            assert len(components(exterior)) <= 1
            checked_random_kernels += 1
        if not shared and len(seed_truth) == 1:
            for kernel, dim in subspaces(dimension):
                exterior = tuple(label for label in labels if not (kernel >> label) & 1)
                if dim < expected_threshold:
                    assert len(components(exterior)) == 1
                    exhaustive_thresholds += 1
            nonzero_mask = sum(1 << value for value in labels)
            density = max(
                Fraction((nonzero_mask & space).bit_count(), dim)
                for space, dim in subspaces(dimension) if dim > 0
            )
            assert density == Fraction(16, 7)
        outcomes.append({
            "shared": shared, "width": width, "s_seed": len(source),
            "s": len(labels), "d": len(basis), "first_disconnect": expected_threshold,
            "restriction_codimension": restriction_codim,
        })
    return outcomes, checked_random_kernels, exhaustive_thresholds


if __name__ == "__main__":
    audit_count, density_cases = check_products()
    results, random_kernels, exhaustive_thresholds = verify_parents()
    print(json.dumps({
        "python": sys.version.split()[0], "seed": SEED,
        "product_exact_kernels": audit_count, "product_exact_density_cases": density_cases,
        "parent_records": results, "parents_random_kernel_tests": random_kernels,
        "character_seed_all_kernels_below_threshold": exhaustive_thresholds,
    }, indent=2))
