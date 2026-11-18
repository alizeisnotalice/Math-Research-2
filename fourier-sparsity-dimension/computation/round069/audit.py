"""Exact Fourier and parity certificates for the affine independence threshold."""
import argparse
from collections import Counter, defaultdict
from itertools import combinations, combinations_with_replacement
import json
from math import factorial, prod
from pathlib import Path
import platform
import time


def walsh(values):
    values = list(values)
    stride = 1
    while stride < len(values):
        for base in range(0, len(values), stride * 2):
            for i in range(base, base + stride):
                a, b = values[i], values[i + stride]
                values[i], values[i + stride] = a + b, a - b
        stride *= 2
    return {i: value for i, value in enumerate(values) if value}


def character(product):
    return 1 - 2 * (product.bit_count() & 1)


def xor_total(items):
    point = 0
    for item in items:
        point ^= item
    return point


def evaluate(n, coefficients):
    return [sum(c * character(point & label) for label, c in coefficients.items())
            for point in range(1 << n)]


def relations(labels, affine):
    found = []
    for length in range(2 if affine else 1, len(labels) + 1, 2 if affine else 1):
        for subset in combinations(labels, length):
            if xor_total(subset) == 0:
                found.append(tuple(subset))
    return found


def multiply(left, right):
    coefficients = Counter()
    for x, a in left.items():
        for y, b in right.items():
            coefficients[x ^ y] += a * b
    return {x: a for x, a in coefficients.items() if a}


def polynomial_from_roots(roots):
    p = {0: 1}
    for root in roots:
        following = Counter()
        for i, c in p.items():
            following[i + 1] += c
            following[i] -= root * c
        p = {i: c for i, c in following.items() if c}
    return p


def exact_annihilator(spectrum, polynomial):
    zero = Counter()
    power = {0: 1}
    for degree in range(max(polynomial) + 1):
        for label, coefficient in power.items():
            zero[label] += coefficient * polynomial.get(degree, 0)
        power = multiply(power, spectrum)
    assert not any(zero.values())


def circuit(length):
    assert length >= 4 and length % 2 == 0
    n = length - 2
    labels = [0, *[1 << i for i in range(n)], (1 << n) - 1]
    coefficients = {label: -1 if label == 0 else 1 for label in labels}
    affine_relations = relations(labels, affine=True)
    assert affine_relations == [tuple(labels)]
    values = evaluate(n, coefficients)
    size = 1 << n
    assert walsh(values) == {label: value * size for label, value in coefficients.items()}
    histogram = Counter(values)
    assert sorted(histogram) == list(range(2 - length, length - 1, 4))
    assert len(histogram) == length // 2

    # Affine sign lift attaches one independent sign to all characters.
    # Every sign pattern with an odd number of negative adjusted signs occurs once.
    augmented = {label | size: value for label, value in coefficients.items()}
    lifted = evaluate(n + 1, augmented)
    assert lifted == [*values, *[-value for value in values]]
    adjusted_patterns = {tuple(c * character(x & label) for label, c in augmented.items())
                         for x in range(size * 2)}
    assert len(adjusted_patterns) == size * 2
    assert all(pattern.count(-1) % 2 == 1 for pattern in adjusted_patterns)

    # Independently verify the only degree-top collisions: complementary subsets.
    degree = length // 2
    classes = defaultdict(list)
    for subset in combinations(augmented, degree):
        classes[xor_total(subset)].append(subset)
    assert all(len(rows) == 2 for rows in classes.values())
    products = {}
    for frequency, (first, second) in classes.items():
        assert set(first).isdisjoint(second)
        first_coeff = factorial(degree) * prod(augmented[label] for label in first)
        second_coeff = factorial(degree) * prod(augmented[label] for label in second)
        assert first_coeff == -second_coeff
        products[frequency] = first_coeff + second_coeff
    assert not any(products.values())
    polynomial = polynomial_from_roots(histogram)
    assert all(deg % 2 == degree % 2 for deg in polynomial)
    exact_annihilator(augmented, polynomial)
    return {"length": length, "physical_dimension": n, "labels": labels,
            "coefficients_before_normalization": coefficients,
            "value_histogram_before_normalization": dict(sorted(histogram.items())),
            "affine_relation_sizes": [len(row) for row in affine_relations],
            "affine_odd_sign_patterns": len(adjusted_patterns),
            "top_collision_classes": len(classes), "annihilator_coefficients": polynomial,
            "highest_independence_order_of_counterexample": length - 1}


def boolean_rows(example):
    labels = example["coefficients_before_normalization"]
    n = example["physical_dimension"]
    address_size = 1 << n
    g_times_eight = evaluate(n, labels)
    table = [1 - 2 * (data < 2 - value // 4)
             for data in range(4) for value in g_times_eight]
    local_profiles = Counter()
    for point in range(address_size):
        row = [table[point + address_size * j] for j in range(4)]
        assert 2 * sum(row) == g_times_eight[point]
        local = walsh(row)
        local_profiles[(g_times_eight[point], tuple(local))] += 1
    full = walsh(table)
    zero_band = {freq: value for freq, value in full.items() if freq < address_size}
    assert zero_band == {label: address_size * val // 2 for label, val in labels.items()}
    return {"truth_entries": len(table), "zero_band_normalized_denominator": 8,
            "zero_band_numerators": labels,
            "actual_data_projection": sorted({freq >> n for freq in full}),
            "full_fourier_sparsity": len(full),
            "slice_profiles": {f"{m},{support}": count for (m, support), count in local_profiles.items()}}


def bounded_weights():
    records = Counter()
    for weights in combinations_with_replacement(range(1, 5), 10):
        total = sum(weights)
        levels = [set(), set()]
        for mask in range(1 << 10):
            negative = sum(w for i, w in enumerate(weights) if (mask >> i) & 1)
            levels[mask.bit_count() & 1].add(total - 2 * negative)
        for parity, values in enumerate(levels):
            records["cases"] += 1
            if len(values) > 5:
                continue
            records["at_most_five_values"] += 1
            assert parity == 1 and len(set(weights)) == 1 and len(values) == 5
    return records


def asymmetric_example():
    n = 8
    labels = [*[1 << i for i in range(n)], (1 << n) - 1]
    values = evaluate(n, dict.fromkeys(labels, 1))
    assert not relations(labels, affine=True)
    assert relations(labels, affine=False) == [tuple(labels)]
    assert sorted(set(values)) == [-7, -3, 1, 5, 9]
    assert walsh(values) == dict.fromkeys(labels, 1 << n)
    return {"dimension": n, "labels": labels, "affine_relation_sizes": [],
            "linear_relation_sizes": [len(labels)], "values": sorted(set(values)),
            "histogram": dict(sorted(Counter(values).items()))}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    start = time.perf_counter()
    examples = [circuit(size) for size in (6, 8, 10, 12)]
    result = {"schema": "a28-r069-affine-wise-circuits-v1",
              "python": platform.python_version(), "circuits": examples,
              "genuine_boolean_rows": boolean_rows(examples[2]),
              "weighted_10_circuit_census": bounded_weights(),
              "asymmetric_five_value_counterexample": asymmetric_example(),
              "elapsed_seconds": time.perf_counter() - start}
    args.output.write_text(json.dumps(result, indent=2) + "\n")
    print([(r["length"], r["top_collision_classes"]) for r in examples])
    print(result["weighted_10_circuit_census"], result["elapsed_seconds"])
