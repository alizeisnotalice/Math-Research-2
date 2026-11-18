#!/usr/bin/env python3
"""Exact spectra and weighted-rank checks for shared-label gluing examples."""

import argparse
from collections import Counter, defaultdict
from fractions import Fraction
from hashlib import sha256
from itertools import combinations
import json
from pathlib import Path
import platform
import sys
import time


def parity(word):
    return word.bit_count() & 1


def prefix_value(prefix):
    return sum(int(bit) << index for index, bit in enumerate(prefix))


def word_at(value, length):
    return "".join(str((value >> index) & 1) for index in range(length))


def rank(points):
    pivots = {}
    for point in points:
        while point:
            pivot = point.bit_length() - 1
            if pivot not in pivots:
                pivots[pivot] = point
                break
            point ^= pivots[pivot]
    return len(pivots)


def affrank(points):
    origin = min(points)
    return rank(point ^ origin for point in points)


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


def square(first, second, third):
    fourth = first ^ second ^ third
    assert len({first, second, third, fourth}) == 4
    return {first: 1, second: 1, third: 1, fourth: -1}


def central_family(subdivision_depth):
    tags = {}

    def fresh(name):
        tag = 1 << len(tags)
        tags[name] = tag
        return tag

    beta0, beta1, beta2 = (fresh(f"beta{i}") for i in range(3))
    beta3 = beta0 ^ beta1 ^ beta2
    gamma12, gamma23, gamma31 = (fresh(f"gamma{name}") for name in ("12", "23", "31"))
    blocks = {
        "0000": square(beta0, beta1, beta2),
        "010": square(beta1, gamma12, gamma31),
        "10": square(beta2, gamma12, gamma23),
        # The shared beta3 is the negative fourth point on both incident blocks.
        "11": square(gamma23, gamma31, beta3 ^ gamma23 ^ gamma31),
        "000100": square(beta0, fresh("small_private_0"), fresh("small_private_1")),
    }
    for index in range(1 << subdivision_depth):
        prefix = "011" + word_at(index, subdivision_depth)
        blocks[prefix] = square(beta1, fresh(f"refined_{index}_0"), fresh(f"refined_{index}_1"))
    for prefix in ("000101", "00011", "001"):
        blocks[prefix] = square(*(fresh(f"filler_{prefix}_{coordinate}") for coordinate in range(3)))
    named_shared = dict(zip(("beta0", "beta1", "beta2", "beta3", "gamma12", "gamma23", "gamma31"),
                            (beta0, beta1, beta2, beta3, gamma12, gamma23, gamma31)))
    return blocks, named_shared, tags


def subcubic_graph(address_bits, edges):
    count = 1 << address_bits
    tags = {}
    incident = [[] for _ in range(count)]

    def fresh(name):
        tag = 1 << len(tags)
        tags[name] = tag
        return tag

    named_shared = {}
    for first, second in edges:
        name = f"edge_{first}_{second}"
        tag = fresh(name)
        named_shared[name] = tag
        incident[first].append(tag)
        incident[second].append(tag)
    blocks = {}
    for vertex, labels in enumerate(incident):
        assert len(labels) <= 3
        while len(labels) < 3:
            labels.append(fresh(f"private_{vertex}_{len(labels)}"))
        blocks[word_at(vertex, address_bits)] = square(*labels)
    return blocks, named_shared, tags


def components(blocks, incidence):
    adjacent = {
        prefix: {other for label in local for other in incidence[label] if other != prefix}
        for prefix, local in blocks.items()
    }
    pending = set(blocks)
    result = []
    while pending:
        stack = [min(pending)]
        group = set()
        while stack:
            prefix = stack.pop()
            if prefix in pending:
                pending.remove(prefix)
                group.add(prefix)
                stack.extend(adjacent[prefix])
        if len(group) > 1:
            shared_all = set.intersection(*(set(blocks[prefix]) for prefix in group))
            old_core_passes = bool(shared_all) and all(
                set(blocks[a]) & set(blocks[b]) == shared_all for a, b in combinations(group, 2)
            )
            result.append({"prefixes": sorted(group), "common_labels": sorted(shared_all),
                           "old_common_core_passes": old_core_passes})
    return result


def bound_rows(blocks, incidence, names, data_dimension):
    output = {}
    lifting = 1 << data_dimension
    for prefix, local in blocks.items():
        mass = Fraction(1, 1 << len(prefix))
        shared = {beta: mass / sum(Fraction(1, 1 << len(other)) for other in incidence[beta])
                  for beta in local if len(incidence[beta]) > 1}
        private = 4 - len(shared)
        local_effective = private + sum(weight ** 2 for weight in shared.values())
        choices = []
        for size in range(len(shared) + 1):
            for subset in combinations(shared, size):
                if rank(label | lifting for label in subset) == size:
                    refund = sum((1 - shared[label] for label in subset), Fraction(0))
                    choices.append((refund, subset))
        refund, subset = max(choices)
        margin = 4 * local_effective - (4 - refund) ** 2
        assert 0 <= refund <= 4 and margin >= 0
        output[prefix] = {
            "mass": str(mass),
            "weights": {names[beta]: str(weight) for beta, weight in shared.items()},
            "private_labels": private,
            "effective_size": str(local_effective),
            "max_independent_refund": str(refund),
            "maximizing_independent_labels": [names[label] for label in subset],
            "squared_margin": str(margin),
            "old_R41_rank_condition": rank(label | lifting for label in shared) == len(shared),
        }
    return output


def exact_spectrum(blocks, address_bits, details):
    address_size = 1 << address_bits
    active = []
    for point in range(address_size):
        sites = [prefix for prefix in blocks
                 if point & ((1 << len(prefix)) - 1) == prefix_value(prefix)]
        assert len(sites) == 1
        active.append(blocks[sites[0]])
    incidence = defaultdict(list)
    directly_assembled = {}
    for prefix, local in blocks.items():
        assert affrank(local) == 2 and sum(local.values()) == 2
        address = prefix_value(prefix)
        width = 1 << len(prefix)
        for label, sign in local.items():
            incidence[label].append(prefix)
            spectrum = directly_assembled.setdefault(label, Counter())
            for frequency in range(width):
                coefficient = sign * (1 - 2 * parity(frequency & address)) * (address_size // width)
                spectrum[frequency] += coefficient

    result = {}
    for label in incidence:
        spectrum = walsh(local.get(label, 0) for local in active)
        direct = directly_assembled[label]
        assert all(coefficient == direct[frequency] for frequency, coefficient in enumerate(spectrum))
        result[label] = {frequency: coefficient for frequency, coefficient in enumerate(spectrum)
                         if coefficient}
        assert result[label]
    assert set(result) == set(incidence)
    if details:
        parent = [(frequency | (label << address_bits), coefficient)
                  for label, band in result.items() for frequency, coefficient in band.items()]
        square_coefficients = Counter()
        for i, (point, coefficient) in enumerate(parent):
            for j in range(i):
                other_point, other_coefficient = parent[j]
                square_coefficients[point ^ other_point] += coefficient * other_coefficient
        assert all(value == 0 for value in square_coefficients.values())
        assert sum(coefficient ** 2 for _, coefficient in parent) == (2 * address_size) ** 2
    return result, incidence


def verify(blocks, named_shared, tags, with_details):
    addresses = max(map(len, blocks))
    spectrum, incidence = exact_spectrum(blocks, addresses, with_details)
    inverse_shared = {label: name for name, label in named_shared.items()}
    actual_shared = {label for label, sites in incidence.items() if len(sites) >= 2}
    assert actual_shared == set(named_shared.values())
    rows = bound_rows(blocks, incidence, inverse_shared, len(tags))
    points = {frequency | (label << addresses) for label, band in spectrum.items()
              for frequency in band}
    result = {
        "address_dimension": addresses,
        "data_dimension": rank(spectrum),
        "blocks": len(blocks),
        "shared_label_multiplicities": {name: len(incidence[beta]) for name, beta in named_shared.items()},
        "shared_band_sizes": {name: len(spectrum[beta]) for name, beta in named_shared.items()},
        "private_sparsity": sum(len(band) for label, band in spectrum.items() if label not in actual_shared),
        "sparsity": len(points),
        "linear_dimension": rank(points),
        "affine_dimension": affrank(points),
        "minimum_local_squared_margin": str(min(Fraction(row["squared_margin"]) for row in rows.values())),
        "old_R41_rank_failures": sorted(prefix for prefix, row in rows.items()
                                       if not row["old_R41_rank_condition"]),
        "old_common_core_passes": all(group["old_common_core_passes"]
                                     for group in components(blocks, incidence)),
        "center_condition": rows.get("0000"),
        "full_square_convolution_checked": with_details,
    }
    if with_details:
        result["local_condition_rows"] = rows
        result["connected_overlaps"] = components(blocks, incidence)
        result["input_blocks"] = blocks
        result["independent_data_labels"] = tags
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    started = time.perf_counter()
    family = []
    for depth in range(6):
        copies = 1 << depth
        blocks, labels, tags = central_family(depth)
        result = verify(blocks, labels, tags, depth == 0)
        assert result["shared_label_multiplicities"]["beta1"] == copies + 2
        assert result["sparsity"] == 754 + 24 * copies ** 2
        linear = 17 + 2 * copies + max(6, 3 + depth)
        assert result["linear_dimension"] == linear and result["affine_dimension"] == linear - 1
        assert result["old_R41_rank_failures"] == ["0000"] and not result["old_common_core_passes"]
        assert result["center_condition"]["max_independent_refund"] == "12/5"
        assert result["center_condition"]["squared_margin"] == "12/25"
        family.append({"subdivision_depth": depth, **result})

    graphs = {}
    cases = {"four_cycle": ((0, 1), (1, 2), (2, 3), (0, 3)),
             "complete_four": tuple(combinations(range(4), 2))}
    for name, edges in cases.items():
        result = verify(*subcubic_graph(2, edges), True)
        assert result["sparsity"] == 64 - 6 * len(edges)
        assert result["linear_dimension"] == 14 - len(edges)
        assert result["affine_dimension"] == 13 - len(edges)
        assert not result["old_common_core_passes"]
        graphs[name] = result

    output = {
        "schema": "a28-r050-weighted-overlap-examples-v1",
        "central_family": family,
        "graphs": graphs,
        "python": sys.version,
        "platform": platform.platform(),
        "arithmetic": "exact integers and fractions; no randomness",
        "source_sha256": sha256(Path(__file__).read_bytes()).hexdigest(),
        "elapsed_seconds": time.perf_counter() - started,
    }
    args.output.write_text(json.dumps(output, indent=2) + "\n")
    print(f"verified {len(family)} weighted families and {len(graphs)} graph examples")


if __name__ == "__main__":
    main()
