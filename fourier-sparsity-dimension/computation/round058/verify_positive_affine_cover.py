#!/usr/bin/env python3
"""Exact checks for positive affine covers and a Boolean overlap family."""

import argparse
from collections import Counter, defaultdict
from hashlib import sha256
from itertools import combinations
import json
from pathlib import Path
import platform
import sys
import time


def parity(word):
    return word.bit_count() & 1


def character(word):
    return 1 - 2 * parity(word)


def span(vectors):
    result = {0}
    for vector in vectors:
        result |= {point ^ vector for point in result}
    return frozenset(result)


def rank(vectors):
    pivots = {}
    for vector in vectors:
        while vector:
            position = vector.bit_length() - 1
            if position not in pivots:
                pivots[position] = vector
                break
            vector ^= pivots[position]
    return len(pivots)


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


def affine_domains(dimension):
    size = 1 << dimension
    directions = {span(basis) for length in range(dimension + 1)
                  for basis in combinations(range(1, size), length)}
    return sorted({tuple(sorted(point ^ offset for point in direction))
                   for direction in directions for offset in range(size)})


def cover_census(total_weight):
    domains = affine_domains(2)
    assert len(domains) == 11
    normals = {
        domain: tuple(normal for normal in range(1, 4)
                      if all(not parity(normal & (point ^ domain[0])) for point in domain))
        for domain in domains
    }
    cover_number = 0
    overlapping = 0
    joint_ranks = Counter()

    def search(index, residual, active):
        nonlocal cover_number, overlapping
        if not any(residual):
            cover_number += 1
            dimension = rank(normal for domain in active for normal in normals[domain])
            assert dimension <= len(active) - 1
            joint_ranks[(len(active), dimension)] += 1
            overlapping += any(set(first) & set(second) for first, second in combinations(active, 2))
            return
        if index == len(domains):
            return
        domain = domains[index]
        for weight in range(1 + min(residual[point] for point in domain)):
            rest = list(residual)
            for point in domain:
                rest[point] -= weight
            search(index + 1, rest, active + [domain] if weight else active)

    search(0, [total_weight] * 4, [])
    return {
        "target_weight": total_weight,
        "exact_covers": cover_number,
        "overlapping_covers": overlapping,
        "profiles": {f"{domains},{normals}": count
                     for (domains, normals), count in sorted(joint_ranks.items())},
    }


def graph_truth(b):
    mask = (1 << b) - 1
    dimension = 3 * b + 2
    result = []
    for point in range(1 << dimension):
        first = point & mask
        second = (point >> b) & mask
        y1 = (point >> (2 * b)) & 1
        y2 = (point >> (2 * b + 1)) & 1
        z = point >> (2 * b + 2)
        result.append(character((y1 ^ parity(first & z)) * (y2 ^ parity(second & z))))
    return result


def expected_bands(b):
    size = 1 << b
    result = {0: (None, 1)}
    for data in range(size):
        result[1 | (data << 2)] = ((0, data), 1)
        result[2 | (data << 2)] = ((1, data), 1)
        result[3 | (data << 2)] = ((2, data), -1)
    return result


def verify_graph(b):
    queries = 1 << (2 * b)
    local_size = 1 << (b + 2)
    normal_size = 1 << b
    truth = graph_truth(b)
    bands = expected_bands(b)
    counts = Counter()
    active = defaultdict(set)
    for u in range(queries):
        transform = walsh(truth[u | (point << (2 * b))] for point in range(local_size))
        first = u & (normal_size - 1)
        second = u >> b
        expected = {0: 1, 1 | (first << 2): 1,
                    2 | (second << 2): 1, 3 | ((first ^ second) << 2): -1}
        actual = {label: value for label, value in enumerate(transform) if value}
        assert actual == {label: sign * local_size // 2 for label, sign in expected.items()}
        # Integer check of the total weighted coverage, with each weight 1/4.
        assert sum(value * value for value in actual.values()) == local_size ** 2
        counts.update(actual.keys())
        for label in actual:
            active[label].add(u)

    expected_spectrum = {0: len(truth) // 2}
    normal_rank_inputs = []
    for label, (description, sign) in bands.items():
        if description is None:
            assert len(active[label]) == queries
            continue
        direction, data = description
        assert counts[label] == normal_size
        normals = []
        for v in range(normal_size):
            if direction == 0:
                alpha = v
            elif direction == 1:
                alpha = v << b
            else:
                alpha = v | (v << b)
            normals.append(alpha)
            frequency = alpha | (label << (2 * b))
            expected_spectrum[frequency] = sign * character(v & data) * (len(truth) >> (b + 1))
        origin = min(active[label])
        assert all(not parity(normal & (origin ^ address))
                   for normal in normals for address in active[label])
        assert rank(normals) == b
        normal_rank_inputs.extend(normals)
    transform = walsh(truth)
    actual = {frequency: value for frequency, value in enumerate(transform) if value}
    assert actual == expected_spectrum
    assert len(actual) == 1 + 3 * 4 ** b and rank(actual) == 3 * b + 2
    assert rank(normal_rank_inputs) == 2 * b and len({frozenset(domain) for domain in active.values()}) == len(bands)
    return {
        "b": b, "truth_table_values": len(truth), "address_rows_checked": queries,
        "data_labels": len(bands), "distinct_affine_activity_domains": len(bands),
        "nonzero_label_repetitions": normal_size, "zero_label_repetitions": queries,
        "query_normal_rank": rank(normal_rank_inputs),
        "sparsity": len(actual), "linear_rank": rank(actual),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    start = time.perf_counter()
    covers = [cover_census(weight) for weight in (1, 2, 3)]
    graphs = [verify_graph(b) for b in range(1, 6)]
    result = {
        "schema": "a28-r056-positive-affine-cover-audit-v1",
        "geometric_covers_dimension_two": covers,
        "boolean_graphs": graphs,
        "python": sys.version,
        "platform": platform.platform(),
        "arithmetic": "exact integers; no randomness",
        "source_sha256": sha256(Path(__file__).read_bytes()).hexdigest(),
        "elapsed_seconds": time.perf_counter() - start,
    }
    args.output.write_text(json.dumps(result, indent=2) + "\n")
    print(f"verified {sum(row['exact_covers'] for row in covers)} affine covers and five Boolean sizes")


if __name__ == "__main__":
    main()
