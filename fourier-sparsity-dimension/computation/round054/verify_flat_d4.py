#!/usr/bin/env python3
"""Independent exact tests of the four-dimensional flat-spectrum argument."""

import argparse
from collections import Counter
from hashlib import sha256
from itertools import combinations, product
import json
from pathlib import Path
import platform
import sys
import time


def parity(word):
    return word.bit_count() & 1


def character(word):
    return 1 - 2 * parity(word)


def walsh(values):
    output = list(values)
    width = 1
    while width < len(output):
        for start in range(0, len(output), 2 * width):
            for index in range(start, start + width):
                x, y = output[index], output[index + width]
                output[index], output[index + width] = x + y, x - y
        width *= 2
    return output


def signed_rows(zeros):
    active = [index for index in range(8) if index not in zeros]
    for signs in product((-1, 1), repeat=len(active)):
        result = [0] * 8
        for index, sign in zip(active, signs):
            result[index] = sign
        yield result


def row_census():
    lines = {frozenset((a, b, a ^ b)) for a, b in combinations(range(1, 8), 2)}
    odd_without_line = 0
    for row in signed_rows({0}):
        coefficients = list(map(abs, walsh(row)))
        large = {label for label in range(1, 8) if coefficients[label] == 3}
        if max(coefficients[1:]) <= 3 and not any(line <= large for line in lines):
            assert coefficients[0] == 7
            odd_without_line += 1
    even_candidates = 0
    for row in signed_rows({0, 1}):
        coefficients = list(map(abs, walsh(row)))
        if max(coefficients[1:]) <= 2:
            assert coefficients[0] == 6
            assert all(coefficients[label] == 2 for label in (2, 4, 6))
            even_candidates += 1
    three_zero_candidates = sum(max(map(abs, walsh(row)[1:])) <= 1 for row in signed_rows({0, 1, 2}))
    assert (odd_without_line, even_candidates, three_zero_candidates) == (2, 2, 0)
    return {"one_zero_sign_patterns": 128, "one_zero_no_line_candidates": odd_without_line,
            "two_zero_sign_patterns": 64, "two_zero_small_label_candidates": even_candidates,
            "three_zero_sign_patterns": 32, "three_zero_small_label_candidates": three_zero_candidates}


def core_census():
    core = [point for point in range(16)
            if ((point & 1) * ((point >> 1) & 1)) ^ (((point >> 2) & 1) * ((point >> 3) & 1))]
    accepted = 0
    for signs in product((-1, 1), repeat=6):
        values = [0] * 16
        for point, sign in zip(core, signs):
            values[point] = sign
        transform = walsh(values)
        if set(transform) <= {-4, 0, 4}:
            dual = [point for point, value in enumerate(transform) if value]
            assert len(dual) == 6 and len({a ^ b for a, b in combinations(dual, 2)}) == 15
            indicator_transform = walsh(int(point in dual) for point in range(16))
            assert indicator_transform[0] == 6 and set(map(abs, indicator_transform[1:])) == {2}
            accepted += 1
    assert accepted == 32
    return {"core_sign_patterns": 64, "accepted_phases": accepted, "dual_indicator_norm": "9/4"}


def addressed_values(norm):
    address_bits = (norm // 2).bit_length() - 1
    dimension = address_bits + 3 * (norm // 2)
    values = []
    for point in range(1 << dimension):
        cell = point & ((1 << address_bits) - 1)
        private = point >> (address_bits + 3 * cell)
        values.append(character((private & 1) ^ (((private >> 1) & 1) * ((private >> 2) & 1))))
    return values, address_bits


def paired_line_values():
    groups = {point: group for group, line in enumerate(((0, 1), (2, 4), (3, 7), (5, 6)))
              for point in line}
    return [character((point >> (3 + groups[point & 7])) & 1) for point in range(128)]


def independent_basis(vectors):
    pivots = {}
    result = []
    for vector in vectors:
        residue = vector
        while residue:
            high = residue.bit_length() - 1
            if high not in pivots:
                pivots[high] = residue
                result.append(vector)
                break
            residue ^= pivots[high]
    return result


def span(basis):
    code = [0]
    for direction in basis:
        code += [word ^ direction for word in code]
    return code


def rref_spaces(dimension, rank):
    for pivot_columns in combinations(range(dimension), rank):
        remaining = set(range(dimension)) - set(pivot_columns)
        writable = [(row, column) for row, pivot in enumerate(pivot_columns)
                    for column in sorted(remaining) if column > pivot]
        for state in range(1 << len(writable)):
            basis = [1 << column for column in pivot_columns]
            for bit, (row, column) in enumerate(writable):
                if state >> bit & 1:
                    basis[row] |= 1 << column
            yield basis


def evaluation_masks(values, norm):
    transform = walsh(values)
    assert all(not value or abs(value) * norm == len(values) for value in transform)
    support = [index for index, value in enumerate(transform) if value]
    signs = [1 if transform[index] > 0 else -1 for index in support]
    assert len(support) == norm ** 2
    dimension = len(values).bit_length() - 1
    masks = [sum(((frequency >> coordinate) & 1) << bit for bit, frequency in enumerate(support))
             for coordinate in range(dimension)]
    constant = (1 << len(support)) - 1
    return support, signs, masks, constant


def verify_full_model(name, values):
    norm = 4
    support, signs, masks, constant = evaluation_masks(values, norm)
    basis = independent_basis([*masks, constant])
    code = span(basis)
    assert len(basis) == 7
    for word in code:
        # All affine evaluation words, including constant phase, have Walsh
        # transform in {0, ±N}; degree≤1 moments in this test must be even.
        active = [position for position in range(len(support)) if word >> position & 1]
        half = [sum(signs[position] * character(point & support[position]) for position in active)
                for point in range(len(values))]
        assert set(half) <= {-norm, 0, norm}
        assert len(active) % 2 == 0
        assert all(sum((support[position] >> coordinate) & 1 for position in active) % 2 == 0
                   for coordinate in range(len(masks)))
    counts = Counter()
    equality_minimum_words = Counter()
    lowest_by_number_minimum = {}
    for plane in rref_spaces(len(basis), 4):
        representatives = [code[vector] for vector in plane]
        union = (representatives[0] | representatives[1] | representatives[2] | representatives[3]).bit_count()
        assert union >= 3 * norm
        counts[union] += 1
        minima = sum(word.bit_count() == norm for word in span(representatives))
        lowest_by_number_minimum[minima] = min(union, lowest_by_number_minimum.get(minima, union))
        if union == 3 * norm:
            equality_minimum_words[minima] += 1
    assert sum(counts.values()) == 11811 and min(counts) == 12
    return {"name": name, "N": norm, "dimension": len(basis), "four_dimensional_spaces": sum(counts.values()),
            "all_support_counts": dict(sorted(counts.items())),
            "sharp_minimum_word_counts": dict(sorted(equality_minimum_words.items())),
            "lowest_support_by_number_of_minimum_words": dict(sorted(lowest_by_number_minimum.items())),
            "affine_words_checked": len(code)}


def check_construction(norm):
    values, address_bits = addressed_values(norm)
    support, _, masks, constant = evaluation_masks(values, norm)
    selected = [masks[index] for index in (address_bits, address_bits + 1, address_bits + 2, address_bits + 4)]
    assert len(independent_basis(selected)) == 4
    union = selected[0] | selected[1] | selected[2] | selected[3]
    assert union.bit_count() == 3 * norm
    assert len(independent_basis([*masks, constant])) == len(masks)
    return {"N": norm, "truth_table_size": len(values), "sparsity": len(support),
            "witness_support_size": union.bit_count()}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    started = time.perf_counter()
    result = {
        "schema": "a28-r055-flat-four-dimensional-audit-v1",
        "parameter_row_census": row_census(),
        "sidon_core_census": core_census(),
        "full_small_models": [
            verify_full_model("two addressed bent phases", addressed_values(4)[0]),
            verify_full_model("four nonparallel address lines", paired_line_values()),
        ],
        "equality_witnesses": [check_construction(norm) for norm in (4, 8)],
        "python": sys.version,
        "platform": platform.platform(),
        "source_sha256": sha256(Path(__file__).read_bytes()).hexdigest(),
        "elapsed_seconds": time.perf_counter() - started,
    }
    args.output.write_text(json.dumps(result, indent=2) + "\n")
    print("verified the local Walsh cases and 23,622 four-dimensional subspaces")


if __name__ == "__main__":
    main()
