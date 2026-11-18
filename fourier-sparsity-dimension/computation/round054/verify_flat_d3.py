#!/usr/bin/env python3
"""Exact Walsh rows and small three-dimensional evaluation code checks."""

import argparse
from collections import Counter
from hashlib import sha256
from itertools import combinations
import json
from pathlib import Path
import platform
import sys
import time


def parity(value):
    return value.bit_count() & 1


def walsh(values):
    values = list(values)
    width = 1
    while width < len(values):
        for start in range(0, len(values), 2 * width):
            for index in range(start, start + width):
                left, right = values[index], values[index + width]
                values[index], values[index + width] = left + right, left - right
        width *= 2
    return values


def row_census():
    profiles = Counter()
    for mask in range(256):
        row = [1 - 2 * ((mask >> index) & 1) for index in range(8)]
        coefficients_in_quarters = sorted(abs(value) // 2 for value in walsh(row) if value)
        profiles[tuple(coefficients_in_quarters)] += 1
    assert profiles == {(4,): 16, (2, 2, 2, 2): 112, (1, 1, 1, 1, 1, 1, 1, 3): 128}
    return {",".join(map(str, profile)): count for profile, count in sorted(profiles.items())}


def independent_basis(vectors):
    pivots = {}
    actual = []
    for vector in vectors:
        reduced = vector
        while reduced:
            pivot = reduced.bit_length() - 1
            if pivot not in pivots:
                pivots[pivot] = reduced
                actual.append(vector)
                break
            reduced ^= pivots[pivot]
    return actual


def affine_rank(points):
    origin = min(points)
    return len(independent_basis(point ^ origin for point in points))


def spaces_of_rank_three(dimension):
    if dimension < 3:
        return
    for pivots in combinations(range(dimension), 3):
        unused = set(range(dimension)) - set(pivots)
        entries = [(row, column) for row, pivot in enumerate(pivots)
                   for column in sorted(unused) if column > pivot]
        for mask in range(1 << len(entries)):
            basis = [1 << coordinate for coordinate in pivots]
            for index, (row, column) in enumerate(entries):
                if (mask >> index) & 1:
                    basis[row] |= 1 << column
            yield basis


def span(basis):
    result = [0]
    for word in basis:
        result += [old ^ word for old in result]
    return result


def sharp_values(norm):
    address_bits = (norm // 2).bit_length() - 1
    address_mask = (1 << address_bits) - 1
    dimension = address_bits + 3 * (norm // 2)
    values = []
    for point in range(1 << dimension):
        cell = point & address_mask
        private = point >> (address_bits + 3 * cell)
        phase = (private & 1) ^ (((private >> 1) & 1) * ((private >> 2) & 1))
        values.append(1 - 2 * phase)
    return values


def check_affine_phase(values, parent):
    """Nonzero constant phase is included in the direct Fourier check."""
    size = len(values)
    for shift, constant in ((0, 1), (size - 1, 1), (size // 2, 0)):
        transform = walsh((1 - 2 * constant) * values[point ^ shift] for point in range(size))
        assert transform == [
            coefficient * (1 - 2 * (constant ^ parity(frequency & shift)))
            for frequency, coefficient in enumerate(parent)
        ]


def equality_checks(support, signed, evaluations, norm, dimension):
    labels = {}
    for bit, frequency in enumerate(support):
        label = sum(((word >> bit) & 1) << index for index, word in enumerate(evaluations))
        if label:
            labels.setdefault(label, []).append(frequency)
    assert len(labels) == 4 and all(len(points) == norm // 2 for points in labels.values())
    assert sorted(word.bit_count() for word in span(evaluations)[1:]) == [norm] * 6 + [2 * norm]
    directions = {frozenset(point ^ points[0] for point in points) for points in labels.values()}
    assert len(directions) == 1
    assert affine_rank(sum(labels.values(), [])) == (2 * norm).bit_length() - 1

    physical_supports = []
    for points in labels.values():
        partial = [sum(signed[frequency] * (1 - 2 * parity(frequency & point)) for frequency in points)
                   for point in range(1 << dimension)]
        assert set(partial) <= {-norm // 2, 0, norm // 2}
        physical_supports.append({point for point, value in enumerate(partial) if value})
    assert all(domain == physical_supports[0] for domain in physical_supports)
    domain = physical_supports[0]
    assert len(domain) * norm == 2 * (1 << dimension)
    assert affine_rank(domain) == len(domain).bit_length() - 1
    x0 = min(domain)
    constants = []
    for points in labels.values():
        phases = {signed[frequency] * (1 - 2 * parity(frequency & x0)) for frequency in points}
        assert len(phases) == 1
        constants.append(phases.pop())
    assert sum(sign < 0 for sign in constants) % 2 == 1


def audit_model(norm, values, complete):
    dimension = len(values).bit_length() - 1
    parent = walsh(values)
    assert all(not value or abs(value) * norm == len(values) for value in parent)
    support = [frequency for frequency, value in enumerate(parent) if value]
    assert len(support) == norm ** 2
    signed = {frequency: 1 if parent[frequency] > 0 else -1 for frequency in support}
    check_affine_phase(values, parent)
    columns = [point | (1 << dimension) for point in support]
    evaluation_axes = independent_basis(
        sum(((column >> coordinate) & 1) << bit for bit, column in enumerate(columns))
        for coordinate in range(dimension + 1)
    )
    code = span(evaluation_axes) if complete else None
    support_counts = Counter()
    direct_geometry_checks = 0
    if complete:
        for basis in spaces_of_rank_three(len(evaluation_axes)):
            evaluations = [code[row] for row in basis]
            union = (evaluations[0] | evaluations[1] | evaluations[2]).bit_count()
            assert union >= 2 * norm
            support_counts[union] += 1
            if union == 2 * norm:
                equality_checks(support, signed, evaluations, norm, dimension)
                direct_geometry_checks += 1
        if norm >= 2:
            assert min(support_counts) == 2 * norm
    else:
        address_bits = (norm // 2).bit_length() - 1
        first_cell = [
            sum(((point >> coordinate) & 1) << index for index, point in enumerate(support))
            for coordinate in range(address_bits, address_bits + 3)
        ]
        assert (first_cell[0] | first_cell[1] | first_cell[2]).bit_count() == 2 * norm
        equality_checks(support, signed, first_cell, norm, dimension)
        direct_geometry_checks = 1
    return {
        "N": norm,
        "physical_dimension": dimension,
        "sparsity": len(support),
        "evaluation_rank": len(evaluation_axes),
        "complete_three_dimensional_enumeration": complete,
        "three_dimensional_subspaces_checked": sum(support_counts.values()),
        "support_size_histogram": dict(sorted(support_counts.items())),
        "equality_geometry_and_phase_checks": direct_geometry_checks,
        "phase_restriction_checks": 3,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    started = time.perf_counter()
    result = {
        "schema": "a28-flat-three-dimensional-audit-v1",
        "all_parameter_row_profiles_in_quarters": row_census(),
        "models": [audit_model(1, [1, -1], True)]
        + [audit_model(norm, sharp_values(norm), norm <= 4) for norm in (2, 4, 8)],
        "python": sys.version,
        "platform": platform.platform(),
        "source_sha256": sha256(Path(__file__).read_bytes()).hexdigest(),
        "elapsed_seconds": time.perf_counter() - started,
    }
    args.output.write_text(json.dumps(result, indent=2) + "\n")
    print("verified all three-parameter Boolean rows and low-N evaluation models")


if __name__ == "__main__":
    main()
