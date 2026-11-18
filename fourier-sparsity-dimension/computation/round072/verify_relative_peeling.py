#!/usr/bin/env python3
"""Exact verification of relative affine peeling and the Boolean graph witness."""

from fractions import Fraction
import hashlib
import json
from pathlib import Path
import platform
import random
import sys
import time


def parity(value):
    return value.bit_count() & 1


def span(vectors):
    result = {0}
    for vector in vectors:
        if vector not in result:
            result |= {vector ^ entry for entry in result}
    return result


def rank(vectors):
    basis = {}
    for vector in vectors:
        while vector and vector.bit_length() in basis:
            vector ^= basis[vector.bit_length()]
        if vector:
            basis[vector.bit_length()] = vector
    return len(basis)


def walsh(values):
    values = list(values)
    width = 1
    while width < len(values):
        for start in range(0, len(values), 2 * width):
            for left in range(start, start + width):
                right = left + width
                first, second = values[left], values[right]
                values[left], values[right] = first + second, first - second
        width *= 2
    return {index: value for index, value in enumerate(values) if value}


def affine_membership(query, anchor, normal):
    return all(not parity((query ^ anchor) & direction) for direction in normal)


def verify_peeling(dimension, normals, anchors, old, margins):
    current = {index: normal for index, normal in enumerate(normals)}
    codimensions = [rank(normal) for normal in normals]
    assigned = [0] * len(normals)
    fixed = []
    total = rank(old | set().union(*normals)) - rank(old)
    telescoped = Fraction(0)
    expected_increments = Fraction(0)
    depth = slices_checked = 0
    while current:
        target = span(old | set().union(*current.values()))
        if target == old:
            break
        for index in tuple(current):
            alternatives = [normal for other, normal in current.items() if other != index]
            if span(old | set().union(*alternatives)) == target:
                del current[index]
        next_layer = {}
        extra = {}
        previous_assigned = assigned.copy()
        for index, normal in current.items():
            alternatives = [item for other, item in current.items() if other != index]
            residual = normal & span(old | set().union(*alternatives))
            partial_span = residual.copy()
            directions = []
            for direction in sorted(normal):
                if direction not in partial_span:
                    directions.append(direction)
                    partial_span |= {entry ^ direction for entry in partial_span}
            assert directions and partial_span == normal
            next_layer[index] = residual
            extra[index] = len(directions)
            fixed.extend((direction, parity(direction & anchors[index])) for direction in directions)
            assigned[index] += len(directions)
        reduced = span(old | set().union(*next_layer.values()))
        fixed_rank = rank(direction for direction, _ in fixed)
        assert fixed_rank == len(fixed)
        assert rank(reduced | {direction for direction, _ in fixed}) == rank(reduced) + fixed_rank
        affine_slice = [
            point for point in range(1 << dimension)
            if all(parity(direction & point) == value for direction, value in fixed)
        ]
        assert len(affine_slice) == 1 << (dimension - fixed_rank)
        for index, residual in next_layer.items():
            p = extra[index]
            h = previous_assigned[index]
            measure = Fraction(sum(affine_membership(query, anchors[index], normals[index]) for query in affine_slice), len(affine_slice))
            assert measure == Fraction(1, 1 << rank(residual))
            assert codimensions[index] == h + p + rank(residual)
            assert h >= depth
            weight = Fraction(2 * (h + margins[index]) * p + p * p, 1 << (h + p))
            telescoped += measure * weight
            weight_directions_only = Fraction(2 * h * p + p * p, 1 << (h + p))
            assert weight_directions_only <= (Fraction(9, 8) if not depth else Fraction(depth + 1, 1 << depth))
            expected_increments += ((margins[index] + h + p)**2 - (margins[index] + h)**2) / (1 << codimensions[index])
            slices_checked += 1
        current = next_layer
        depth += 1
    assert sum(assigned) == total
    assert telescoped == expected_increments
    assert telescoped == sum(
        ((margin + count)**2 - margin**2) / (1 << codimension)
        for margin, count, codimension in zip(margins, assigned, codimensions)
    )
    return depth, slices_checked, total


def check_affine_families(seed, cases):
    rng = random.Random(seed)
    levels = conditioned_checks = nonzero_old = total_dimensions = 0
    for _ in range(cases):
        dimension = rng.randrange(2, 7)
        ambient_size = 1 << dimension
        bands = rng.randrange(2, 10)
        normals = [span(rng.choices(range(ambient_size), k=rng.randrange(0, dimension + 1))) for _ in range(bands)]
        anchors = rng.choices(range(ambient_size), k=bands)
        old = span(rng.choices(range(ambient_size), k=rng.randrange(0, dimension + 1)))
        margins = [Fraction(rng.randrange(5), 4) for _ in range(bands)]
        depth, checks, total = verify_peeling(dimension, normals, anchors, old, margins)
        levels += depth
        conditioned_checks += checks
        nonzero_old += len(old) > 1
        total_dimensions += total
    return {
        "seed": seed,
        "families": cases,
        "families_with_nonzero_prepaid_space": nonzero_old,
        "levels_checked": levels,
        "layer_affine_measures_checked": conditioned_checks,
        "relative_direction_dimensions_summed": total_dimensions,
    }


def graph_witness(h):
    m = 1 << h
    mask = m - 1
    query_count = 1 << m

    def value(query, data):
        address = data & mask
        other = (data >> h) & mask
        switch = (data >> (2 * h)) & 1
        return 1 - 2 * (parity(address & other) ^ (((query >> address) & 1) * switch))

    actual = walsh(value(query, data) for data in range(2 * m * m) for query in range(query_count))
    expected = {}
    quotient_rows = {}
    amplitude = 1 << (m + h)
    for s in range(2):
        for q in range(m):
            for p in range(m):
                data = p | (q << h) | (s << (2 * h))
                quotient_rows[data] = q
                sign = 1 - 2 * parity(p & q)
                expected[data << m] = amplitude * sign
                expected[(data << m) | (1 << q)] = amplitude * sign * (1 - 2 * s)
    assert actual == expected
    assert len(actual) == 4 * m**2
    assert rank(actual) == m + 2 * h + 1
    actual_local_ranks = []
    for query in range(query_count):
        local = walsh(value(query, data) for data in range(2 * m * m))
        assert len(local) == m * m
        assert all(data >> (2 * h) == (query >> quotient_rows[data]) & 1 for data in local)
        actual_local_ranks.append(rank(local))
    assert max(actual_local_ranks) == 2 * h + 1
    anchors = [((data >> (2 * h)) & 1) << direction for data, direction in quotient_rows.items()]
    depth, _, directions = verify_peeling(m, [{0, 1 << direction} for direction in quotient_rows.values()], anchors, {0}, [Fraction(0)] * len(quotient_rows))
    assert depth == 1 and directions == m
    direction_marginal = Fraction(1, 2 * m)
    direction_squared_cost = len(quotient_rows) * direction_marginal**2 / 2
    assert direction_squared_cost == Fraction(1, 4)
    return {
        "h": h,
        "truth_table_size": query_count * 2 * m * m,
        "true_sparsity": len(actual),
        "full_rank": rank(actual),
        "maximum_real_local_representative_rank": max(actual_local_ranks),
        "real_local_sparsity": m * m,
        "peeled_coactive_bands": m,
        "mean_direction_squared_cost": str(direction_squared_cost),
    }


def main():
    started = time.monotonic()
    result = {
        "schema": "a28-r072-relative-peeling-v1",
        "affine_families": check_affine_families(20261118, 2000),
        "boolean_graphs": [graph_witness(h) for h in (1, 2, 3)],
        "python": sys.version,
        "platform": platform.platform(),
        "source_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "arithmetic": "exact integers and fractions",
        "elapsed_seconds": time.monotonic() - started,
    }
    Path(sys.argv[1]).write_text(json.dumps(result, ensure_ascii=False, indent=2) + "\n")


if __name__ == "__main__":
    main()
