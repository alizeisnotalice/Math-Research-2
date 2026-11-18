"""Exact Fourier, density certificate, and greedy-prefix audits for A28."""

import json
import sys
from collections import Counter, defaultdict
from fractions import Fraction
from itertools import combinations, product
from math import isqrt


def echelon(vectors):
    pivots = {}
    for index, value in enumerate(vectors):
        remainder, used = reduce_mod(value, pivots)
        if remainder:
            pivots[remainder.bit_length()] = (remainder, used ^ (1 << index))
    return pivots


def reduce_mod(value, pivots):
    used = 0
    for bit in sorted(pivots, reverse=True):
        if value & (1 << (bit - 1)):
            row, mask = pivots[bit]
            value ^= row
            used ^= mask
    return value, used


def rank(vectors):
    return len(echelon(vectors))


def span(vectors):
    output = [0]
    for vector in vectors:
        output += [vector ^ previous for previous in output]
    return output


def signed_walsh(table):
    output = list(table)
    width = 1
    while width < len(output):
        for start in range(0, len(output), 2 * width):
            for left in range(start, start + width):
                right = left + width
                a, b = output[left], output[right]
                output[left], output[right] = a + b, a - b
        width *= 2
    return {index: value for index, value in enumerate(output) if value}


def audit_walsh(spectrum, dimension, exact_value):
    table = [exact_value(x) for x in range(1 << dimension)]
    actual = signed_walsh(table)
    scale = (1 << dimension) // isqrt(len(spectrum))
    assert actual == {index: value * scale for index, value in spectrum.items()}


def fiber_weights(spectrum, queries):
    quotient = echelon(queries)
    assert len(quotient) == len(queries)
    fibers = defaultdict(list)
    for frequency, sign in spectrum.items():
        remainder, mask = reduce_mod(frequency, quotient)
        fibers[remainder].append((mask, sign))
    weights = {}
    for label, entries in fibers.items():
        anchor, first_sign = entries[0]
        difference_rank = rank(mask ^ anchor for mask, _ in entries)
        assert len(entries) == 1 << difference_rank
        assert rank(
            ((mask ^ anchor) << 1) | int(sign != first_sign)
            for mask, sign in entries
        ) == difference_rank
        weights[label] = Fraction(1, len(entries))
    if len(queries) <= 5:
        for label, entries in fibers.items():
            occupied = sum(
                bool(sum(sign * (-1) ** ((mask & u).bit_count() % 2)
                         for mask, sign in entries))
                for u in range(1 << len(queries))
            )
            assert weights[label] == Fraction(occupied, 1 << len(queries))
    return weights, quotient


def maximum_weight_basis(weights):
    chosen = []
    for label in sorted(weights, key=lambda point: (-weights[point], point)):
        if rank([*chosen, label]) > len(chosen):
            chosen.append(label)
    return sum(weights[label] for label in chosen)


def verify_prefix(spectrum, queries, expected, next_query=None):
    weights, quotient = fiber_weights(spectrum, queries)
    mean_size = sum(weights.values())
    basis_weight = maximum_weight_basis(weights)
    assert mean_size - basis_weight == expected
    reduction = None
    if next_query is not None:
        collisions = Counter(a ^ b for a, b in combinations(weights, 2))
        representative, _ = reduce_mod(next_query, quotient)
        assert representative
        reduction = collisions[representative]
        assert reduction == max(collisions.values())
    return {
        "query_rank": len(queries),
        "query_vector": next_query,
        "projected_size": len(weights),
        "mean_slice_size": str(mean_size),
        "max_basis_weight": str(basis_weight),
        "tau": str(expected),
        "greedy_reduction": reduction,
    }


def verify_uniform_bases(bands, assignments):
    support = set().union(*map(set, bands))
    target_rank = rank(support)
    appearances = Counter()
    total = 0
    for assignment in assignments:
        for shift in range(len(bands[0])):
            basis = []
            for block, band in enumerate(bands):
                anchor = band[shift]
                basis.append(anchor)
                if block in assignment:
                    partner = anchor ^ assignment[block]
                    assert partner in band
                    basis.append(partner)
            assert rank(basis) == len(basis) == target_rank
            appearances.update(basis)
            total += 1
    probability = Fraction(target_rank, len(support))
    assert all(Fraction(appearances[point], total) == probability for point in support)
    density = 1 / probability
    colors = (density.numerator + density.denominator - 1) // density.denominator
    return {"bases_in_uniform_certificate": total,
            "maximum_density": str(density), "c": colors}


def ordinary(k):
    leaves = 1 << k
    addresses = [1 << (leaves + i) for i in range(k)]
    bands = [[(1 << i) ^ direction for direction in span(addresses)]
             for i in range(leaves)]
    assignments = [
        {(start + i) % leaves: address for i, address in enumerate(addresses)}
        for start in range(leaves)
    ]
    certificate = verify_uniform_bases(bands, assignments)
    assert certificate["c"] == (leaves - k + 1 if k else 1)
    spectrum = {
        frequency: (-1) ** (((frequency >> leaves) & tag).bit_count() % 2)
        for tag, band in enumerate(bands) for frequency in band
    }
    if k <= 3:
        def value(x):
            return (-1) ** ((x >> (x >> leaves)) & 1)
        audit_walsh(spectrum, leaves + k, value)
    rows = []
    for h in range(k + 1):
        width = leaves >> h
        expected = width * width - Fraction(leaves + k - h, 1 << h)
        next_query = addresses[h] if h < k else None
        rows.append(verify_prefix(spectrum, addresses[:h], expected, next_query))
    return {"k": k, **certificate, "prefixes": rows}


def four_line(private_depth):
    block_size = 1 << private_depth
    leaf_count = 4 * block_size
    directions = [1, 6, 4, 3]
    origins = [0, 2, 3, 5]
    address_planes = [
        [linear_form << leaf_count for linear_form in range(8)
         if ((linear_form & direction).bit_count() % 2) == 0]
        for direction in directions
    ]
    private_nodes = {
        (block, heap_index): 1 << (leaf_count + 3
                                   + block * (block_size - 1) + heap_index - 1)
        for block in range(4) for heap_index in range(1, block_size)
    }
    bands = []
    spectrum = {}
    for block in range(4):
        for leaf in range(block_size):
            tag = block * block_size + leaf
            heap_index = block_size + leaf
            paths = []
            while heap_index > 1:
                paths.append((private_nodes[block, heap_index // 2], heap_index % 2))
                heap_index //= 2
            phase = origins[block] << leaf_count
            for direction, outcome in paths:
                if outcome:
                    phase ^= direction
            band = [(1 << tag) ^ address ^ private
                    for address in address_planes[block]
                    for private in span(direction for direction, _ in paths)]
            bands.append(band)
            spectrum.update((point, (-1) ** ((phase & point).bit_count() % 2))
                            for point in band)

    def private_assignment(block, node, omitted, lo, hi):
        if hi - lo == 1:
            return {}
        midpoint = (lo + hi) // 2
        left_skip = omitted if omitted < midpoint else lo
        right_skip = omitted if omitted >= midpoint else midpoint
        assigned = private_assignment(block, 2 * node, left_skip, lo, midpoint)
        assigned.update(private_assignment(block, 2 * node + 1, right_skip, midpoint, hi))
        spare = right_skip if omitted < midpoint else left_skip
        assigned[block * block_size + spare] = private_nodes[block, node]
        return assigned

    assignments = []
    for omitted_tag in range(leaf_count):
        excluded_block, excluded_leaf = divmod(omitted_tag, block_size)
        chosen_blocks = [block for block in range(4) if block != excluded_block]
        address_basis = next(
            values for values in product(*(address_planes[i][1:] for i in chosen_blocks))
            if rank(values) == 3
        )
        allocation = {}
        for block in range(4):
            omitted = excluded_leaf if block == excluded_block else 0
            allocation.update(private_assignment(block, 1, omitted, 0, block_size))
        allocation.update({block * block_size: vector
                           for block, vector in zip(chosen_blocks, address_basis)})
        assert set(allocation) == set(range(leaf_count)) - {omitted_tag}
        assignments.append(allocation)
    certificate = verify_uniform_bases(bands, assignments)
    assert certificate["c"] == leaf_count // 2 + 1

    if private_depth <= 1:
        def value(x):
            address = (x >> leaf_count) & 7
            block = next(i for i in range(4) if address in (origins[i], origins[i] ^ directions[i]))
            node = 1
            while node < block_size:
                node = 2 * node + bool(x & private_nodes[block, node])
            return (-1) ** ((x >> (block * block_size + node - block_size)) & 1)
        audit_walsh(spectrum, 2 * leaf_count - 1, value)

    q1 = next(v for v in range(1 << leaf_count, 8 << leaf_count, 1 << leaf_count)
              if sum(v in plane for plane in address_planes) == 2)
    plane = next(plane for plane in address_planes if q1 in plane)
    q2 = next(v for v in plane if v not in (0, q1))
    q3 = next(v << leaf_count for v in range(1, 8) if v << leaf_count not in plane)
    address_queries = [q1, q2, q3]
    expected = [(leaf_count - 1) ** 2,
                10 * block_size ** 2 - 6 * block_size + 1,
                Fraction(13 * block_size ** 2 - 14 * block_size + 5, 4)]
    rows = [
        verify_prefix(spectrum, address_queries[:h], expected[h], address_queries[h])
        for h in range(3)
    ]
    queries = list(address_queries)
    frontier = [(block, 1, block_size) for block in range(4)]
    while True:
        expected = Fraction(sum(size * (size - 1) ** 2 for _, _, size in frontier),
                            leaf_count)
        block, node, largest = max(frontier, key=lambda item: item[2])
        query = private_nodes[block, node] if largest > 1 else None
        rows.append(verify_prefix(spectrum, queries, expected, query))
        if query is None:
            assert len(queries) == leaf_count - 1
            break
        queries.append(query)
        frontier.remove((block, node, largest))
        frontier.extend([(block, 2 * node, largest // 2),
                         (block, 2 * node + 1, largest // 2)])
    return {"private_depth": private_depth, **certificate, "prefixes": rows}


def pseudospectrum(m):
    common = 1 << m
    first = [1 << i for i in range(m)]
    second = [point ^ common for point in first]
    support = set(first + second)
    assert rank(first) == rank(second) == m
    assert rank(support) == m + 1 < len(support)
    differences = Counter(a ^ b for a, b in combinations(support, 2))
    assert min(differences.values()) >= 2
    for coefficients in product((-1, 1), repeat=2 * m) if m <= 6 else ():
        rows = [
            [a + address_sign * b for a, b in zip(coefficients[:m], coefficients[m:])]
            for address_sign in (-1, 1)
        ]
        assert max(sum(value != 0 for value in row) for row in rows) >= 2
    return {"m": m, "rank": m + 1, "size": 2 * m, "c": 2,
            "maximum_density": str(Fraction(2 * m, m + 1)),
            "minimum_pair_multiplicity": min(differences.values())}


def main():
    result = {
        "python": sys.version.split()[0],
        "ordinary": [ordinary(k) for k in range(6)],
        "four_line": [four_line(k) for k in range(4)],
        "pseudospectrum": [pseudospectrum(m) for m in [3, 4, 5, 6, 16, 128]],
    }
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
