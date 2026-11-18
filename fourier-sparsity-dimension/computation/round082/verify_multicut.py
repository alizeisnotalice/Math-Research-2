"""Exact rank-loss and low-overlap tests for simultaneous Boolean queries."""

import json
import sys
from collections import Counter
from functools import cache


@cache
def rank(vectors):
    pivots = {}
    for value in vectors:
        while value:
            bit = value.bit_length()
            if bit not in pivots:
                pivots[bit] = value
                break
            value ^= pivots[bit]
    return len(pivots)


def rank_of(vectors):
    return rank(tuple(sorted(set(vectors))))


@cache
def character_masks(size):
    return [
        sum(((point & frequency).bit_count() % 2) << point for point in range(size))
        for frequency in range(size)
    ]


@cache
def support(table, size):
    return tuple(
        frequency for frequency, mask in enumerate(character_masks(size))
        if (table ^ mask).bit_count() * 2 != size
    )


def lifted(points, dim):
    return tuple((1 << dim) | point for point in points)


def fibre_kernel(fibres, active, dim, query_dim):
    """Independent calculation from active-fibre differences and lift relations."""
    differences = []
    quotient_pivots = {}
    relation_lifts = []
    projected_height = 1 << (dim - query_dim)
    for label in active:
        fibre = fibres[label]
        representative = fibre[0]
        differences.extend(point ^ representative for point in fibre[1:])
        image = projected_height | label
        lift = (1 << dim) | representative
        while image:
            bit = image.bit_length()
            if bit not in quotient_pivots:
                quotient_pivots[bit] = (image, lift)
                break
            pivot_image, pivot_lift = quotient_pivots[bit]
            image ^= pivot_image
            lift ^= pivot_lift
        if not image:
            assert lift < 1 << query_dim
            relation_lifts.append(lift)
    difference_rank = rank_of(differences)
    full_kernel_rank = rank_of(differences + relation_lifts)
    return difference_rank, full_kernel_rank - difference_rank, full_kernel_rank


def weighted_basis_cost(labels, projected_dim, weights, default=0):
    selected = []
    for beta in sorted(labels, key=lambda label: (-weights[label], label)):
        vector = (1 << projected_dim) | beta
        if rank_of(selected + [vector]) > len(selected):
            selected.append(vector)
    return sum(weights[vector ^ (1 << projected_dim)] - default for vector in selected)


@cache
def analyse(spectrum, children, dim, query_dim):
    m = 1 << query_dim
    data_dim = dim - query_dim
    lifted_parent = lifted(spectrum, dim)
    parent_rank = rank_of(lifted_parent)
    fibres = {}
    for frequency in spectrum:
        fibres.setdefault(frequency >> query_dim, []).append(frequency)
    all_labels = set(fibres)
    assert all_labels == set().union(*map(set, children))
    projected_rank = rank_of(lifted(all_labels, data_dim))
    h = parent_rank - projected_rank
    assert 0 <= h <= query_dim

    child_ranks = []
    split = Counter()
    for row in children:
        parent_retained = lifted(tuple(e for b in row for e in fibres[b]), dim)
        row_rank = rank_of(lifted(row, data_dim))
        d_rank, relations, kernel_rank = fibre_kernel(fibres, row, dim, query_dim)
        assert kernel_rank == rank_of(parent_retained) - row_rank
        deleted_rank = parent_rank - rank_of(parent_retained)
        assert parent_rank - row_rank == deleted_rank + d_rank + relations
        assert parent_rank - row_rank == h + (projected_rank - row_rank)
        child_ranks.append(row_rank)
        split["rows_with_deletion_rank"] += deleted_rank > 0
        split["rows_with_fibre_direction_rank"] += d_rank > 0
        split["rows_with_relation_only_residual"] += relations > 0

    active_counts = Counter(label for row in children for label in row)
    assert all(n * active_counts[b] >= m for b,n in
               ((label, len(fibre)) for label,fibre in fibres.items()))
    overlap = {label for label,count in active_counts.items() if count >= 2}
    overlap_rank = rank_of(lifted(overlap, data_dim))
    shared_weight = weighted_basis_cost(overlap, data_dim, active_counts, default=1)
    weighted = weighted_basis_cost(all_labels, data_dim, active_counts)
    assert sum(child_ranks) >= weighted
    repeat_rank = sum(child_ranks) - projected_rank
    core_vectors = lifted(overlap, data_dim)
    common_intersections = sum(
        row_rank + overlap_rank - rank_of(core_vectors + lifted(row, data_dim))
        for row_rank,row in zip(child_ranks,children)
    )
    assert repeat_rank >= common_intersections - overlap_rank >= shared_weight >= overlap_rank

    row_sets = [set(row) for row in children]
    hierarchical_rank = 0
    repeated_counts = Counter()
    while len(row_sets) > 1:
        merged_rows = []
        for offset in range(0, len(row_sets), 2):
            left, right = row_sets[offset:offset+2]
            union = left | right
            hierarchical_rank += (rank_of(lifted(left, data_dim))
                                  + rank_of(lifted(right, data_dim))
                                  - rank_of(lifted(union, data_dim)))
            repeated_counts.update(left & right)
            merged_rows.append(union)
        row_sets = merged_rows
    assert hierarchical_rank == repeat_rank
    assert repeated_counts == Counter({b: t-1 for b,t in active_counts.items() if t > 1})

    parent_size = len(spectrum)
    child_size_total = sum(map(len, children))
    U = sum(child_ranks) + m
    delta = U - parent_rank - 1
    deficit = m * child_size_total - parent_size
    assert delta >= m - h - 1 + shared_weight
    assert deficit == sum(m * active_counts[b] - len(fibre) for b,fibre in fibres.items())
    assert deficit <= (m*m-1)*len(overlap)
    assert 2*deficit <= 3*m*sum(count-1 for count in active_counts.values())
    true_improvement = any(
        parent_size * (r+1)**2 >= len(row)*(parent_rank+1)**2
        for row,r in zip(children,child_ranks)
    )
    low_overlap = (
        m * child_size_total * (m-h-1+shared_weight) >=
        (m*m-1)*len(overlap)*U
    )
    common_identity = (
        m*child_size_total*(2*U*delta-delta*delta) >= deficit*U*U
    )
    squared_sum = (
        parent_size * sum((r+1)**2 for r in child_ranks) >=
        child_size_total * (parent_rank+1)**2
    )
    assert not low_overlap or common_identity
    assert not common_identity or squared_sum
    assert not squared_sum or true_improvement
    return split, {
        "verified_low_overlap_sufficient": low_overlap,
        "verified_aggregate_sufficient": common_identity,
        "verified_squared_sum_sufficient": squared_sum,
        "actual_non_decreasing_child": true_improvement,
    }


def cut_table(table, dim, query_dim):
    m = 1 << query_dim
    data_size = 1 << (dim-query_dim)
    return tuple(
        support(sum(((table >> (address + m*x)) & 1) << x
                    for x in range(data_size)), data_size)
        for address in range(m)
    )


def main():
    results = {"python": sys.version.split()[0], "four_variable_layouts": []}
    for query_dim in (1, 2, 3):
        counts = Counter()
        for table in range(1 << 16):
            rows = cut_table(table, 4, query_dim)
            split, decisions = analyse(support(table, 16), rows, 4, query_dim)
            counts.update(split)
            counts.update({key:int(value) for key,value in decisions.items()})
        results["four_variable_layouts"].append({
            "query_dim": query_dim,
            "tables": 1 << 16,
            "rows": (1 << 16)*(1 << query_dim),
            **counts,
        })

    table = 0
    for x in range(16):
        x1,x2,x3,x4 = [(x >> j) & 1 for j in range(4)]
        table |= ((x1*x2) ^ (x3*x4)) << (2*x)
        table |= ((x1*x3) ^ (x2*x4)) << (2*x+1)
    points = support(table, 32)
    rows = cut_table(table, 5, 1)
    split, decisions = analyse(points,rows,5,1)
    assert len(points) == 16 and len({e >> 1 for e in points}) == 16
    assert rank_of(lifted(points,5)) == 6
    assert [rank_of(lifted(row,4)) for row in rows] == [5,5]
    assert split["rows_with_relation_only_residual"] == 2
    assert not decisions["actual_non_decreasing_child"]
    results["bent_switch"] = {
        "parent_support": points,
        "parent_rank": 6,
        "child_ranks": [5,5],
        "both_child_support_sizes": 16,
        "singleton_fibres": 16,
    }
    print(json.dumps(results, indent=2))


if __name__ == "__main__":
    main()
