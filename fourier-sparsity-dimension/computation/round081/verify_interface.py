"""Exact audit of Boolean slices and random-basis budgets across a common interface."""

import json
import sys
from collections import Counter
from fractions import Fraction
from functools import cache


@cache
def span(vectors):
    generated = {0}
    for vector in vectors:
        if vector not in generated:
            generated |= {previous ^ vector for previous in generated}
    return frozenset(generated)


def rank(vectors):
    return len(span(tuple(sorted(vectors)))).bit_length() - 1


@cache
def characters(size):
    return [
        sum(((argument & frequency).bit_count() % 2) << argument
            for argument in range(size))
        for frequency in range(size)
    ]


@cache
def support(table, size):
    return tuple(
        frequency
        for frequency, mask in enumerate(characters(size))
        if (table ^ mask).bit_count() * 2 != size
    )


@cache
def basis_cost(points):
    """Squared minimum norm via maximal densest successive rank layers."""
    unassigned = tuple(point for point in points if point)
    assigned = ()
    cost = Fraction(0)
    while unassigned:
        old_rank = rank(assigned)
        best_key = (Fraction(0), 0)
        best_layer = ()
        best_rank = old_rank
        for mask in range(1, 1 << len(unassigned)):
            layer = tuple(point for index, point in enumerate(unassigned)
                          if mask >> index & 1)
            new_rank = rank(assigned + layer)
            increment = new_rank - old_rank
            assert increment > 0
            key = (Fraction(len(layer), increment), len(layer))
            if key > best_key:
                best_key = key
                best_layer = layer
                best_rank = new_rank
        cost += Fraction((best_rank - old_rank)**2, len(best_layer))
        assigned += best_layer
        selected = set(best_layer)
        unassigned = tuple(point for point in unassigned if point not in selected)
    return cost


def residual_cost(core, side_interfaces, activities):
    options = []
    for order in ((0, 1), (1, 0)):
        current = tuple(core)
        total = Fraction(0)
        for side in order:
            larger = current + tuple(side_interfaces[side])
            total += (1 + activities[side]) * (rank(larger) - rank(current))
            current = larger
        options.append(total)
    return min(options)


def audit_layout(core_dim, first_dim, second_dim):
    dimension = core_dim + first_dim + second_dim
    site_count = 1 << dimension
    address_count = 1 << core_dim
    data_size = 1 << (first_dim + second_dim)
    first_mask = (1 << first_dim) - 1
    stats = Counter()
    for table in range(1 << site_count):
        points = support(table, site_count)
        quotient = tuple(point >> core_dim for point in points)
        if any((value & first_mask) and (value >> first_dim) for value in quotient):
            continue

        sides = [
            tuple(point for point, value in zip(points, quotient)
                  if value & first_mask),
            tuple(point for point, value in zip(points, quotient)
                  if value >> first_dim),
        ]
        core_points = tuple(point for point in points if point >> core_dim == 0)
        outer_ranks = [rank(tuple(point >> core_dim for point in side))
                       for side in sides]
        interface = span(points) & frozenset(range(address_count))
        side_interfaces = [span(side) & interface for side in sides]
        assert span(core_points + tuple(side_interfaces[0]) + tuple(side_interfaces[1])) == interface

        costs = [Fraction(0), Fraction(0)]
        activities = [Fraction(0), Fraction(0)]
        for address in range(address_count):
            row_table = sum(
                ((table >> (address + address_count * datum)) & 1) << datum
                for datum in range(data_size)
            )
            row_spectrum = support(row_table, data_size)
            first = any(value & first_mask for value in row_spectrum)
            second = any(value >> first_dim for value in row_spectrum)
            assert not (first and second)
            if first or second:
                side = int(second)
                costs[side] += basis_cost(row_spectrum) / address_count
                activities[side] += Fraction(1, address_count)
            else:
                assert row_spectrum == (0,)
        assert sum(activities) <= 1
        assert all(r * r <= len(side) * q
                   for r, side, q in zip(outer_ranks, sides, costs))

        bound = sum(costs) + basis_cost(core_points)
        bound += residual_cost(core_points, side_interfaces, activities)
        assert basis_cost(points) <= bound
        stats["supported_functions"] += 1
        stats["true_restrictions"] += address_count
        if all(sides):
            stats["both_global_sides"] += 1
        if rank(core_points) < rank(interface):
            stats["interface_not_spanned_by_actual_core_points"] += 1
        if basis_cost(points) == bound:
            stats["full_budget_equalities"] += 1
    return dict(stats)


def main():
    verify = {
        "python": sys.version.split()[0],
        "parameters": ["core_dim", "first_dim", "second_dim"],
        "layouts": [],
    }
    for layout in ((0, 1, 3), (0, 2, 2), (1, 1, 2), (2, 1, 1)):
        verify["layouts"].append({"layout": layout, **audit_layout(*layout)})
    print(json.dumps(verify, indent=2))


if __name__ == "__main__":
    main()
