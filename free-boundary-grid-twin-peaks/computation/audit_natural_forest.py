#!/usr/bin/env python3
"""Exact small-grid audit for the higher-neighbour natural forest.

For a fixed labelled graph, every nonpeak independently selects one of its
higher-labelled neighbours uniformly.  The resulting directed forest flows to
the local maxima.  This script records an exact n=2 basin decomposition and
the n=3 conditional entropy integrand without using any previous peak DP.
"""

from collections import Counter, defaultdict
from fractions import Fraction
from itertools import permutations, product
import json
from math import factorial
from pathlib import Path
import sys


def grid(n):
    adj = []
    for y in range(n):
        for x in range(n):
            q = []
            for dx, dy in ((1, 0), (-1, 0), (0, 1), (0, -1)):
                xx, yy = x + dx, y + dy
                if 0 <= xx < n and 0 <= yy < n:
                    q.append(yy * n + xx)
            adj.append(tuple(q))
    return adj


def roots_and_choices(lab, adj):
    roots = tuple(v for v, ns in enumerate(adj)
                  if all(lab[v] > lab[u] for u in ns))
    choices = [tuple(u for u in ns if lab[u] > lab[v]) for v, ns in enumerate(adj)]
    return roots, choices


def basin_key(roots, parent, adj):
    """Return all basin masks, with internal coverage/connectedness checks."""
    vertices = set(roots) | set(parent)
    full = sum(1 << v for v in vertices)
    masks = {r: 1 << r for r in roots}
    for v in parent:
        u = v
        seen = set()
        while u not in masks:
            assert u not in seen
            seen.add(u)
            u = parent[u]
        masks[u] |= 1 << v
    union = 0
    for r, mask in masks.items():
        assert mask & (1 << r)
        assert not (union & mask)
        union |= mask
        reached = {r}
        while True:
            grown = reached | {w for v in reached for w in adj[v] if mask & (1 << w)}
            if grown == reached:
                break
            reached = grown
        assert reached == {v for v in vertices if mask & (1 << v)}
    assert union == full
    return tuple(sorted(masks.items()))


def n2_forest_decomposition():
    adj = grid(2)
    mass = defaultdict(Fraction)
    forest_mass = defaultdict(Fraction)
    per_label = []
    conditional_partitions = []
    k2 = 0
    for lab in permutations(range(4)):
        roots, choices = roots_and_choices(lab, adj)
        if len(roots) != 2:
            continue
        k2 += 1
        nonroots = [v for v in range(4) if v not in roots]
        total = Fraction()
        conditional = defaultdict(Fraction)
        for selected in product(*(choices[v] for v in nonroots)):
            parent = {v: u for v, u in zip(nonroots, selected)}
            denom = 1
            for v in nonroots:
                denom *= len(choices[v])
            p = Fraction(1, denom)
            total += p
            key = basin_key(roots, parent, adj)
            mass[key] += p
            forest_mass[(roots, tuple(sorted(parent.items())))] += p
            conditional[key] += p
        assert total == 1
        assert len(conditional) == 4
        assert set(conditional.values()) == {Fraction(1, 4)}
        conditional_partitions.append(conditional)
        per_label.append({"labels": list(lab), "roots": list(roots), "choice_mass": str(total)})
    assert k2 == 8 and sum(mass.values(), Fraction()) == 8
    assert len(mass) == 8 and set(mass.values()) == {Fraction(1)}
    # A second implementation uses random local priority lists.  For a fixed
    # directed forest F, a priority profile contributes all labelings obeying
    # the associated strict inequalities; its count is a linear-extension
    # count.  This is the finite n=2 audit of the general poset formula in
    # natural-region-weights.md.
    def extensions(inequalities):
        return sum(
            all(lab[hi] > lab[lo] for hi, lo in inequalities)
            for lab in permutations(range(4))
        )

    poset_mass = {}
    for (roots, items), direct in forest_mass.items():
        parent = dict(items)
        nonroots = [v for v in range(4) if v not in roots]
        denominator = 1
        order_choices = []
        for v in nonroots:
            denominator *= factorial(len(adj[v]))
            order_choices.append(list(permutations(adj[v])))
        numerator = 0
        for profile in product(*order_choices):
            inequalities = []
            for r in roots:
                inequalities.extend((r, u) for u in adj[r])
            for v, order in zip(nonroots, profile):
                p = parent[v]
                inequalities.append((p, v))
                for u in order[:order.index(p)]:
                    inequalities.append((v, u))
            numerator += extensions(inequalities)
        poset_mass[(roots, items)] = Fraction(numerator, denominator)
        assert poset_mass[(roots, items)] == direct

    return {
        "K2_labellings": k2,
        "weighted_forest_mass": str(sum(mass.values(), Fraction())),
        "basins": [
            {"root_masks": [[r, m] for r, m in key], "weight": str(value)}
            for key, value in sorted(mass.items())
        ],
        "per_labelling_normalization": per_label,
        "preference_poset_formula": {
            "all_directed_forests_agree": True,
            "forest_count": len(forest_mass),
            "total": str(sum(poset_mass.values(), Fraction())),
        },
        "entropy_direction_audit": {
            "rooted_partitions": len(mass),
            "marginal_distribution": "uniform on 8 rooted partitions",
            "conditional_distribution_given_each_K2_labelling": "uniform on 4 rooted partitions",
            "H_Pi": "log(8)",
            "H_Pi_given_L": "log(4)",
            "I_Pi_L": "log(2)",
            "identity": "H(Pi)=H(Pi|L)+I(Pi;L)",
        },
    }


def n3_entropy_integrand():
    adj = grid(3)
    hist = Counter()
    profile = Counter()
    distances = defaultdict(lambda: [0, 0, 0])  # count, sum e2, sum e3
    k2 = 0
    for lab in permutations(range(9)):
        roots, choices = roots_and_choices(lab, adj)
        if len(roots) != 2:
            continue
        k2 += 1
        e2 = e3 = c3 = c4 = 0
        for v in range(9):
            if v in roots:
                continue
            d = len(choices[v])
            assert 1 <= d <= 4
            c3 += (d == 3)
            c4 += (d == 4)
            while d % 2 == 0:
                e2 += 1
                d //= 2
            while d % 3 == 0:
                e3 += 1
                d //= 3
            assert d == 1
        hist[e2, e3] += 1
        profile[c3, c4] += 1
        # For P_3, |E|-(N-2)=12-7=5.  This checks the exact edge budget.
        assert 2 ** e2 * 3 ** e3 * 4 ** c3 * 2 ** c4 == 2 ** 5 * 3 ** c3
        a, b = roots
        rho = abs(a % 3 - b % 3) + abs(a // 3 - b // 3)
        z = distances[rho]
        z[0] += 1
        z[1] += e2
        z[2] += e3
    assert k2 == 143808
    sum2 = sum(e2 * c for (e2, _), c in hist.items())
    sum3 = sum(e3 * c for (_, e3), c in hist.items())
    return {
        "K2_labellings": k2,
        "histogram_product_higher_neighbour_counts": [
            {"power_2": a, "power_3": b, "count": c}
            for (a, b), c in sorted(hist.items())
        ],
        "mean_entropy_coefficients": {
            "log_2": [sum2, k2],
            "log_3": [sum3, k2],
        },
        "edge_budget": {
            "sum_higher_degrees": 12,
            "nonpeaks": 7,
            "exponent_bound": 5,
            "maximum_product": max(2 ** a * 3 ** b for a, b in hist),
            "maximum_count": sum(c for (a, b), c in hist.items() if 2 ** a * 3 ** b == 32),
            "c3_c4_histogram": [
                {"c3": a, "c4": b, "count": c}
                for (a, b), c in sorted(profile.items())
            ],
        },
        "by_manhattan_peak_distance": {
            str(r): {"count": c, "sum_log_2": x, "sum_log_3": y}
            for r, (c, x, y) in sorted(distances.items())
        },
    }


def main():
    ans = {
        "definition": "For each nonpeak v select a higher neighbour uniformly; roots are local maxima.",
        "n2_exact_forest_partition": n2_forest_decomposition(),
        "n3_conditional_entropy_integrand": n3_entropy_integrand(),
    }
    out = Path(sys.argv[1]) if len(sys.argv) == 2 else Path("natural-forest-audit.json")
    out.write_text(json.dumps(ans, indent=2) + "\n")
    print(json.dumps({"status": "passed", "output": str(out), "n2_mass": 8, "n3_K2": 143808}))


if __name__ == "__main__":
    main()
