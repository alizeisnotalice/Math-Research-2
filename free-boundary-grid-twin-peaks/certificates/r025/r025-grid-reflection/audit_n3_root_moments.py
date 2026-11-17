#!/usr/bin/env python3
"""Independent exact n=3 audit for the root-pair moment normalizations.

Enumerates all 9! vertex labelings (the tuple stores labels in row-major
vertex order) and prints only integer numerators.  It is deliberately small:
the calculation is a cross-check of the ordered-root factors in the analytic
first-merge identities, not asymptotic evidence.
"""
from fractions import Fraction
from itertools import permutations
import json


N = 3
VERTICES = [(x, y) for y in range(N) for x in range(N)]
INDEX = {v: i for i, v in enumerate(VERTICES)}
NEIGHBORS = []
for x, y in VERTICES:
    NEIGHBORS.append(
        [INDEX[z] for z in ((x - 1, y), (x + 1, y), (x, y - 1), (x, y + 1))
         if z in INDEX]
    )


def main():
    pairs = {}
    a2 = 0
    product_sum = 0  # Sum_{p<q} K(p,q)t_x(p)t_x(q).
    square_sum = 0   # Sum_{p<q} K(p,q)(t_x(p)^2+t_x(q)^2).
    for labels in permutations(range(N * N)):
        peaks = [i for i, ns in enumerate(NEIGHBORS)
                 if all(labels[i] > labels[j] for j in ns)]
        if len(peaks) != 2:
            continue
        a2 += 1
        p, q = sorted(peaks)
        pairs[(p, q)] = pairs.get((p, q), 0) + 1
        tp = 2 * VERTICES[p][0] - (N - 1)
        tq = 2 * VERTICES[q][0] - (N - 1)
        product_sum += tp * tq
        square_sum += tp * tp + tq * tq

    # The ordered-root law has E[t(X)t(Y)]=product_sum/A2 and
    # E[t(X)^2]=square_sum/(2*A2).
    p00, p22, p02 = INDEX[(0, 0)], INDEX[(2, 2)], INDEX[(0, 2)]
    out = {
        "n": N,
        "A2": a2,
        "K_00_22": pairs[tuple(sorted((p00, p22)))],
        "K_00_02": pairs[tuple(sorted((p00, p02)))],
        "tx_product_numerator": product_sum,
        "E_tx_ty": str(Fraction(product_sum, a2)),
        "tx_square_numerator": square_sum,
        "E_tx_squared": str(Fraction(square_sum, 2 * a2)),
        "checks": {
            "A2_expected": a2 == 143808,
            "reflection_witness_expected": (
                pairs[tuple(sorted((p00, p22)))] == 7872
                and pairs[tuple(sorted((p00, p02)))] == 2688
            ),
            "product_expected": Fraction(product_sum, a2) == Fraction(-170, 107),
            "square_expected": Fraction(square_sum, 2 * a2) == Fraction(2110, 749),
        },
    }
    print(json.dumps(out, indent=2, sort_keys=True))
    assert all(out["checks"].values())


if __name__ == "__main__":
    main()
