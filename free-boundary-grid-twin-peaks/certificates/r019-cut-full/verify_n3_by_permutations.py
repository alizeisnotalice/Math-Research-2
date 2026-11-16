#!/usr/bin/env python3
"""Independent 3x3 full-permutation audit for cut covariance numerators."""
from itertools import permutations
import json
from pathlib import Path


def main():
    n = 3
    N = n * n
    cuts = n - 1
    adj = [[] for _ in range(N)]
    for y in range(n):
        for x in range(n):
            v = y * n + x
            for X, Y in ((x - 1, y), (x + 1, y), (x, y - 1), (x, y + 1)):
                if 0 <= X < n and 0 <= Y < n:
                    adj[v].append(Y * n + X)
    A2 = 0
    L = [0] * cuts
    H = [[0] * cuts for _ in range(cuts)]
    for labels in permutations(range(N)):
        peaks = [v for v in range(N) if all(labels[v] > labels[u] for u in adj[v])]
        if len(peaks) != 2:
            continue
        A2 += 1
        p, q = peaks
        for j in range(cuts):
            fp = -1 if p % n <= j else 1
            fq = -1 if q % n <= j else 1
            L[j] += fp + fq
            for k in range(cuts):
                gp = -1 if p % n <= k else 1
                gq = -1 if q % n <= k else 1
                H[j][k] += fp * gq + fq * gp
    den = 4 * A2 * A2
    cov = [[2 * A2 * H[j][k] - L[j] * L[k] for k in range(cuts)] for j in range(cuts)]
    out = {"n": n, "method": "all 9! label permutations", "A2": str(A2),
           "cut_marginal_numerators": list(map(str, L)),
           "ordered_pair_cross_numerators": [[str(x) for x in row] for row in H],
           "covariance_common_denominator": str(den),
           "covariance_numerators": [[str(x) for x in row] for row in cov]}
    target = Path(__file__).with_name("n3-permutation-audit.json")
    target.write_text(json.dumps(out, indent=2) + "\n")
    print(json.dumps(out, indent=2))


if __name__ == "__main__":
    main()
