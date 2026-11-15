"""Exact first two peak factorial moments on free square grids.

The calculation uses the iid-U[0,1] representation of a uniform labelling:
a vertex is a peak when its value exceeds every neighbour.
"""
from fractions import Fraction
import json
from math import factorial


def grid(n):
    adj = []
    for y in range(n):
        for x in range(n):
            adj.append({Y*n + X for X, Y in ((x-1, y), (x+1, y),
                       (x, y-1), (x, y+1)) if 0 <= X < n and 0 <= Y < n})
    return adj


def two_peak_probability(adj, u, v):
    if v in adj[u]:
        return Fraction()
    c = len(adj[u] & adj[v])
    a = len(adj[u]) - c
    b = len(adj[v]) - c
    # Integral over values s,t at u,v, split at s=t:
    # int s^a t^b min(s,t)^c ds dt.
    return (Fraction(1, a+c+1) + Fraction(1, b+c+1)) / (a+b+c+2)


def main():
    rows = []
    for n in range(2, 16):
        adj = grid(n)
        mean = sum((Fraction(1, len(a)+1) for a in adj), Fraction())
        pair = sum((two_peak_probability(adj, u, v)
                    for u in range(n*n) for v in range(u+1, n*n)), Fraction())
        if n >= 4:
            closed = (Fraction(n**4, 50) + Fraction(n**3, 25)
                      - Fraction(11*n*n, 450) - Fraction(n, 15)
                      - Fraction(13, 315))
            assert pair == closed
        rows.append({"n": n, "vertices": n*n,
                     "E_K": str(mean), "E_choose_K_2": str(pair),
                     "E_choose_K_2_over_N2": str(pair/(n*n)**2)})
    print(json.dumps(rows, indent=2))


if __name__ == '__main__':
    main()
