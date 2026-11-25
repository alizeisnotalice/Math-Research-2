"""Exact small-block experiment for cancellation of all prime logarithms."""

from functools import reduce
from math import gcd

import sympy as sp
from intervals import certify
from sondow import case


def inspect(N):
    rows = [case(n) for n in range(N, 2 * N + 1)]
    primes = sorted(set().union(*(r["prime_exponents"] for r in rows)))
    M = sp.Matrix([[r["prime_exponents"].get(p, 0) for r in rows] for p in primes])
    lattice = []
    for x in M.nullspace():
        v = [int(z * sp.ilcm(1, *[a.q for a in x])) for z in x]
        g = reduce(gcd, v)
        lattice.append([z // g for z in v])
    if not lattice:
        return {"N": N, "rank": M.cols - len(lattice), "nullity": 0}
    # This is a sublattice; basis vectors give candidates, without a shortest-vector claim.
    v = min(lattice, key=lambda row: sum(abs(z) for z in row))
    if M * sp.Matrix(v) != sp.zeros(len(primes), 1):
        raise ArithmeticError("failed cancellation")
    b = sum(z * r["D"] for z, r in zip(v, rows))
    a = sum(z * r["dA"] for z, r in zip(v, rows))
    return {
        "N": N,
        "rank": M.cols - len(lattice),
        "nullity": len(lattice),
        "coeff": v,
        "b": b,
        "a": a,
        "form": certify(b, -a, {}, start=128),
        "sum_abs": sum(map(abs, v)),
    }


if __name__ == "__main__":
    import json

    print(json.dumps([inspect(n) for n in [3, 8, 32]], indent=2))
