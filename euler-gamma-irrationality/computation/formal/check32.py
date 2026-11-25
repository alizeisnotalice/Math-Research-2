"""Independent check for the 33-by-2 integer basis; uses the original first-round columns."""

import json
import sys
from pathlib import Path
from math import gcd
from itertools import combinations
import sympy as s

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from sondow import case


def verify(cert):
    n = cert["N"]
    rows = [case(i) for i in range(n, 2 * n + 1)]
    primes = sorted(set().union(*(r["prime_exponents"] for r in rows)))
    E = s.Matrix([[r["prime_exponents"].get(p, 0) for r in rows] for p in primes])
    assert len(E.nullspace()) == cert["nullity"]
    if n != 32:
        assert cert["nullity"] == 0
        return
    x, y = [r["v"] for r in cert["vectors"]]
    assert E * s.Matrix(x) == E * s.Matrix(y) == s.zeros(len(primes), 1)
    minors = (x[i] * y[j] - x[j] * y[i] for i, j in combinations(range(len(x)), 2))
    assert gcd(*minors) == 1  # index in saturation is the last determinantal divisor
    for vector, rec in zip((x, y), cert["vectors"]):
        for key, dest in [("D", "B"), ("dA", "A")]:
            assert sum(c * r[key] for c, r in zip(vector, rows)) == rec[dest]
    a, b = cert["vectors"]
    assert a["A"] * b["B"] - a["B"] * b["A"] == cert["image_det"] != 0


if __name__ == "__main__":
    for cert in json.loads(Path(__file__).with_name("kernel.json").read_text()):
        verify(cert)
    print("verified kernels, saturation and image determinant")
