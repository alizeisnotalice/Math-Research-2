"""Exact integer-polynomial audit of the L=3 tensor carry obstruction."""
from collections import defaultdict


def add(*polys):
    out = defaultdict(int)
    for p in polys:
        for k, v in p.items():
            out[k] += v
    return {k: v for k, v in out.items() if v}


def mul(p, q):
    out = defaultdict(int)
    for (i, j), a in p.items():
        for (k, ell), b in q.items():
            out[i + k, j + ell] += a * b
    return {k: v for k, v in out.items() if v}


def main():
    A = {(0, 0): 1, (1, 0): 1}
    B = {(0, 0): 1, (2, 0): 1, (1, 1): 1}
    P_Q = {(0, 0): 1, (1, 0): 1, (2, 0): 1,
           (0, 1): 1, (1, 1): 1, (2, 1): 1}
    product = mul(A, B)
    diff = add(product, {(3, 0): -1}, {(0, 1): 1})
    assert product == {(0, 0): 1, (1, 0): 1, (2, 0): 1,
                       (3, 0): 1, (1, 1): 1, (2, 1): 1}
    assert diff == P_Q
    # The difference product - P_Q is u^3-v and vanishes after v=x^3.
    assert add(product, {k: -v for k, v in P_Q.items()}) == {(3, 0): 1, (0, 1): -1}
    print({"status": "passed", "L": 3,
           "difference": "u^3-v", "canonical_factor_lifts_nonnegative": True})


if __name__ == "__main__":
    main()
