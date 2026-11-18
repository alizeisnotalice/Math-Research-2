"""Recompute the degree-16 saturation obstruction over F_163."""

import argparse
from itertools import product
import json
from math import isqrt
from pathlib import Path
from time import perf_counter


MODULUS = 163
T = -6
H = (1, -4, 6, -7, 4, -2, 1)


def require(condition, detail):
    if not condition:
        raise ValueError(detail)


def evaluate(coefficients, argument):
    result = 0
    for coefficient in coefficients:
        result = result * argument + coefficient
    return result


def constant(value):
    return (value % MODULUS, 0, 0)


def subtract(left, right):
    return tuple((a - b) % MODULUS for a, b in zip(left, right))


def residual(a, b, coefficient):
    index = len(a)
    require(index == len(b), "unequal prefix lengths")
    result = constant(coefficient)
    for i in range(1, index):
        j = index - i
        left, right = a[i], b[j]
        if i < 8:
            require(left[1:] == (0, 0), "unexpected variable in A prefix")
            term = tuple(left[0] * entry for entry in right)
        else:
            require(j < 8 and right[1:] == (0, 0), "nonlinear product")
            term = tuple(right[0] * entry for entry in left)
        result = subtract(result, term)
    return result


def main():
    start = perf_counter()
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    require(all(MODULUS % divisor for divisor in range(2, isqrt(MODULUS) + 1)), "nonprime")
    require(evaluate(H, T) % MODULUS == 0, "not a root of H")
    m = 1 - T + T * T
    n = 1 - T * m
    p = 1 - (1 - T) * n
    v = T * T * m
    q = (1 - n * v) % MODULUS
    require(q == 49, "incorrect degree-8 row")

    rows = []
    for u, bits in product((0, q), product((0, 1), repeat=4)):
        a = list(map(constant, (1, T, 0, n, p, 0, 0, 0, u)))
        b = list(map(constant, (1, 1 - T, m, 0, 0, v, 0, 0, q - u)))
        for k, value in enumerate((1, 1, 1, *bits[:3]), start=9):
            row = residual(a, b, value)
            if k == 9:
                left = (0, 1, 0)
            elif k == 10:
                left = (0, 0, 1)
            elif k in (12, 13):
                left = row
            else:
                left = constant(0)
            a.append(left)
            b.append(subtract(row, left))
        reduced = residual(a, b, bits[3])
        expected = (66 if u == 0 else 142) + sum(a * b for a, b in zip((11, 78, 6, 1), bits))
        require(reduced == constant(expected), (u, bits, reduced, expected))
        require(reduced[0] != 0, (u, bits, reduced))
        rows.append({"A8": u, "C12_through_C15": bits, "R15": reduced})
    require(len(rows) == 32, "incomplete coverage")
    report = {"prime": MODULUS, "t_image": T, "H_at_t_integer": evaluate(H, T),
              "q_image": q, "count": len(rows), "rows": rows, "status": "verified"}
    with args.out.open("x") as output:
        json.dump(report, output, indent=2)
        output.write("\n")
    print(json.dumps({"status": "verified", "count": len(rows), "seconds": perf_counter() - start}))


if __name__ == "__main__":
    main()
