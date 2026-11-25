"""Outward dyadic enclosures for integer forms in 1, gamma and log(k)."""

import argparse
import json
from fractions import Fraction as F
from functools import cache
from math import ceil, comb, floor


def outward(lo, hi, bits):
    scale = 1 << bits
    return F(floor(lo * scale), scale), F(ceil(hi * scale), scale)


@cache
def harmonic(n):
    return sum((F(1, k) for k in range(1, n + 1)), F(0))


def log_series(x, bits):
    t = (x - 1) / (x + 1)
    total, power, j = F(0), t, 1
    while True:
        total += 2 * power / j
        power *= t * t
        j += 2
        tail = 2 * power / (j * (1 - t * t))
        if tail <= F(1, 1 << bits):
            return outward(total, total + tail, bits)


@cache
def log_interval(k, bits):
    if k < 1:
        raise ValueError("positive log argument required")
    e = k.bit_length() - 1
    a, b = log_series(F(2), bits)
    c, d = log_series(F(k, 1 << e), bits)
    return e * a + c, e * b + d


@cache
def bernoulli(n):
    if n == 0:
        return F(1)
    return -sum(comb(n + 1, k) * bernoulli(k) for k in range(n)) / (n + 1)


@cache
def gamma_interval(bits):
    n = 1 << bits.bit_length()
    a, b = log_interval(n, bits + 8)
    center = harmonic(n) - F(1, 2 * n)
    for k in range(1, n + 1):
        term = bernoulli(2 * k) / (2 * k * n ** (2 * k))
        if abs(term) < F(1, 1 << bits):
            return outward(center - b + min(term, 0), center - a + max(term, 0), bits)
        center += term
    raise ArithmeticError("remainder target unmet")


def certify(g, c, logs, mode="sign", start=32, limit=512):
    if start < 4 or limit < start or mode not in ("sign", "floor"):
        raise ValueError("invalid precision range or decision")
    bits = start
    while True:
        lo = hi = F(c)
        terms = [(v, log_interval(k, bits)) for k, v in logs.items() if v]
        if g:
            terms.append((g, gamma_interval(bits)))
        for coeff, (a, b) in terms:
            left, right = coeff * a, coeff * b
            lo += min(left, right)
            hi += max(left, right)
        if mode == "floor":
            decision = floor(lo) if floor(lo) == floor(hi) else None
        else:
            decision = "positive" if lo > 0 else "negative" if hi < 0 else None
        if decision is not None or bits == limit:
            return {"bits": bits, "lo": str(lo), "hi": str(hi), "decision": decision}
        bits = min(limit, 2 * bits)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("gamma", type=int)
    parser.add_argument("constant", type=int)
    parser.add_argument(
        "--log", nargs=2, type=int, action="append", default=[], metavar=("K", "COEFF")
    )
    parser.add_argument("--mode", choices=["sign", "floor"], default="sign")
    parser.add_argument("--start", type=int, default=32)
    parser.add_argument("--limit", type=int, default=512)
    args = parser.parse_args()
    logs = {}
    for k, v in args.log:
        logs[k] = logs.get(k, 0) + v
    print(
        json.dumps(
            certify(args.gamma, args.constant, logs, args.mode, args.start, args.limit)
        )
    )
