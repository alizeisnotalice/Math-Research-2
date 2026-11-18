"""Directed bounds for optimally chosen support vertices and the odd-Reuleaux cap obstruction."""

import argparse
from fractions import Fraction as Q
from hashlib import sha256
import json
from pathlib import Path
import platform

import flint
from flint import arb, ctx


def rational(q):
    q = Q(q)
    return arb(q.numerator) / q.denominator


def one_row(n):
    m = n + 1 if n % 2 == 0 else n + 2
    pi = arb.pi()
    h = pi / m
    radius = 1 / (2 * (h / 2).cos())
    lower = 2 * (3 * h / 4).sin() * (h / 4).sin()
    assert (lower - radius + (1 - radius) * h.cos()).contains(0)
    assert lower > 0 and radius * (1 - (2 * h).cos()) > lower
    c = 3 * pi / (2 * n)
    # Unique root of sin(a)-2sin(2a-c) on c/2 < a < c.
    lo, hi = Q(1, 2), Q(1)
    for _ in range(110):
        mid = (lo + hi) / 2
        a = c * rational(mid)
        derivative = a.sin() - 2 * (2 * a - c).sin()
        if derivative > 0:
            lo = mid
        else:
            assert derivative < 0
            hi = mid
    a = c * rational(lo).union(rational(hi))
    # Exact cos(2a-c) - cos(a), in a form without loss of leading digits.
    upper = 2 * ((3 * a - c) / 2).sin() * ((c - a) / 2).sin()
    disk = (pi / (2 * n)).sin() ** 2
    simple = 3 * pi**2 / (8 * n * n)
    assert upper < simple and upper > lower
    return dict(
        N=n,
        odd_m=m,
        sharp_uniform_lower=str(lower),
        finite_upper=str(upper),
        quadratic_upper=str(simple),
        disk_error=str(disk),
        scaled_lower=str(n * n * lower),
        scaled_upper=str(n * n * upper),
    )


def main(output):
    assert __debug__
    ctx.prec = 256
    result = dict(
        python=platform.python_version(),
        python_flint=flint.__version__,
        precision_bits=ctx.prec,
        script_sha256=sha256(Path(__file__).read_bytes()).hexdigest(),
        rows=[
            one_row(n) for n in (3, 4, 5, 8, 15, 16, 31, 32, 63, 64, 1000, 18028, 10**6)
        ],
        uniform_leading=str(3 * arb.pi() ** 2 / 8),
        fixed_class_leading=str(arb.pi() ** 2 / 4),
    )
    output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("output", type=Path)
    main(p.parse_args().output)
