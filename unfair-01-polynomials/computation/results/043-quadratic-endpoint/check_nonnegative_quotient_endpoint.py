#!/usr/bin/env python3
"""Finite exact regression for C=R*B^2 in the rational terminal branch.

This is only a small check of the endpoint proof, not a general enumeration.
"""
from __future__ import annotations

import itertools
import json
import platform
import sys
from pathlib import Path


def mul(left: tuple[int, ...], right: tuple[int, ...]) -> tuple[int, ...]:
    out = [0] * (len(left) + len(right) - 1)
    for i, a in enumerate(left):
        for j, b in enumerate(right):
            out[i + j] += a * b
    return tuple(out)


def trim(poly: tuple[int, ...]) -> tuple[int, ...]:
    last = len(poly)
    while last > 1 and poly[last - 1] == 0:
        last -= 1
    return poly[:last]


def main() -> None:
    # Every binary R,B through degree 4 with constant coefficient one.
    examined = accepted = 0
    for rbits in itertools.product((0, 1), repeat=4):
        R = (1,) + rbits
        for bbits in itertools.product((0, 1), repeat=4):
            B = (1,) + bbits
            examined += 1
            C = mul(R, mul(B, B))
            if all(c in (0, 1) for c in C):
                accepted += 1
                assert trim(B) == (1,)
                assert trim(C) == trim(R)
    out = {
        "python": sys.version.split()[0],
        "platform": platform.platform(),
        "range": "binary R,B with constant 1 and displayed length 5 (degrees at most 4)",
        "examined": examined,
        "accepted": accepted,
        "status": "all exact integer convolution checks passed",
    }
    Path(__file__).with_name("check_nonnegative_quotient_endpoint.json").write_text(
        json.dumps(out, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(out, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
