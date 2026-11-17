#!/usr/bin/env python3
"""Independent exact AB prefix tree for the S234 degree-seven core.

The fixed core is P=[1,L,a,L*a+k,d,e,f,1].  A=G P and B=Q Pbar are
required coefficientwise in [0,1], and C=A B is binary.  Coefficients of
G,Q are arbitrary integers, so a nonrational row has exactly one possible
integer translate, while a rational row retains both 0/1 endpoints.

This is a necessary-prefix engine only.  A nonempty tree at a finite depth is
recorded as unknown, not as a counterexample.  It deliberately does not
import the parent's HGQ or annulus programs.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import sys
import time
from collections import Counter
from fractions import Fraction
from pathlib import Path

Pair = tuple[Fraction, Fraction]  # x+y sqrt(D)
ZERO: Pair = (Fraction(0), Fraction(0))
ONE: Pair = (Fraction(1), Fraction(0))


def add(x: Pair, y: Pair) -> Pair:
    return (x[0] + y[0], x[1] + y[1])


def sub(x: Pair, y: Pair) -> Pair:
    return (x[0] - y[0], x[1] - y[1])


def mul(x: Pair, y: Pair, D: int) -> Pair:
    return (x[0] * y[0] + D * x[1] * y[1], x[0] * y[1] + x[1] * y[0])


def scale(n: int, x: Pair) -> Pair:
    return (n * x[0], n * x[1])


def conjugate(x: Pair) -> Pair:
    return (x[0], -x[1])


def sign(x: Pair, D: int) -> int:
    """Exact sign of x0+x1 sqrt(D); only rational arithmetic is used."""
    p, q = x
    if q == 0:
        return (p > 0) - (p < 0)
    if p >= 0 and q > 0:
        return 1
    if p <= 0 and q < 0:
        return -1
    delta = p * p - D * q * q
    if p > 0:
        return (delta > 0) - (delta < 0)
    return (delta < 0) - (delta > 0)


def floor_pair(x: Pair, D: int) -> int:
    lo, hi = -1, 1
    while sign(sub(x, (Fraction(lo), Fraction(0))), D) < 0:
        hi, lo = lo, 2 * lo
    while sign(sub(x, (Fraction(hi), Fraction(0))), D) >= 0:
        lo, hi = hi, 2 * hi
    while hi - lo > 1:
        mid = (lo + hi) // 2
        if sign(sub(x, (Fraction(mid), Fraction(0))), D) >= 0:
            lo = mid
        else:
            hi = mid
    return lo


def omega(m: int, n: int, D: int) -> Pair:
    return ((Fraction(2 * m + n, 2), Fraction(n, 2)) if D % 4 == 1
            else (Fraction(m), Fraction(n)))


def pair_text(x: Pair) -> list[str]:
    return [str(x[0]), str(x[1])]


def candidates(base: Pair, D: int, forced: int | None) -> list[tuple[int, Pair]]:
    """All integer g making base+g lie in [0,1], respecting a force."""
    if forced is not None:
        z = add(base, (Fraction(forced), Fraction(0)))
        return [(forced, z)] if sign(z, D) >= 0 and sign(sub(ONE, z), D) >= 0 else []
    if base[1] == 0:
        assert base[0].denominator == 1, base
        b = int(base[0])
        return [(j, add(base, (Fraction(j), Fraction(0)))) for j in (-b, 1 - b)]
    j = -floor_pair(base, D)
    z = add(base, (Fraction(j), Fraction(0)))
    assert sign(z, D) >= 0 and sign(sub(ONE, z), D) > 0
    return [(j, z)]


def unit_candidates(base: Pair, D: int) -> list[tuple[int, Pair]]:
    """All translates putting base in [0,1], before an endpoint force."""
    return candidates(base, D, None)


def canonical_state_hash(states: list[tuple[tuple[int, ...], tuple[int, ...], tuple[Pair, ...], tuple[Pair, ...]]]) -> str:
    """Same G/Q-only canonical serialization as the independent HGQ tree."""
    h = hashlib.sha256()
    for G, Q, _A, _B in sorted(states, key=lambda s: (s[0], s[1])):
        h.update(json.dumps([G, Q], separators=(",", ":")).encode() + b"\n")
    return h.hexdigest()


def coeff_base(P: list[Pair], prior: tuple[int, ...], k: int, D: int) -> Pair:
    out = ZERO
    for i in range(1, min(k, len(P) - 1) + 1):
        out = add(out, scale(prior[k - i], P[i]))
    return out


def c_coefficient(A: tuple[Pair, ...], B: tuple[Pair, ...], k: int, D: int) -> int | None:
    z = ZERO
    for i in range(k + 1):
        z = add(z, mul(A[i], B[k - i], D))
    if z[1] != 0 or z[0].denominator != 1:
        return None
    return int(z[0]) if z[0] in (0, 1) else None


def run_core(row: list[int], depth: int, cap: int, seconds: float | None = None) -> dict:
    """Return a compact exact transcript for one [D,L,am,an,k,dm,dn,e,f] row."""
    D, L, am, an, k0, dm, dn, e, f = row
    a, d = omega(am, an, D), omega(dm, dn, D)
    P: list[Pair] = [ONE, (Fraction(L), Fraction(0)), a, add(scale(L,a),(Fraction(k0),Fraction(0))), d,
                     (Fraction(e), Fraction(0)), (Fraction(f), Fraction(0)), ONE]
    Pb = [conjugate(x) for x in P]
    # state holds all integer cofactor terms and A/B through its current depth.
    states: list[tuple[tuple[int, ...], tuple[int, ...], tuple[Pair, ...], tuple[Pair, ...]]] = [
        ((1,), (1,), (ONE,), (ONE,))
    ]
    levels = []
    rejections: Counter[str] = Counter()
    deadline = None if seconds is None else time.monotonic() + seconds
    for k in range(1, depth + 1):
        next_states = []
        forced = k in (1,3)
        level_rejected: Counter[str] = Counter(forced_endpoint=0, nonbinary_C=0)
        for G, Q, A, B in states:
            if deadline is not None and time.monotonic() > deadline:
                return {
                    "status": "unknown_time_limit", "row": row, "depth_reached": k - 1,
                    "levels": levels, "live": len(states), "rejections": dict(rejections),
                }
            base_a = coeff_base(P, G, k, D)
            base_b = coeff_base(Pb, Q, k, D)
            choices_a = unit_candidates(base_a, D)
            choices_b = unit_candidates(base_b, D)
            for gk, Ak in choices_a:
                for qk, Bk in choices_b:
                    if forced and (Ak != ZERO or Bk != ZERO):
                        level_rejected["forced_endpoint"] += 1
                        rejections["forced_endpoint"] += 1
                        continue
                    ck = c_coefficient(A + (Ak,), B + (Bk,), k, D)
                    if ck is None:
                        level_rejected["nonbinary_C"] += 1
                        rejections["nonbinary_C"] += 1
                        continue
                    next_states.append((G + (gk,), Q + (qk,), A + (Ak,), B + (Bk,)))
                    if len(next_states) > cap:
                        return {
                            "status": "unknown_state_cap", "row": row, "depth_reached": k,
                            "levels": levels, "live": len(next_states), "rejections": dict(rejections),
                        }
        states = next_states
        digest = canonical_state_hash(states)
        levels.append({"k": k, "live": len(states), "state_sha256": digest,
                       "rejected": dict(level_rejected)})
        if not states:
            return {"status": "dead", "row": row, "depth_reached": k,
                    "levels": levels, "live": 0, "rejections": dict(rejections)}
    sample = []
    for G, Q, A, B in states[:3]:
        sample.append({"G": list(G), "Q": list(Q),
                       "A": [pair_text(z) for z in A], "B": [pair_text(z) for z in B]})
    return {"status": "unknown_depth", "row": row, "depth_reached": depth,
            "levels": levels, "live": len(states), "rejections": dict(rejections), "sample_states": sample}


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


