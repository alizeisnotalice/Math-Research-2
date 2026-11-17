#!/usr/bin/env python3
"""Exact independent-peak IE coefficients h_1,h_2,h_3 on punctured grids.

This is deliberately distinct from R028's prefix-tail DP.  It evaluates the
proposed fixed-prefix formulas using C_H(I), then compares with its exact L
histograms for every hole on n=3,4.
"""
from __future__ import annotations

import argparse
import json
from collections import defaultdict
from math import factorial
from pathlib import Path


def graph(n: int):
    N = n * n
    adj = [0] * N
    for y in range(n):
        for x in range(n):
            v = y * n + x
            for dx, dy in ((1, 0), (-1, 0), (0, 1), (0, -1)):
                xx, yy = x + dx, y + dy
                if 0 <= xx < n and 0 <= yy < n:
                    adj[v] |= 1 << (yy * n + xx)
    return adj


def independent_bags(n: int, allowed: int):
    rows = [r for r in range(1 << n) if not r & (r << 1)]
    allowed_rows = [(allowed >> (n * y)) & ((1 << n) - 1) for y in range(n)]
    bags = defaultdict(list)

    def walk(y: int, prev: int, I: int):
        if y == n:
            bags[I.bit_count()].append(I)
            return
        for r in rows:
            if not r & prev and not r & ~allowed_rows[y]:
                walk(y + 1, r, I | (r << (n * y)))

    walk(0, 0, 0)
    return bags


def closure(mask: int, adj: list[int]) -> int:
    ans, q = mask, mask
    while q:
        b = q & -q
        q ^= b
        ans |= adj[b.bit_length() - 1]
    return ans


def falling(M: int, r: int) -> int:
    ans = 1
    for j in range(r):
        ans *= M - j
    return ans


def directed_edges(P: int, U: int, adj: list[int]) -> int:
    ans, q = 0, P
    while q:
        b = q & -q
        q ^= b
        ans += (adj[b.bit_length() - 1] & U).bit_count()
    return ans


def eta_words(length: int, U: int, Nv: int, adj: list[int]):
    """Words eta grown from N(v), excluding a root p by later closure test."""
    if length == 0:
        yield 0
        return

    def rec(S: int, used: int, left: int):
        if left == 0:
            yield S
            return
        avail = U & ~used & (Nv | closure(S, adj))
        # For the first point closure(empty)=empty, so this is exactly Nv.
        q = avail
        while q:
            b = q & -q
            q ^= b
            yield from rec(S | b, used | b, left - 1)

    # A word's order matters.  The set S is enough for admissibility only
    # after the order is recursively generated; duplicate final masks from
    # distinct valid orders are deliberately yielded separately.
    yield from rec(0, 0, length)


def c_h_values(n: int, v: int, adj: list[int]):
    N, M = n * n, n * n - 1
    full = (1 << N) - 1
    H = full ^ (1 << v)
    Nv = adj[v]
    Z = H & ~Nv
    bags = independent_bags(n, Z)
    C = {0: factorial(M)}
    cls = {0: 0}
    for k in range(1, max(bags) + 1):
        for I in bags[k]:
            cl = closure(I, adj)
            # I lies in Z, hence v is not in cl.  This is cl_H(I).
            assert not (cl >> v) & 1
            numerator = sum(C[I ^ b] for b in (1 << j for j in range(N) if (I >> j) & 1))
            assert numerator % cl.bit_count() == 0
            C[I] = numerator // cl.bit_count()
            cls[I] = cl
    return H, Nv, Z, bags, C, cls


def h_for_v(n: int, v: int, max_t: int = 3):
    N, M = n * n, n * n - 1
    full = (1 << N) - 1
    adj = graph(n)
    H, Nv, Z, bags, C, cls = c_h_values(n, v, adj)
    nums = [0] * (max_t + 1)
    for I, count in C.items():
        k, clI = I.bit_count(), cls[I]
        U = H & ~clI
        P = Z & U
        sign = -1 if k & 1 else 1
        for t in range(1, max_t + 1):
            den = falling(M, t + 1)
            subtotal = 0
            for eta in eta_words(t - 1, U, Nv, adj):
                # P\closed(eta) enforces that root p itself and every eta
                # entry are distinct and non-adjacent: no p-neighbour occurs
                # before the terminal x.
                subtotal += directed_edges(P & ~closure(eta, adj), U, adj)
            if subtotal:
                # Any contributing endpoint realizes a fixed ordered prefix
                # of t+1 vertices in U, so its deletion/insertion fibre gives
                # this separate divisibility.  If no prefix exists (e.g. U
                # has fewer than t+1 vertices), no such division is used.
                assert count % den == 0, (n, v, I, t, count, den)
                nums[t] += sign * (count // den) * subtotal
    h = []
    for t in range(1, max_t + 1):
        h.append(nums[t])
    return {"v": v, "xy": [v % n, v // n], "h": h,
            "ie_signed_sums_after_prefix_division": nums[1:],
            "independent_sets": sum(map(len, bags.values()))}


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--tail-dir", type=Path, required=True)
    p.add_argument("--output", type=Path, required=True)
    a = p.parse_args()
    out = {"status": "passed", "checks": []}
    for n in (3, 4):
        tail = json.loads((a.tail_dir / f"tail-n{n}.json").read_text())
        got = [h_for_v(n, v) for v in range(n * n)]
        for x, ref in zip(got, tail["vertices"]):
            assert x["h"] == ref["histogram_L"][:3], (n, x["v"], x["h"], ref["histogram_L"][:3])
        out["checks"].append({"n": n, "all_vertices": True, "values": got})
        print(f"n={n}: all {n*n} punctures h1,h2,h3 matched R028 tail DP")
    a.output.write_text(json.dumps(out, indent=2) + "\n")


if __name__ == "__main__":
    main()
