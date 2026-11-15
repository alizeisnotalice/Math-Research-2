#!/usr/bin/env python3
"""Exact Bernstein branch-and-prune certificate for the k=4,6 tail guards.

All arithmetic is fractions.Fraction.  A leaf is discarded only by an upper
Bernstein bound < 0 for a required >=0 guard, or a lower Bernstein bound >0
for F=h^2 q, where h>0.  The generated JSON contains the full binary cover.
"""
from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction as Q
from math import comb
import json
from pathlib import Path
from typing import Dict, Tuple, List

Mon = Tuple[int, int]
Poly = Dict[Mon, Q]


def clean(p: Poly) -> Poly:
    return {m: c for m, c in p.items() if c}


def add(p: Poly, q: Poly) -> Poly:
    out = dict(p)
    for m, c in q.items():
        out[m] = out.get(m, Q(0)) + c
    return clean(out)


def neg(p: Poly) -> Poly:
    return {m: -c for m, c in p.items()}


def sub(p: Poly, q: Poly) -> Poly:
    return add(p, neg(q))


def scale(c: Q, p: Poly) -> Poly:
    return clean({m: c * v for m, v in p.items()})


def mul(p: Poly, q: Poly) -> Poly:
    out: Poly = {}
    for (i, j), a in p.items():
        for (r, t), b in q.items():
            m = (i + r, j + t)
            out[m] = out.get(m, Q(0)) + a * b
    return clean(out)


def power(p: Poly, n: int) -> Poly:
    out: Poly = {(0, 0): Q(1)}
    base = p
    while n:
        if n & 1:
            out = mul(out, base)
        base = mul(base, base)
        n //= 2
    return out


ONE: Poly = {(0, 0): Q(1)}
S: Poly = {(1, 0): Q(1)}
U: Poly = {(0, 1): Q(1)}


def degree(p: Poly) -> Tuple[int, int]:
    return max(i for i, _ in p), max(j for _, j in p)


def affine_power(a: Q, d: Q, n: int) -> List[Q]:
    return [Q(comb(n, r)) * a ** (n - r) * d ** r for r in range(n + 1)]


def local_power(p: Poly, slo: Q, shi: Q, ulo: Q, uhi: Q) -> Poly:
    ds, du = shi - slo, uhi - ulo
    out: Poly = {}
    for (i, j), coef in p.items():
        xs, xu = affine_power(slo, ds, i), affine_power(ulo, du, j)
        for r, cr in enumerate(xs):
            for t, ct in enumerate(xu):
                m = (r, t)
                out[m] = out.get(m, Q(0)) + coef * cr * ct
    return clean(out)


def bernstein_bounds(p: Poly, box: Tuple[Q, Q, Q, Q]) -> Tuple[Q, Q]:
    n, m = degree(p)
    loc = local_power(p, *box)
    vals: List[Q] = []
    for r in range(n + 1):
        for t in range(m + 1):
            v = Q(0)
            for (i, j), c in loc.items():
                if i <= r and j <= t:
                    v += c * Q(comb(r, i), comb(n, i)) * Q(comb(t, j), comb(m, j))
            vals.append(v)
    return min(vals), max(vals)


def fmt(x: Q) -> str:
    return str(x.numerator) if x.denominator == 1 else f"{x.numerator}/{x.denominator}"


def build_guards(k: int, r: int) -> Tuple[Dict[str, Poly], Poly]:
    h = add(ONE, S)
    a = sub(scale(Q(2), S), ONE)
    y = power(S, k + 1)
    def v(j: int) -> Poly:
        return add(ONE, add(scale(Q(j + 2), power(S, j + 1)), scale(Q(j + 1), power(S, j + 2))))
    V = v(k)
    Vminus = v(k - 2)
    Vplus = v(k + 2)
    E = r + k + 1
    g1 = add(sub(power(S, E), ONE), U)
    g2 = sub(add(power(S, E), ONE), U)
    # Multiply g3 by h>0 and fm/fp by h^2>0, after W_j=V_j/h^2.
    g3 = sub(V, mul(U, h))
    ss = mul(S, S)
    h2 = mul(h, h)
    fm = add(mul(sub(ss, ONE), V), mul(U, sub(V, mul(ss, Vminus))))
    fp = sub(mul(sub(ONE, ss), V), mul(U, sub(Vplus, mul(ss, V))))
    # F=h^2 q, a positive rescaling of q on the whole box.
    F = add(
        sub(mul(mul(S, a), mul(V, V)), mul(mul(mul(a, h), V), U)),
        mul(mul(sub(a, scale(Q(3), mul(ss, y))), h), mul(U, U)),
    )
    return {"g1": g1, "g2": g2, "g3h": g3, "fm_h2": fm, "fp_h2": fp}, F


@dataclass
class Node:
    box: Tuple[Q, Q, Q, Q]
    depth: int


def make_tree(k: int, r: int, max_depth: int = 42) -> dict:
    guards, F = build_guards(k, r)
    delta = Q(1, 4 * k * (k + 1))
    root = Node((delta, Q(1), Q(0), Q(2)), 0)
    leaves: List[dict] = []
    stack: List[Node] = [root]
    internal = 0
    max_seen = 0
    while stack:
        node = stack.pop()
        max_seen = max(max_seen, node.depth)
        bounds = {name: bernstein_bounds(poly, node.box) for name, poly in guards.items()}
        fb = bernstein_bounds(F, node.box)
        reason = None
        for name, (_, hi) in bounds.items():
            if hi < 0:
                reason = {"kind": "guard_upper_negative", "guard": name, "bound": fmt(hi)}
                break
        if reason is None and fb[0] > 0:
            reason = {"kind": "F_lower_positive", "bound": fmt(fb[0])}
        rec = {"box": [fmt(x) for x in node.box], "depth": node.depth}
        if reason is not None:
            rec["leaf"] = reason
            leaves.append(rec)
            continue
        if node.depth >= max_depth:
            raise RuntimeError(f"unresolved box at depth {node.depth}: {rec}; F={tuple(map(fmt,fb))}; guards={{{', '.join(f'{n}:{tuple(map(fmt,b))}' for n,b in bounds.items())}}}")
        slo, shi, ulo, uhi = node.box
        # Split the larger normalized side; exact midpoint preserves a binary cover.
        s_norm = (shi - slo) / (Q(1) - delta)
        u_norm = (uhi - ulo) / Q(2)
        if u_norm >= s_norm:
            mid = (ulo + uhi) / 2
            stack.append(Node((slo, shi, mid, uhi), node.depth + 1))
            stack.append(Node((slo, shi, ulo, mid), node.depth + 1))
            rec["split"] = "u"
        else:
            mid = (slo + shi) / 2
            stack.append(Node((mid, shi, ulo, uhi), node.depth + 1))
            stack.append(Node((slo, mid, ulo, uhi), node.depth + 1))
            rec["split"] = "s"
        internal += 1
    return {
        "format": "rational-bernstein-cover-v1",
        "k": k,
        "r": r,
        "delta": fmt(delta),
        "root": [fmt(delta), "1", "0", "2"],
        "polynomial_degrees": {**{n: degree(p) for n, p in guards.items()}, "F_h2q": degree(F)},
        "internal_count": internal,
        "leaf_count": len(leaves),
        "max_depth": max_seen,
        "leaves": leaves,
    }


def main() -> None:
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", type=Path, required=True)
    outdir = parser.parse_args().out
    assert Path("/Volumes/data").is_mount()
    outdir.mkdir(parents=True, exist_ok=True)
    for k, r in ((4, 7), (6, 9)):
        cert = make_tree(k, r)
        path = outdir / f"bernstein_tail_k{k}_r{r}.json"
        path.write_text(json.dumps(cert, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
        print(f"k={k} r={r} leaves={cert['leaf_count']} internal={cert['internal_count']} max_depth={cert['max_depth']} path={path}")


if __name__ == "__main__":
    main()
