#!/usr/bin/env python3
"""Exact finite search for a counterexample to the off-p overlap lemma.

For R=Q(x^p), divisibility R|U splits independently into p residue rows:
each binary row polynomial is divisible by Q.  This program enumerates row
polynomials exactly, combines them into binary U,V, then requires U*V and
H=(U*V)/R to be binary.  A printed hit is independently rechecked with
integer convolution and long division.  An empty run is only finite evidence.
"""
from __future__ import annotations
from itertools import product
from pathlib import Path
from time import perf_counter
import argparse
import json

HERE = Path(__file__).resolve().parent

QS = {
    "1-y+y3": [1, -1, 0, 1],
    "1-y2+y3": [1, 0, -1, 1],
    "1+y2+y3": [1, 0, 1, 1],
    "1+y+y3": [1, 1, 0, 1],
    "Phi6": [1, -1, 1],
}

def trim(a):
    a = a[:]
    while len(a) > 1 and a[-1] == 0:
        a.pop()
    return a

def mul(a, b):
    c = [0] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        if x:
            for j, y in enumerate(b):
                if y:
                    c[i+j] += x*y
    return trim(c)

def divide_exact(a, d):
    a = a[:]
    while len(a) >= len(d):
        q = a[-1]  # d is monic
        if q:
            j = len(a) - len(d)
            for i, x in enumerate(d):
                a[i+j] -= q*x
        a.pop()
    return a == []

def quotient_exact(a, d):
    a = a[:]
    if a == [0]:
        return [0]
    q = [0] * max(1, len(a) - len(d) + 1)
    while len(a) >= len(d):
        z = a[-1]
        j = len(a) - len(d)
        q[j] = z
        for i, x in enumerate(d):
            a[i+j] -= z*x
        a = trim(a)
    if a != [0]:
        return None
    return trim(q)

def row_options(q, maxdeg, forced_constant=None):
    """All binary W of degree <=maxdeg with Q|W, with quotient included."""
    ans = []
    for mask in range(1 << (maxdeg + 1)):
        if forced_constant is not None and (mask & 1) != forced_constant:
            continue
        w = [(mask >> i) & 1 for i in range(maxdeg + 1)]
        z = quotient_exact(trim(w), q)
        if z is not None:
            ans.append((tuple(w), tuple(z)))
    return ans

def combine_rows(rows, p, cap):
    u = [0] * (cap + 1)
    e = [0] * (cap + 1)
    for r, (w, z) in enumerate(rows):
        for j, v in enumerate(w):
            if r + p*j <= cap:
                u[r+p*j] = v
        for j, v in enumerate(z):
            if r + p*j <= cap:
                e[r+p*j] = v
    return trim(u), trim(e)

def is_binary(a):
    return all(x in (0, 1) for x in a)

def support(a):
    return {i for i, v in enumerate(a) if v}

def unique_sum(sa, sb):
    seen = set()
    for i in sa:
        for j in sb:
            k = i+j
            if k in seen:
                return False
            seen.add(k)
    return True

def candidates(p, q, cap):
    opts = []
    for r in range(p):
        d = (cap-r)//p
        opts.append(row_options(q, d, forced_constant=1 if r == 0 else None))
    out = []
    for rows in product(*opts):
        u, e = combine_rows(rows, p, cap)
        if u[0] != 1 or not is_binary(u):
            raise AssertionError("row construction error")
        out.append((u, e))
    return out, [len(x) for x in opts]

def run_one(p, qname, q, cap, hit_limit):
    r = [0] * (p*(len(q)-1)+1)
    for j, x in enumerate(q):
        r[p*j] = x
    us, row_counts = candidates(p, q, cap)
    t0 = perf_counter()
    tested_pairs = 0
    unique_pairs = 0
    hits = []
    for i, (u, e) in enumerate(us):
        su = support(u)
        for v, d in us[i:]:
            tested_pairs += 1
            if not unique_sum(su, support(v)):
                continue
            unique_pairs += 1
            uv = mul(u, v)
            assert is_binary(uv)
            h = quotient_exact(uv, r)
            if h is None or not is_binary(h):
                continue
            for n in range(1, min(len(e), len(d))):
                if n % p and e[n] == d[n] == 1:
                    # Full independent integer recheck before recording.
                    assert mul(r, e) == u and mul(r, d) == v
                    assert mul(r, h) == uv and is_binary(u) and is_binary(v)
                    hits.append({"n": n, "U": u, "V": v, "E": e, "D": d, "H": h})
                    if len(hits) >= hit_limit:
                        break
            if len(hits) >= hit_limit:
                break
        if len(hits) >= hit_limit:
            break
    return {
        "p": p, "Q_name": qname, "Q_ascending": q, "degree_cap_U_V": cap,
        "row_option_counts": row_counts, "U_candidate_count": len(us),
        "tested_unordered_pairs": tested_pairs, "unique_sum_pairs": unique_pairs,
        "elapsed_seconds": perf_counter()-t0, "hits": hits,
        "claim_boundary": "Exact finite enumeration only; an empty hit list is not an exclusion beyond this cap."
    }

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--p", type=int, default=2)
    ap.add_argument("--cap", type=int, default=24)
    ap.add_argument("--q", choices=sorted(QS), default=None)
    ap.add_argument("--hit-limit", type=int, default=3)
    ap.add_argument("--output", type=Path, default=HERE/"offp_overlap_search_output.json")
    args = ap.parse_args()
    todo = [(args.q, QS[args.q])] if args.q else list(QS.items())
    reports = [run_one(args.p, name, q, args.cap, args.hit_limit) for name, q in todo]
    args.output.write_text(json.dumps({"reports": reports}, indent=2)+"\n")
    for x in reports:
        print(x["Q_name"], "rows", x["row_option_counts"], "U", x["U_candidate_count"],
              "pairs", x["tested_unordered_pairs"], "unique", x["unique_sum_pairs"],
              "hits", len(x["hits"]), f"{x['elapsed_seconds']:.3f}s")

if __name__ == "__main__":
    main()
