#!/usr/bin/env python3
"""Exact necessary-prefix enumeration for P=1+a*x+c*x^2+d*x^3+e*x^4+f*x^5+x^6.

All field elements use the integer basis p+q*a with a^2=T*a-N.
The actual embedding is selected by orient=+1/-1; every comparison is reduced
to an integer square comparison with sqrt(D), never a floating decision.
"""
from collections import Counter
from pathlib import Path
import json, math

OUT = Path(__file__).with_name("core6_single_results.json")
ZERO, ONE = (0, 0), (1, 0)


def add(x, y): return (x[0] + y[0], x[1] + y[1])
def sub(x, y): return (x[0] - y[0], x[1] - y[1])
def scale(k, x): return (k * x[0], k * x[1])


def mul(x, y, T, N):
    p, q = x; r, s = y
    return (p*r - N*q*s, p*s + q*r + T*q*s)


def sign(x, T, D, orient):
    """Sign of p+q*a where a=(T+orient*sqrt(D))/2."""
    p, q = x
    A, B = 2*p + q*T, q*orient
    if B == 0: return (A > 0) - (A < 0)
    if A == 0: return (B > 0) - (B < 0)
    if (A > 0) == (B > 0): return 1 if A > 0 else -1
    diff = A*A - D*B*B
    assert diff != 0  # D nonsquare and B nonzero.
    if A > 0 and B < 0: return 1 if diff > 0 else -1
    return 1 if diff < 0 else -1


def floor_expr(x, T, D, orient):
    p, q = x
    guess = math.floor(p + q * (T + orient * math.sqrt(D)) / 2)
    while sign(sub(x, (guess, 0)), T, D, orient) < 0:
        guess -= 1
    while sign(sub(x, (guess + 1, 0)), T, D, orient) >= 0:
        guess += 1
    return guess


def in_unit(x, T, D, orient):
    return sign(x, T, D, orient) >= 0 and sign(sub(ONE, x), T, D, orient) >= 0


def bit(x): return x == ZERO or x == ONE


def eligible_tn():
    for T in range(-14, 15):
        for N in range(-56, 57):
            D = T*T - 4*N
            if D > 0 and math.isqrt(D)**2 != D and D < (15 - abs(T))**2:
                yield T, N, D


def e_choices(K, T, D, orient):
    center = -floor_expr(K, T, D, orient)
    return {e for e in (center - 1, center, center + 1)
            if -21 <= e <= 21 and in_unit(add(K, (e, 0)), T, D, orient)}


def conv_at(A, B, k, T, N):
    total = ZERO
    for i in range(k + 1):
        total = add(total, mul(A[i], B[k-i], T, N))
    return total


def main():
    pairs = list(eligible_tn())
    counts = Counter(tn_pairs=len(pairs), oriented_roots=2*len(pairs))
    candidates = []
    for T, N, D in pairs:
        for orient in (-1, 1):
            aa, abar = (0, 1), (T, -1)
            g1 = -floor_expr(aa, T, D, orient)
            h1 = -floor_expr(abar, T, D, orient)
            A1, B1 = add(aa, (g1, 0)), add(abar, (h1, 0))
            if not (in_unit(A1, T, D, orient) and in_unit(B1, T, D, orient) and bit(add(A1, B1))):
                continue
            counts["pass_c1"] += 1
            for c in range(-21, 22):
                g2 = -floor_expr(add((c, 0), scale(g1, aa)), T, D, orient)
                h2 = -floor_expr(add((c, 0), scale(h1, abar)), T, D, orient)
                A2 = add(add((c, 0), scale(g1, aa)), (g2, 0))
                B2 = add(add((c, 0), scale(h1, abar)), (h2, 0))
                if not (in_unit(A2, T, D, orient) and in_unit(B2, T, D, orient) and
                        bit(add(add(A2, B2), mul(A1, B1, T, N)))):
                    continue
                counts["pass_c2"] += 1
                for d in range(-30, 31):
                    g3 = -floor_expr(add(add((d, 0), scale(g1, (c, 0))), scale(g2, aa)), T, D, orient)
                    h3 = -floor_expr(add(add((d, 0), scale(h1, (c, 0))), scale(h2, abar)), T, D, orient)
                    A3 = add(add(add((d, 0), scale(g1, (c, 0))), scale(g2, aa)), (g3, 0))
                    B3 = add(add(add((d, 0), scale(h1, (c, 0))), scale(h2, abar)), (h3, 0))
                    if not (in_unit(A3, T, D, orient) and in_unit(B3, T, D, orient)):
                        continue
                    Apre, Bpre = [ONE, A1, A2, A3], [ONE, B1, B2, B3]
                    if not bit(conv_at(Apre, Bpre, 3, T, N)):
                        continue
                    counts["pass_c3"] += 1
                    KA = add(add(scale(g3, aa), (g2*c, 0)), (g1*d, 0))
                    KB = add(add(scale(h3, abar), (h2*c, 0)), (h1*d, 0))
                    for e in sorted(e_choices(KA, T, D, orient) & e_choices(KB, T, D, orient)):
                        A4, B4 = add(KA, (e, 0)), add(KB, (e, 0))
                        Apre4, Bpre4 = Apre + [A4], Bpre + [B4]
                        if not bit(conv_at(Apre4, Bpre4, 4, T, N)):
                            continue
                        counts["pass_c4"] += 1
                        c5 = conv_at(Apre4 + [ZERO], Bpre4 + [ZERO], 5, T, N)
                        if not bit(c5):
                            continue
                        counts["pass_c5"] += 1
                        for f in range(-7, 8):
                            g5 = -(f + g1*e + g2*d + g3*c)
                            h5 = -(f + h1*e + h2*d + h3*c)
                            candidates.append({
                                "T": T, "N": N, "D": D, "orientation": orient,
                                "c": c, "d": d, "e": e, "f": f,
                                "G1_to_G5": [g1, g2, g3, 0, g5],
                                "Q1_to_Q5": [h1, h2, h3, 0, h5],
                            })
    payload = {
        "scope": "Necessary low/high prefix filter only; no existence claim for full A,B.",
        "bounds": {"T": [-14, 14], "N": [-56, 56], "c_e": [-21, 21], "d": [-30, 30], "f": [-7, 7]},
        "counts": dict(counts), "candidate_count": len(candidates), "candidates": candidates,
    }
    OUT.write_text(json.dumps(payload, indent=2) + "\n")
    print("PASS exact core6 enumeration counts:", dict(counts))
    print("PASS exact core6 candidate count:", len(candidates), "artifact:", OUT)


if __name__ == "__main__":
    main()
