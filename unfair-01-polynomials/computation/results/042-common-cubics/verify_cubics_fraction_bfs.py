#!/usr/bin/env python3
"""Pure-Python exact audit of the round-42 cubic candidate/BFS certificate.

No Wolfram code, WXF payload, pilot state list, or reported Boolean is read.
Uses only integer/Fraction arithmetic: a Sturm chain for the coefficient box
and isolated-root Fraction intervals for all automaton boundary decisions.
"""
from fractions import Fraction as F
from hashlib import sha256
from pathlib import Path
import json


def trim(a):
    a = list(a)
    while len(a) > 1 and a[-1] == 0:
        a.pop()
    return a


def derivative(a):
    return trim([i * a[i] for i in range(1, len(a))] or [F(0)])


def divrem(a, b):
    a, b = trim(a), trim(b)
    q = [F(0)] * max(1, len(a) - len(b) + 1)
    while len(a) >= len(b) and a != [0]:
        c = a[-1] / b[-1]
        d = len(a) - len(b)
        q[d] += c
        for i, u in enumerate(b):
            a[d + i] -= c * u
        a = trim(a)
    return trim(q), trim(a)


def rem_monic(a, b):
    return divrem([F(x) for x in a], [F(x) for x in b])[1]


def sturm(a):
    z = [trim([F(x) for x in a]), derivative([F(x) for x in a])]
    while z[-1] != [0]:
        r = divrem(z[-2], z[-1])[1]
        if r == [0]:
            break
        z.append(trim([-u for u in r]))
    return z


def peval(a, x):
    out = F(0)
    for u in reversed(a):
        out = out * x + u
    return out


def sign_at(a, x):
    z = peval(a, x)
    return (z > 0) - (z < 0)


def sign_infty(a, positive):
    s = (a[-1] > 0) - (a[-1] < 0)
    return s if positive or (len(a) - 1) % 2 == 0 else -s


def variations(signs):
    signs = [s for s in signs if s]
    return sum(a != b for a, b in zip(signs, signs[1:]))


def roots_between(chain, left, right):
    sl = [sign_infty(a, False) if left == "-inf" else sign_at(a, left) for a in chain]
    sr = [sign_infty(a, True) if right == "inf" else sign_at(a, right) for a in chain]
    return variations(sl) - variations(sr)


def iadd(I, J):
    return I[0] + J[0], I[1] + J[1]


def imul(I, J):
    z = (I[0]*J[0], I[0]*J[1], I[1]*J[0], I[1]*J[1])
    return min(z), max(z)


def pinterval(a, I):
    out = (F(a[-1]), F(a[-1]))
    for u in reversed(a[:-1]):
        out = iadd(imul(out, I), (F(u), F(u)))
    return out


def fmt(x):
    return str(x.numerator) if x.denominator == 1 else f"{x.numerator}/{x.denominator}"


def main():
    # Candidate enumeration: monic cubic, constant 1. Reducible iff ±1 root.
    candidates = []
    for c1 in range(-5, 6):
        for c2 in range(-5, 6):
            q = [F(1), F(c1), F(c2), F(1)]
            irreducible = peval(q, F(1)) != 0 and peval(q, F(-1)) != 0
            chain = sturm(q)
            positive = roots_between(chain, F(0), "inf")
            real = roots_between(chain, "-inf", "inf")
            in_band = roots_between(chain, F(-2), F(-1, 2))
            if irreducible and positive == 0 and in_band == real and real in (1, 3):
                candidates.append((c1, c2, real))
    expected = [(-1, 0, 1), (-1, 1, 1), (0, -1, 1), (0, 1, 1),
                (1, -1, 1), (1, 0, 1), (1, 2, 1), (2, 1, 1)]
    assert candidates == expected, candidates

    # One representative from each reciprocal pair, with rational isolation
    # interval, B and D (not their squares) deliberately wider than the tail.
    jobs = [
        ("m10", -1, 0, F(-4, 3), F(-29, 22), F(16, 5), F(10)),
        ("m11", -1, 1, F(-2), F(-11, 6), F(5, 4), F(3, 2)),
        ("01", 0, 1, F(-3, 2), F(-10, 7), F(3), F(6)),
        ("12", 1, 2, F(-2), F(-7, 4), F(4, 3), F(2)),
    ]
    results = []
    for tag, c1, c2, lo, hi, B, D in jobs:
        q = [F(1), F(c1), F(c2), F(1)]
        q2 = [F(0)] * 7
        for i, u in enumerate(q):
            for j, v in enumerate(q):
                q2[i+j] += u*v
        assert peval(q, lo) < 0 < peval(q, hi)
        # q' positive throughout each declared interval; interval Horner certifies.
        qprimeI = pinterval(derivative(q), (lo, hi))
        assert qprimeI[0] > 0
        # Narrow root isolate by exact bisection, retaining q(lo)<0<q(hi).
        for _ in range(192):
            mid = (lo + hi) / 2
            if peval(q, mid) < 0:
                lo = mid
            else:
                hi = mid
        I = (lo, hi)
        # Tail bounds are strict because |rho| > -declared right endpoint.
        # The bisection right endpoint remains strictly to the right of the
        # negative root, so -hi is a strict lower bound for |rho|.  This is
        # needed in the m12 row, where the coarse endpoint only gives equality.
        Rmin = -hi
        assert F(1, 1)/(Rmin - 1) < B
        assert F(1, 1)/(Rmin - 1)**2 < D

        def step(a, bit):
            return tuple(rem_monic([F(bit)] + [F(u) for u in a], q2) + [F(0)] * 6)[:6]

        def status(a):
            # r(rho), r'(rho) are evaluated by r mod Q and (r') mod Q.
            v = rem_monic(a, q)
            dv = rem_monic(derivative(a), q)
            VI, DI = pinterval(v, I), pinterval(dv, I)
            if VI[0] >= -B and VI[1] <= B and DI[0] >= -D and DI[1] <= D:
                return True, {"v": VI, "dv": DI, "reason": "inside"}
            if VI[0] > B or VI[1] < -B:
                return False, {"v": VI, "dv": DI, "reason": "value"}
            if DI[0] > D or DI[1] < -D:
                return False, {"v": VI, "dv": DI, "reason": "derivative"}
            raise AssertionError((tag, a, v, dv, VI, DI, B, D))

        start = (F(1), F(0), F(0), F(0), F(0), F(0))
        ok, _ = status(start)
        assert ok
        states, seen, head = [start], {start}, 0
        boundary = []
        while head < len(states):
            a = states[head]
            head += 1
            for bit in (0, 1):
                z = step(a, bit)
                ok, cert = status(z)
                if ok:
                    if z not in seen:
                        seen.add(z); states.append(z)
                else:
                    boundary.append({"from": [int(v) for v in a], "bit": bit,
                                     "target": [int(v) for v in z],
                                     "reason": cert["reason"],
                                     "v_interval": [fmt(cert["v"][0]), fmt(cert["v"][1])],
                                     "dv_interval": [fmt(cert["dv"][0]), fmt(cert["dv"][1])]})
        assert head == len(states) and (F(0),)*6 not in seen
        results.append({"tag": tag, "c1": c1, "c2": c2,
                        "root_interval_after_192_bisections": [fmt(lo), fmt(hi)],
                        "declared_root_interval": [fmt(jobs[[z[0] for z in jobs].index(tag)][3]), fmt(jobs[[z[0] for z in jobs].index(tag)][4])],
                        "qprime_interval_on_declared": [fmt(qprimeI[0]), fmt(qprimeI[1])],
                        "bounds": {"value": fmt(B), "derivative": fmt(D)},
                        "state_count": len(states), "boundary_edge_count": len(boundary),
                        "boundary_reason_counts": {k: sum(e["reason"] == k for e in boundary) for k in ("value", "derivative")},
                        "zero_not_reachable": True, "boundary_certificates": boundary})

    expected_sizes = {"m10": 180, "m11": 2, "01": 12, "12": 2}
    assert {r["tag"]: r["state_count"] for r in results} == expected_sizes
    output = {"status": "passed_independent_fraction_cubic_enumeration_and_bfs",
              "method": "pure Python Fraction Sturm/root intervals/polynomial remainders; no Wolfram or pilots read",
              "candidate_box": [-5, 5], "candidate_count": len(candidates),
              "candidates": [{"c1": a, "c2": b, "real_roots": n} for a,b,n in candidates],
              "representative_bfs": results,
              "claim_boundary": "This verifies the stated degree-three common-factor obstruction; it does not prove a general S-integrality or degree>=4 result."}
    here = Path(__file__).resolve(); output["source_sha256"] = sha256(here.read_bytes()).hexdigest()
    here.with_suffix(".json").write_text(json.dumps(output, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"status": output["status"], "candidate_count": len(candidates),
                      "states": {r["tag"]: r["state_count"] for r in results},
                      "boundary_edges": {r["tag"]: r["boundary_edge_count"] for r in results},
                      "source_sha256": output["source_sha256"]}, sort_keys=True))


if __name__ == "__main__":
    main()
