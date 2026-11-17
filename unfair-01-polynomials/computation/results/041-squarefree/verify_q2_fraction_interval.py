#!/usr/bin/env python3
"""Independent Fraction certificate for the Q^2 Newman-state obstruction.

Does not import or execute the Wolfram verifier, pilot, or its reported edges.
All polynomial and interval operations below are implemented with Python ints
and fractions.Fraction.
"""
from fractions import Fraction
from hashlib import sha256
from pathlib import Path
import json


def trim(a):
    a = list(a)
    while len(a) > 1 and a[-1] == 0:
        a.pop()
    return a


def rem_monic(a, m):
    a = list(a)
    while len(a) >= len(m):
        c = a[-1]
        if c:
            d = len(a) - len(m)
            for i, u in enumerate(m):
                a[d + i] -= c * u
        a.pop()
    return trim(a or [0])


def add(I, J):
    return (I[0] + J[0], I[1] + J[1])


def mul(I, J):
    z = (I[0] * J[0], I[0] * J[1], I[1] * J[0], I[1] * J[1])
    return min(z), max(z)


def poly_interval(a, I):
    out = (Fraction(a[-1]), Fraction(a[-1]))
    for c in reversed(a[:-1]):
        out = add(mul(out, I), (Fraction(c), Fraction(c)))
    return out


def qval(t):
    return 1 + t * t + t * t * t


def derivative(a):
    return [i * a[i] for i in range(1, len(a))] or [0]


def fmt(t):
    return f"{t.numerator}/{t.denominator}" if t.denominator != 1 else str(t.numerator)


def main():
    # Ascending coefficient convention. q=x^3+x^2+1 and q^2 are monic.
    q = [1, 0, 1, 1]
    q2 = [1, 0, 2, 2, 1, 2, 1]
    states = [
        (1, 0, 0, 0, 0, 0), (0, 1, 0, 0, 0, 0), (1, 1, 0, 0, 0, 0),
        (0, 0, 1, 0, 0, 0), (1, 0, 1, 0, 0, 0), (0, 1, 1, 0, 0, 0),
        (1, 1, 1, 0, 0, 0), (0, 0, 0, 1, 0, 0), (1, 0, 0, 1, 0, 0),
        (0, 0, 1, 1, 0, 0), (1, 0, 1, 1, 0, 0), (0, 1, 1, 1, 0, 0),
        (1, 1, 1, 1, 0, 0), (0, 0, 0, 1, 1, 0), (1, 0, 0, 1, 1, 0),
        (0, 1, 0, 1, 1, 0), (1, 1, 0, 1, 1, 0),
    ]
    state_set = set(states)

    def step(a, bit):
        t = a[5]
        return (bit - t, a[0], a[1] - 2*t, a[2] - 2*t, a[3] - t, a[4] - 2*t)

    # Coarse isolation and exact bisection. q'>0 on the coarse interval:
    # y<0 and 2+3y<0 there, hence their product is positive.
    lo, hi = Fraction(-22, 15), Fraction(-10, 7)
    assert qval(lo) < 0 < qval(hi)
    qpI = mul((lo, hi), (2 + 3*lo, 2 + 3*hi))
    assert qpI[0] > 0
    for _ in range(96):
        mid = (lo + hi) / 2
        if qval(mid) < 0:
            lo = mid
        else:
            hi = mid
    I = (lo, hi)
    assert qval(lo) < 0 < qval(hi)

    # Exact reduction checks for all 34 edges, independently of supplied list.
    edges = []
    for a in states:
        for bit in (0, 1):
            target = step(a, bit)
            raw = [bit] + list(a)
            assert rem_monic(raw, q2) == trim(target)
            edges.append((a, bit, target))
    outside = [edge for edge in edges if edge[2] not in state_set]
    assert len(edges) == 34 and len(outside) == 18
    assert (0, 0, 0, 0, 0, 0) not in state_set

    certificates = []
    for a, bit, z in outside:
        v = rem_monic(z, q)
        dv = rem_monic(derivative(list(z)), q)
        VI, DI = poly_interval(v, I), poly_interval(dv, I)
        if VI[0] > 3 or VI[1] < -3:
            kind = "value"
        elif DI[0] > 6 or DI[1] < -6:
            kind = "derivative"
        else:
            raise AssertionError((a, bit, z, v, dv, VI, DI))
        certificates.append({
            "from": list(a), "bit": bit, "target": list(z), "v_mod_q": v,
            "dv_mod_q": dv, "v_interval": [fmt(VI[0]), fmt(VI[1])],
            "dv_interval": [fmt(DI[0]), fmt(DI[1])], "separating_bound": kind,
        })

    # Reverse-tail bounds from |rho|>10/7: exact strict comparison.
    R = Fraction(10, 7)
    assert Fraction(1, 1) / (R - 1) < 3
    assert Fraction(1, 1) / (R - 1)**2 < 6
    out = {
        "status": "passed_independent_fraction_interval_q2_audit",
        "method": "pure Python integer polynomial arithmetic plus Fraction Horner intervals",
        "root_interval_after_96_bisections": [fmt(lo), fmt(hi)],
        "coarse_root_interval": ["-22/15", "-10/7"],
        "q_prime_interval_on_coarse_box": [fmt(qpI[0]), fmt(qpI[1])],
        "state_count": len(states), "edge_count": len(edges), "outside_edge_count": len(outside),
        "tail_bounds": {"value_lt": "3", "derivative_lt": "6"},
        "outside_certificates": certificates,
        "scope": "Independently verifies the finite state/edge obstruction used for Q^2; it does not itself prove unrelated crossed-model hypotheses.",
    }
    here = Path(__file__).resolve()
    out["source_sha256"] = sha256(here.read_bytes()).hexdigest()
    rec = here.with_suffix(".json")
    rec.write_text(json.dumps(out, indent=2, sort_keys=True) + "\n")
    print(json.dumps({k: out[k] for k in ("status", "state_count", "edge_count", "outside_edge_count", "source_sha256")}, sort_keys=True))


if __name__ == "__main__":
    main()
