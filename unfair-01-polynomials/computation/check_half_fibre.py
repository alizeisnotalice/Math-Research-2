"""High-precision audit of the half-fibre obstruction for H(x^14)H(x^15).

This is evidence for a *local obstruction only*: K(x^7) selects seven of the
fourteen lifts over every root of H.  It verifies the numerical factorisation
K(z)K(-z)=H(z^2); the exact identities used in the note are independent of
these floating-point values.
"""
from __future__ import annotations

import itertools
import argparse
import json
from pathlib import Path

import mpmath as mp


HERE = Path(__file__).resolve().parent
parser=argparse.ArgumentParser();parser.add_argument('--out',type=Path,required=True);OUT=parser.parse_args().out;assert Path('/Volumes/data').is_mount()
mp.mp.dps = 100


def mul(p, q):
    out = [mp.mpc(0) for _ in range(len(p) + len(q) - 1)]
    for i, x in enumerate(p):
        for j, y in enumerate(q):
            out[i + j] += x * y
    return out


def coeffs_from_roots(roots):
    # Ascending coefficients of the monic product prod(x-r).
    p = [mp.mpc(1)]
    for r in roots:
        p = mul(p, [-r, 1])
    return p


def clean(z):
    return {"re": mp.nstr(mp.re(z), 30), "im": mp.nstr(mp.im(z), 8)}


def main():
    # H(x)=1+x+x^2+x^5+x^6, in descending order for polyroots.
    roots = mp.polyroots([1, 1, 0, 0, 1, 1, 1], maxsteps=1000, error=False)
    upper = sorted((z for z in roots if mp.im(z) > 0), key=lambda z: float(mp.arg(z)))
    assert len(upper) == 3
    h_ascending = [mp.mpc(1), 1, 1, 0, 0, 1, 1]
    target = [h_ascending[i // 2] if i % 2 == 0 else mp.mpc(0)
              for i in range(13)]
    records = []
    for signs in itertools.product((0, 1), repeat=3):
        chosen = []
        for z, bit in zip(upper, signs):
            theta = mp.sqrt(z)
            if bit:
                theta = -theta
            chosen += [theta, theta.conjugate()]
        k = coeffs_from_roots(chosen)
        km = [c * ((-1) ** i) for i, c in enumerate(k)]
        product = mul(k, km)
        residual = max(abs(product[i] - target[i])
                       for i in range(13))
        # K(x^7)'s selected roots are all seven roots over each theta.
        lifted = []
        zeta = mp.e ** (2j * mp.pi / 7)
        for theta in chosen:
            gamma = mp.e ** (mp.log(theta) / 7)
            lifted += [gamma * (zeta ** j) for j in range(7)]
        moments = []
        for power in range(1, 7):
            moments.append(max(abs(sum(r ** power for r in lifted)),
                               abs(sum(r ** (-power) for r in lifted))))
        records.append({
            "signs_on_three_conjugate_pairs": list(signs),
            "K_coefficients_ascending": [clean(c) for c in k],
            "K_constant_residual_from_1": mp.nstr(abs(k[0] - 1), 8),
            "K_real_coefficient_max_imag": mp.nstr(max(abs(mp.im(c)) for c in k), 8),
            "K_times_Kminus_equals_Hz2_max_residual": mp.nstr(residual, 8),
            "selected_42_root_moment_maxima_k_1_to_6": [mp.nstr(x, 8) for x in moments],
        })
    # Exact support test for C=H(x^14)H(x^15).
    support_h = (0, 1, 2, 5, 6)
    support = [14 * i + 15 * j for i in support_h for j in support_h]
    assert len(support) == len(set(support)) == 25
    # Exact coefficient map for K(-x^7) H(x^15): no cancellation can hide
    # an individual coefficient of K(-x), since (i,j) -> 7i+15j is injective.
    b_support = [7 * i + 15 * j for i in range(7) for j in support_h]
    assert len(b_support) == len(set(b_support)) == 35
    result = {
        "precision_decimal_digits": mp.mp.dps,
        "H_support": list(support_h),
        "a": 14,
        "b": 15,
        "C_support_size": len(support),
        "C_binary_exact": True,
        "B_coefficient_map_support_size": len(b_support),
        "B_coefficient_map_injective_exact": True,
        "construction": "For each alpha root, choose theta with theta^2=alpha; K has theta and conjugates. K(x^7) selects the seven roots x^7=theta, exactly half of each x^14=alpha fibre.",
        "balance": "m_a,i=7, m_b,i=0 for all i; b*m_a,i+a*m_b,i=15*7=105.",
        "records": records,
    }
    OUT.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps({
        "output": str(OUT), "records": len(records), "C_binary_exact": True,
        "B_coefficient_map_injective_exact": True,
        "max_factorisation_residual": max(float(r["K_times_Kminus_equals_Hz2_max_residual"])
                                          for r in records),
        "max_moment_residual": max(float(x) for r in records
                                   for x in r["selected_42_root_moment_maxima_k_1_to_6"]),
    }, indent=2))


if __name__ == "__main__":
    main()
