#!/Users/zzh/.venvs/main/bin/python
"""Checkpointable exact necessary screens for a selected b-layer of core6.

The input ledger is conditional on the separately audited root-envelope bounds.
This program applies only further necessary tests: positive values of both real
embeddings at seven rational x>0, and exact integer Schur tests for every root
of H=P*Pbar lying in 1/2<|z|<2.  It never treats a surviving row as a model.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import sys
from pathlib import Path


HERE = Path(__file__).resolve().parent
POINTS = ((1, 4), (1, 2), (3, 4), (1, 1), (4, 3), (2, 1), (4, 1))


def primitive(coeffs: list[int]) -> list[int]:
    g = 0
    for value in coeffs:
        g = math.gcd(g, abs(value))
    if g > 1:
        coeffs = [value // g for value in coeffs]
    if coeffs[-1] < 0:
        coeffs = [-value for value in coeffs]
    return coeffs


def schur_witness(ascending: list[int]) -> dict | None:
    """None exactly when all roots lie strictly inside the unit disk."""
    f = primitive(ascending)
    while len(f) > 1:
        constant, leading = f[0], f[-1]
        if abs(constant) >= abs(leading):
            return {"degree": len(f) - 1, "constant": constant, "leading": leading}
        nxt = [leading * f[i] - constant * f[-1-i] for i in range(len(f))]
        assert nxt[0] == 0
        f = primitive(nxt[1:])
    return None


def scaled(coeffs: list[int]) -> list[int]:
    return [coefficient * (2**i) for i, coefficient in enumerate(coeffs)]


def h_coefficients(b: int, T: int, N: int, ell: int, d: int, e: int) -> list[int]:
    """Ascending coefficients of H=P*Pbar, verified symbolically in self_test."""
    return [
        1,
        2*b,
        b*b + T,
        2*b*T + 2*ell,
        2*d + b*b*T + 2*b*ell + N,
        2*e + 2*b*d + 2*b*N + ell*T,
        2 + 2*b*e + d*T + b*b*N + b*ell*T + ell*ell,
        2*b + e*T + b*d*T + 2*d*ell,
        T + e*b*T + 2*e*ell + d*d,
        b*T + 2*ell + 2*d*e,
        2*d + e*e,
        2*e,
        1,
    ]


def positive_axis_witness(b: int, T: int, delta: int, ell: int, d: int, e: int) -> dict | None:
    """Exact test P(x)>0 and Pbar(x)>0 at every fixed rational point."""
    for p, q in POINTS:
        # q^6 P(p/q)=R+aK, K=p^2q^3(q+bp).
        K = p*p*q**3 * (q + b*p)
        R = q**6 + b*p*q**5 + ell*p**3*q**3 + d*p**4*q**2 + e*p**5*q + p**6
        L = 2*R + T*K
        discriminant_numerator = L*L - delta*K*K
        if L <= 0 or discriminant_numerator <= 0:
            return {"point": [p, q], "L": L, "product_numerator": discriminant_numerator}
    return None


def annulus_witness(h: list[int]) -> dict | None:
    outer = schur_witness(scaled(h))
    if outer is not None:
        return {"orientation": "H(2x)", **outer}
    inner = schur_witness(scaled(list(reversed(h))))
    if inner is not None:
        return {"orientation": "reverse(H)(2x)", **inner}
    return None


def self_test() -> None:
    assert schur_witness([-1, 2]) is None       # root 1/2
    assert schur_witness([-2, 1]) is not None   # root 2
    assert schur_witness([-1, 1]) is not None   # boundary root 1
    # Hand convolution of P=1+x+x^2+(x+2)x^3+3x^4+4x^5+x^6,
    # with T=1,N=0,b=1,ell=2; SymPy independently checks this vector below.
    assert h_coefficients(1, 1, 0, 2, 3, 4) == [
        1, 2, 2, 6, 11, 16, 19, 21, 30, 29, 22, 8, 1
    ]


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--b", type=int, action="append", required=True,
                        help="one or more nonzero b layers; no all-layers default")
    args = parser.parse_args()
    if any(abs(b) > 7 or b == 0 for b in args.b):
        raise SystemExit("each --b must lie in [-7,-1] union [1,7]")
    self_test()
    ledger_path = HERE / "core6_linked_conditional_ledger_corrected.json"
    ledger = json.loads(ledger_path.read_text())
    requested = set(args.b)
    mothers = [row for row in ledger["rows"] if row["b"] in requested]
    counts = {"mothers": len(mothers), "raw_d_e": 0, "positive_axis": 0, "annulus_schur": 0}
    axis_first = {f"{p}/{q}": 0 for p, q in POINTS}
    schur_first = {"H(2x)": 0, "reverse(H)(2x)": 0}
    survivors: list[dict[str, int]] = []
    examples: dict[str, object] = {}
    # Compact, order-sensitive audit trail for every candidate's exact first
    # failure.  This avoids creating a multi-gigabyte JSONL for 22.9m rows.
    failure_digest = hashlib.sha256()
    for row in mothers:
        b, T, N, delta, ell = (row[k] for k in ("b", "T", "N", "delta", "ell"))
        for d in range(-21, 22):
            for e in range(-7, 8):
                counts["raw_d_e"] += 1
                axis = positive_axis_witness(b, T, delta, ell, d, e)
                if axis is not None:
                    axis_first[f"{axis['point'][0]}/{axis['point'][1]}"] += 1
                    examples.setdefault("positive_axis", {"row": [b,T,N,delta,ell,d,e], "witness": axis})
                    failure_digest.update(repr((b,T,N,delta,ell,d,e,"axis",axis)).encode() + b"\n")
                    continue
                counts["positive_axis"] += 1
                h = h_coefficients(b, T, N, ell, d, e)
                schur = annulus_witness(h)
                if schur is not None:
                    schur_first[schur["orientation"]] += 1
                    examples.setdefault("annulus_schur", {"row": [b,T,N,delta,ell,d,e], "witness": schur})
                    failure_digest.update(repr((b,T,N,delta,ell,d,e,"schur",schur)).encode() + b"\n")
                    continue
                counts["annulus_schur"] += 1
                survivors.append({**row, "d": d, "e": e})
                failure_digest.update(repr((b,T,N,delta,ell,d,e,"survivor")).encode() + b"\n")
    assert ((counts["raw_d_e"] - counts["positive_axis"])
            + (counts["positive_axis"] - counts["annulus_schur"])
            + counts["annulus_schur"] == counts["raw_d_e"])
    result = {
        "status": "NECESSARY_SCREENS_ONLY",
        "b_layers": sorted(requested),
        "input_ledger_sha256": hashlib.sha256(ledger_path.read_bytes()).hexdigest(),
        "counts": counts,
        "axis_first_failure": axis_first,
        "schur_first_failure": schur_first,
        "first_failure_record_schema": "(b,T,N,delta,ell,d,e,kind,witness) in deterministic ledger/d/e order",
        "first_failure_stream_sha256": failure_digest.hexdigest(),
        "examples": examples,
        "survivors": survivors,
        "script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "scope": "A survivor is not claimed to define G,Q,A,B or a counterexample.",
    }
    suffix = "_".join(str(b) for b in sorted(requested)).replace("-", "neg")
    output = HERE / f"core6_de_screen_b_{suffix}.json"
    output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps({k: result[k] for k in ("status", "b_layers", "counts", "axis_first_failure", "schur_first_failure")}, indent=2))
    print("output_sha256=" + hashlib.sha256(output.read_bytes()).hexdigest())


if __name__ == "__main__":
    main()
