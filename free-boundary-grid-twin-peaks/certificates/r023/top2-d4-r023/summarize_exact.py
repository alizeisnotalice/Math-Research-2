#!/usr/bin/env python3
"""Create exact rational tau2/K=2 comparisons from n3,...,n7 JSON runs."""
import argparse
import json
from fractions import Fraction
from pathlib import Path


def rational(num, den):
    value = Fraction(int(num), int(den))
    return {"numerator": str(value.numerator), "denominator": str(value.denominator),
            "decimal": format(value.numerator / value.denominator, ".16g")}


def one(path):
    x = json.loads(Path(path).read_text())
    B, A2 = int(x["Btop2"]), int(x["A2"])
    td, tf, tg = (int(x[k]) for k in ("tau2_squared_euclidean_total", "tau2_firstsum", "tau2_secondsum"))
    kd, kf, kg = (int(x[k]) for k in ("squared_euclidean_total", "firstsum", "secondsum"))
    hist = {key: value for key, value in x["tau2_distance_counts"].items()}
    assert sum(map(int, hist.values())) == B
    assert tf - 2 * tg == 4 * td and kf - 2 * kg == 4 * kd
    n2 = x["n"] ** 2
    tau_q, tau_v, tau_c = Fraction(td, n2 * B), Fraction(tf, 8 * B), Fraction(tg, 4 * B)
    k2_q, k2_v, k2_c = Fraction(kd, n2 * A2), Fraction(kf, 8 * A2), Fraction(kg, 4 * A2)
    assert tau_q == 2 * (tau_v - tau_c) / n2
    assert k2_q == 2 * (k2_v - k2_c) / n2
    return {"n": x["n"], "Btop2": str(B), "tau2_distance_counts": hist,
            "E_tau2_L_equals_A2_over_Btop2": rational(A2, B),
            "tau2_E_distance_squared": rational(td, B), "tau2_Q": rational(td, n2 * B),
            "tau2_V": rational(tf, 8 * B), "tau2_C": rational(tg, 4 * B),
            "K2_A2": str(A2), "K2_E_distance_squared": rational(kd, A2),
            "K2_Q": rational(kd, n2 * A2), "K2_V": rational(kf, 8 * A2),
            "K2_C": rational(kg, 4 * A2)}


parser = argparse.ArgumentParser()
parser.add_argument("inputs", nargs="+", help="n3..n7 result JSON files")
parser.add_argument("--out", required=True)
args = parser.parse_args()
out = {"kind": "exact-tau2-top2-D4-summary", "results": [one(path) for path in args.inputs],
       "scope_note": "tau=2 is a birth-rank law; K2 is the law conditioned on exactly two peaks. They are compared, not identified."}
Path(args.out).write_text(json.dumps(out, indent=2) + "\n")
print(json.dumps({"passed": True, "out": args.out, "n": [x["n"] for x in out["results"]]}))
