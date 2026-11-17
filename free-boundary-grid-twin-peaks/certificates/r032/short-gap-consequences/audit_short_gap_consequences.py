#!/usr/bin/env python3
"""Independent integer consequences of the n6 global and n7 centre audits."""
import argparse
from fractions import Fraction
import hashlib
import json
from pathlib import Path


def fraction_record(x):
    return {"numerator": x.numerator, "denominator": x.denominator,
            "decimal_display_only": float(x)}


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--n6", type=Path, required=True)
    p.add_argument("--n7-center", type=Path, required=True)
    p.add_argument("--output", type=Path, required=True)
    a = p.parse_args()
    n6 = json.loads(a.n6.read_text())
    n7 = json.loads(a.n7_center.read_text())
    results = {}
    for label, hs, b, a2 in (
        ("n6_global", n6["aggregate_h"], n6["reference_B"], n6["reference_A2"]),
        ("n7_center_only", n7["h"], n7["reference"]["Bv"], n7["reference"]["A2v"]),
    ):
        h = [int(hs[f"h{i}"]) for i in range(1, 6)]
        b, a2 = int(b), int(a2)
        residual = h[2] - 3*h[3] + 2*h[4]
        tails = {}
        for t in (4, 5, 6):
            count = b - sum(h[:t-1])
            assert count >= 0
            tails[str(t)] = {
                "count": count, "fraction_of_B": fraction_record(Fraction(count, b)),
                "QE8_slack": 8*b - t*(t+1)*count,
            }
        unseen_b = b-sum(h)
        unseen_a2 = a2-sum(t*v for t, v in enumerate(h, 1))
        assert unseen_b > 0 and unseen_a2 >= 6*unseen_b
        results[label] = {
            "h3_minus_3h4_plus_2h5": residual,
            "tails": tails,
            "L_at_least_6_B_fraction": fraction_record(Fraction(unseen_b, b)),
            "L_at_least_6_lifted_fraction": fraction_record(Fraction(unseen_a2, a2)),
            "L_at_least_6_conditional_mean": fraction_record(Fraction(unseen_a2, unseen_b)),
        }
    assert results["n6_global"]["h3_minus_3h4_plus_2h5"] == 0
    assert results["n7_center_only"]["h3_minus_3h4_plus_2h5"] == -13869959586449623922064347091056397498917889736704
    assert all(x["QE8_slack"] >= 0 for x in results["n6_global"]["tails"].values())
    assert all(x["QE8_slack"] < 0 for x in results["n7_center_only"]["tails"].values())
    out = {
        "status": "passed",
        "scope": "Finite n6 global and n7 fixed-highest-root centre consequences. No all-layer or asymptotic assertion; centre failure does not refute the global QE8 candidate.",
        "inputs": {"n6_sha256": hashlib.sha256(a.n6.read_bytes()).hexdigest(),
                   "n7_center_sha256": hashlib.sha256(a.n7_center.read_bytes()).hexdigest()},
        "results": results,
    }
    a.output.parent.mkdir(parents=True, exist_ok=True)
    a.output.write_text(json.dumps(out, indent=2) + "\n")
    print(json.dumps(out, indent=2))


if __name__ == "__main__":
    main()
