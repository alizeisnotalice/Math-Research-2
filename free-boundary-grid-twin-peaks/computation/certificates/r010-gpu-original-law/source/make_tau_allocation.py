#!/usr/bin/env python3
"""Build a new-seed production allocation from tau-r001 pilot only.

This is a variance heuristic, not a confidence procedure.  Pilot observations
are deliberately not included in later production totals.
"""
import argparse
import json
from math import exp, log, sqrt
from pathlib import Path


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("pilot", type=Path)
    ap.add_argument("--n", type=int, default=7)
    ap.add_argument("--batches", type=int, default=512)
    ap.add_argument("--floor", type=int, default=4)
    args = ap.parse_args()
    rows = [json.loads(x) for x in args.pilot.read_text().splitlines()]
    n = args.n
    strata = [x for x in rows if x["kind"] == "stratum" and x["n"] == n]
    total = next(x for x in rows if x["kind"] == "summary" and x["n"] == n)
    assert len(strata) == n * n - 2
    assert args.batches >= args.floor * len(strata)
    q = total["Q_euclidean_hat"]
    scores = []
    for x in strata:
        if x["ESS"] is None:
            score = 0.0
        else:
            M = x["samples"]
            la = x["log_A2_tau_hat"]
            # ESS=(sum W)^2/sum W^2, and Ahat=sum W/M.
            ew2 = exp(2 * (la + log(M)) - log(x["ESS"]) - log(M))
            var = max(0.0, ew2 - exp(2 * la))
            # Delta-method proxy for W*(e2-NQ); +1 retains A2 precision.
            score = sqrt(var) * (1.0 + abs(n * n * (x["Q_tau_euclidean"] - q)))
        scores.append(score)
    extra = args.batches - args.floor * len(strata)
    if sum(scores) == 0:
        additions = [0] * len(scores)
    else:
        raw = [extra * s / sum(scores) for s in scores]
        additions = [int(x) for x in raw]
        for _, i in sorted((raw[i] - additions[i], i) for i in range(len(raw)))[::-1][:extra - sum(additions)]:
            additions[i] += 1
    allocation = [{"tau": x["tau"], "batches": args.floor + additions[i], "pilot_score": scores[i],
                   "pilot_alive": x["alive"], "pilot_ESS": x["ESS"]}
                  for i, x in enumerate(strata)]
    assert sum(x["batches"] for x in allocation) == args.batches
    print(json.dumps({"n": n, "source": str(args.pilot), "pilot_not_combined": True,
                      "seed_domain": "0xA022E*** disjoint from pilot 0xA022D***",
                      "floor_batches": args.floor, "total_batches": args.batches,
                      "score": "sqrt(Var_hat(W))*(1+|n^2(Q_tau-Q_total)|)",
                      "allocation": allocation}, indent=2))


if __name__ == "__main__":
    main()
