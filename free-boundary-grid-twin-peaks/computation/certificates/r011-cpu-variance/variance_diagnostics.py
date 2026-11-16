#!/usr/bin/env python3
"""Aggregate raw forced-tau SIS moments; all uncertainty is descriptive only."""
import json
import math
import sys
from collections import defaultdict
from pathlib import Path

def num(row, key):
    return float(row[key])

def rows_at(path, name):
    return [json.loads(x) for x in (path / name).read_text().splitlines()
            if json.loads(x).get("kind") != "metadata"]

def calculate(rows):
    by = defaultdict(lambda: defaultdict(float))
    for r in rows:
        d = by[r["tau"]]
        d["N"] += r["samples"]
        for k in ("sumW", "sumW2", "sumWf", "sumW2f", "sumW2f2"):
            d[k] += num(r, k)
        # sum of batch maxima is not valid; retain the true maximum instead.
        d["maxW"] = max(d["maxW"], num(r, "A2_top_share") * num(r, "sumW"))
        d["dead"] += r["dead"]
    A = sum(d["sumW"] / d["N"] for d in by.values())
    B = sum(d["sumWf"] / d["N"] for d in by.values())
    q = B / A if A else float("nan")
    second = sum(d["sumW2"] / d["N"]**2 for d in by.values())
    ess = A*A / second if second else 0.0
    maxcontrib = max((d["maxW"] / d["N"] / A for d in by.values()), default=0.0)
    delta_var_num = 0.0
    strata = []
    for tau, d in sorted(by.items()):
        N = d["N"]
        md = (d["sumWf"] - q*d["sumW"]) / N
        md2 = (d["sumW2f2"] - 2*q*d["sumW2f"] + q*q*d["sumW2"]) / N
        # Unbiased sample variance of D=W(f-q), followed by Var(mean D).
        rawvar = md2-md*md
        tolerance = 1e-12 * max(1.0, abs(md2), abs(md*md))
        if rawvar < -tolerance:
            raise ValueError(f"substantive negative variance in tau={tau}: {rawvar}")
        var_d = max(0.0, rawvar * N/(N-1)) if N > 1 else 0.0
        delta_var_num += var_d / N
        strata.append({"tau": tau, "samples": int(N),
                        "mass_share": (d["sumW"]/N/A if A else 0.0),
                        "ess": ((d["sumW"]/N)**2/(d["sumW2"]/N**2)
                                if d["sumW2"] else 0.0),
                        "top_contribution_share": (d["maxW"]/N/A if A else 0.0),
                        "zero_contribution": int(d["dead"]),
                        "dead": int(d["dead"]), "delta_var_component": var_d/N/(A*A) if A else 0.0})
    strata.sort(key=lambda x: x["mass_share"], reverse=True)
    cum = 0.0; mass90=[]
    for x in strata:
        cum += x["mass_share"]; mass90.append(x["tau"])
        if cum >= .90: break
    return {"strata": len(by), "A2_estimate": A, "Q_e2": q,
            "global_ess": ess, "global_max_contribution_share": maxcontrib,
            "delta_method_se_Q": math.sqrt(delta_var_num)/A if A else None,
            "dead_total": int(sum(d["dead"] for d in by.values())),
            "mass90_taus": sorted(mass90), "strata_by_mass": strata}

def main():
    root = Path(sys.argv[1])
    if len(sys.argv) < 3:
        raise SystemExit("usage: variance_diagnostics.py ROOT COMPLETED_DOMAIN [COMPLETED_DOMAIN ...]")
    dirs = [root / x for x in sys.argv[2:]]
    report = {"description": "Raw-moment aggregation. ESS, maximum contribution, delta SE, and across-domain SD are finite-sample diagnostics, not certified error bars.", "domains": {}}
    qs_by_n=defaultdict(list)
    raw_by_n=defaultdict(list)
    for p in dirs:
        if not (p / "batches.jsonl").exists():
            raise ValueError(f"missing batches.jsonl: {p}")
        rs = rows_at(p, "batches.jsonl")
        meta = json.loads((p / "batches.jsonl").read_text().splitlines()[0])
        expected = set(range(2, meta["n"]*meta["n"]+1))
        got = {r["tau"] for r in rs}
        if got != expected:
            raise ValueError(f"incomplete tau set in {p.name}: expected {len(expected)}, got {len(got)}")
        if len({r["seed"] for r in rs}) != len(rs):
            raise ValueError(f"repeated seed in {p.name}")
        for tau in expected:
            rr=[r for r in rs if r["tau"]==tau]
            if len(rr) != meta["batches"] or any(r["samples"] != meta["samples_per_tau_batch"] for r in rr):
                raise ValueError(f"incomplete batch/sample count in {p.name}, tau={tau}")
        ans = calculate(rs)
        ans.update({"n":meta["n"], "seed_domain":meta["seed_domain"], "directory":p.name})
        report["domains"][p.name]=ans; qs_by_n[ans["n"]].append(ans["Q_e2"]); raw_by_n[ans["n"]].extend(rs)
        (p / "variance-diagnostics.json").write_text(json.dumps(ans, indent=2)+"\n")
    report["across_domains_by_n"]={}
    for n, qs in qs_by_n.items():
        mean=sum(qs)/len(qs)
        sd=math.sqrt(sum((x-mean)**2 for x in qs)/(len(qs)-1)) if len(qs)>1 else None
        report["across_domains_by_n"][str(n)]={"count":len(qs),"mean_Q":mean,
          "sample_sd_Q":sd,
          "descriptive_se_of_domain_mean_Q":(sd/math.sqrt(len(qs)) if sd is not None else None),
          "pooled_raw_moment_Q":calculate(raw_by_n[n])["Q_e2"]}
    (root / "all-variance-diagnostics.json").write_text(json.dumps(report, indent=2)+"\n")
    print(json.dumps(report, indent=2))
if __name__ == "__main__": main()
