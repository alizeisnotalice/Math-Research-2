#!/usr/bin/env python3
"""Make the fresh (not-resume) cap-6 total-degree 47/48 input ledger."""
from __future__ import annotations
import argparse, hashlib, itertools, json, math, os
from pathlib import Path

HERE = Path(__file__).resolve().parent
FROZEN_SHA = "d68f1010c61ad1d911749742fef8c790c0a8064bf34512b402e7e65e4703e710"
TOTALS = (47, 48)
CAP = 6

def sha(p: Path) -> str: return hashlib.sha256(p.read_bytes()).hexdigest()

def independent_supports(degree: int, cap: int):
    for size in range(2, min(degree + 1, cap) + 1):
        for x in itertools.combinations(range(1, degree), size - 2):
            yield (0, *x, degree)

def support_count(degree: int, cap: int) -> int:
    return sum(math.comb(degree - 1, k) for k in range(0, min(cap - 2, degree - 1) + 1))

def expected_jobs():
    jobs=[]; coverage=[]
    for total in TOTALS:
        for m in range(1, total // 2 + 1):
            n=total-m; supports=list(independent_supports(m,CAP))
            if len(supports) != support_count(m,CAP): raise RuntimeError("independent enumeration count mismatch")
            coverage.append({"total":total,"m":m,"n":n,"expected_I_count":len(supports),"binomial_count":support_count(m,CAP)})
            jobs.extend({"total":total,"m":m,"n":n,"index":i,"I":list(I)} for i,I in enumerate(supports))
    return jobs,coverage

def validate(plan: dict):
    if plan.get("format") != "a30-total47-48-cap6-v1": raise RuntimeError("plan format")
    if plan.get("parameters") != {"totals":[47,48],"cap":6,"remote_workers":15,"remote_seconds_limit":14400}: raise RuntimeError("plan parameter contract")
    if plan["frozen_support_sha256"] != FROZEN_SHA: raise RuntimeError("frozen binding")
    wanted,coverage=expected_jobs()
    if plan.get("jobs") != wanted or plan.get("coverage") != coverage: raise RuntimeError("jobs are not the independent full enumeration")
    keys=[(x["total"],x["m"],x["n"],x["index"],tuple(x["I"])) for x in plan["jobs"]]
    if len(keys)!=len(set(keys)): raise RuntimeError("duplicate job")
    # Construct each degree's independent list once: per-row reconstruction
    # would turn this validation into a quadratic accidental workload.
    support_table={m:list(independent_supports(m,CAP)) for m in range(1,25)}
    for x in plan["jobs"]:
        if x["total"] not in TOTALS or x["m"]+x["n"] != x["total"] or not(1<=x["m"]<=x["total"]//2): raise RuntimeError("out-of-range degree pair")
        choices=support_table[x["m"]]
        if not(0<=x["index"]<len(choices)) or list(choices[x["index"]]) != x["I"]: raise RuntimeError("wrong I/index")

def main():
    ap=argparse.ArgumentParser(); ap.add_argument("--out",type=Path,required=True); a=ap.parse_args()
    frozen=HERE/"plan/frozen_support.py"
    if sha(frozen)!=FROZEN_SHA: raise SystemExit("frozen classifier hash mismatch")
    if not Path("/Volumes/data").is_mount() or a.out.exists(): raise SystemExit("data unavailable or output already exists")
    jobs,coverage=expected_jobs()
    plan={"format":"a30-total47-48-cap6-v1","scope":"Fresh exact singleton-unit support classification over total degrees 47 and 48 only; no total46 resume and no general theorem.","parameters":{"totals":[47,48],"cap":6,"remote_workers":15,"remote_seconds_limit":14400},"frozen_support_sha256":FROZEN_SHA,"coverage":coverage,"jobs":jobs,"job_total":len(jobs)}
    validate(plan); a.out.mkdir(parents=True); target=a.out/"support48-plan.json"; target.write_text(json.dumps(plan,indent=2,sort_keys=True)+"\n")
    print(json.dumps({"job_total":len(jobs),"coverage":coverage,"plan_sha256":sha(target)},sort_keys=True))
if __name__=="__main__": main()
