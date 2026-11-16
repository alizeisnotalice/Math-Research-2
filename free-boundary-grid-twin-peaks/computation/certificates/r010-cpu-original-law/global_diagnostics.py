#!/usr/bin/env python3
"""Global weight diagnostics across forced-second-birth-rank strata.

Each stratum estimate is S_tau/N_tau.  With raw trajectories considered as
importance contributions W/N_tau, the diagnostic global ESS is
    (sum_tau S_tau/N_tau)^2 / sum_tau S2_tau/N_tau^2.
The JSON summaries store S_tau/N_tau and S_tau^2/S2_tau, so this can be
recomputed without the executable.  This is descriptive, not an error bound.
"""
import json
import sys
from pathlib import Path

out = Path(sys.argv[1])
rows = [json.loads(x) for x in (out / "summary.jsonl").read_text().splitlines()]
parts=[]
for r in rows:
    mean=float(r["A2_mean"]); ess=float(r["A2_ess"]); top=float(r["A2_top_share"])
    sq = mean*mean/ess if ess else 0.0
    parts.append({"tau":r["tau"], "samples":r["samples"], "dead":r["dead"],
                  "mean_A2":mean, "weight_square_scaled":sq,
                  "within_ESS":ess, "within_top_share":top})
total=sum(x["mean_A2"] for x in parts)
den=sum(x["weight_square_scaled"] for x in parts)
for x in parts:
    x["mass_share"]=x["mean_A2"]/total if total else 0.0
    x["global_max_contribution_bound"]=x["within_top_share"]*x["mass_share"]
ordered=sorted(parts,key=lambda x:x["mass_share"],reverse=True)
cum=0.0; major=[]
for x in ordered:
    cum+=x["mass_share"];major.append(x["tau"])
    if cum>=.90:break
report={"n":rows[0]["n"],"total_A2_estimate":total,
        "global_weight_ESS":total*total/den if den else 0.0,
        "global_max_normalized_contribution":max((x["global_max_contribution_bound"] for x in parts),default=0.0),
        "major_mass_tau_to_90pct":major,"major_mass_cumulative_share":cum,
        "strata":parts,
        "missing_for_Q_uncertainty":"stored output has no sum (W*euclidean2)^2, so no ESS or standard error for the Q numerator can be reconstructed",
        "interpretation":"descriptive finite-sample weight diagnostics; not a confidence interval or conditional-law certificate"}
(out/"global-diagnostics.json").write_text(json.dumps(report,indent=2)+"\n")
print(json.dumps({k:v for k,v in report.items() if k!="strata"},indent=2))
