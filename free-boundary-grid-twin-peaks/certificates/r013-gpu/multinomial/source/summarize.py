#!/usr/bin/env python3
"""Finite-repeat five-moment and delta-SE summary; no confidence claim."""
from collections import defaultdict
import json, math, sys

rows = defaultdict(list)
for line in open(sys.argv[1]):
    x=json.loads(line)
    if x.get("kind")=="replicate": rows[x["tau"]].append(x)

def record(label, xs):
    r=len(xs); z=[x["Z"] for x in xs]; u=[x["U"] for x in xs]
    mz=sum(z)/r; mu=sum(u)/r; q=mu/mz if mz else None
    def avg(a): return sum(a)/r
    # Plug-in delta estimate from independent ideal repeats; descriptive only.
    d=[ui-q*zi for zi,ui in zip(z,u)] if q is not None else []
    delta_se=(math.sqrt(sum((x-sum(d)/r)**2 for x in d)/(r-1)/r)/mz if r>1 and mz else None)
    return {"tau":label,"reps":r,"Z_mean":mz,"U_mean":mu,"ratio_numerical_only":q,
            "Z2_mean":avg([x*x for x in z]),"ZU_mean":avg([a*b for a,b in zip(z,u)]),"U2_mean":avg([x*x for x in u]),
            "delta_SE_ideal_repeat_model_only":delta_se,"dead":sum(x["dead"] for x in xs)}

summary=[record(t,xs) for t,xs in sorted(rows.items())]
byrep=defaultdict(lambda:[0.0,0.0])
for xs in rows.values():
    for x in xs: byrep[x["replicate"]][0]+=x["Z"]; byrep[x["replicate"]][1]+=x["U"]
summary.append(record("all_tau_combined", [{"Z":z,"U":u,"dead":False} for _,(z,u) in sorted(byrep.items())]))
for x in summary: print(json.dumps(x,sort_keys=True))
