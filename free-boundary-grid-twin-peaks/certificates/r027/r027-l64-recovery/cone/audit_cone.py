#!/usr/bin/env python3
from fractions import Fraction
from pathlib import Path
import hashlib,json,sys
path=Path(sys.argv[1])
x=json.loads(path.read_text());n=x["n"];assert n==24 and x["R"]==32
h=list(map(Fraction,x["distance_L_weighted_probabilities"]))
assert len(h)==2*n-1 and h[0]==h[1]==0 and all(v>=0 for v in h)
s=sum(h);h=[v/s for v in h]
b=Fraction(x["tau2_normalizer_mean"])/Fraction(x["lift_full_normalizer_mean"])
rhs=sum(Fraction(d,2*(d-1))*h[d] for d in range(2,len(h)))
slack=rhs-b
r=12;tail=sum(h[r+1:]);tail_bound=Fraction(2*r,r-1)*(1-b)
out={"scope":"Exact rational necessary-condition check of the stored finite-decimal plug-in pair; positive slack is not an accuracy certificate", "summary_sha256":hashlib.sha256(path.read_bytes()).hexdigest(),"n":n,"R":32,"B_over_A_hat":str(b),"B_over_A_hat_decimal":float(b),"cone_slack":str(slack),"cone_slack_decimal":float(slack),"violates_necessary_cone":slack<0,"strict_tail_d_gt_12_decimal":float(tail),"necessary_tail_upper_decimal":float(tail_bound),"lift_ESS":x["liftZ_weight"]["ESS"],"lift_maxshare":x["liftZ_weight"]["maxshare"]}
print(json.dumps(out,indent=2))
