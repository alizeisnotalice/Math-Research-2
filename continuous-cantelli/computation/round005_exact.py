"""Exact three-time Gaussian-marginal obstruction. Not a deterministic graph."""
from fractions import Fraction as F
from collections import defaultdict
from math import factorial
from pathlib import Path
import json,sys,os,platform,hashlib
s,d,v,u,a=F(1,2),F(1,2),F(1,8),F(5,8),F(2)
# Entries are time, Gaussian variance, Gaussian mean, weight polynomial in epsilon.
terms=[(s,v,a,0,1),(F(1),F(1),F(0),1,0),(F(1),v,a,0,-1),
       (F(1),u,a,0,-1),(1+d,u,a,0,1)]
caloric=defaultdict(lambda:[F(0),F(0)]);spatial=defaultdict(lambda:[F(0),F(0)])
for t,var,mean,c0,c1 in terms:
 caloric[(mean,(var-t)/2)][0]+=c0;caloric[(mean,(var-t)/2)][1]+=c1
 spatial[(mean,var)][0]+=c0;spatial[(mean,var)][1]+=c1
caloric={k:z for k,z in caloric.items() if z!=[0,0]};spatial={k:z for k,z in spatial.items() if z!=[0,0]}
assert caloric=={(F(0),F(0)):[F(1),F(0)]}
assert spatial=={(F(0),F(1)):[F(1),F(0)]}
# epsilon=(3/5)*exp(-16/3); positivity bound is 73/80.
assert 8<3**2 and F(8,5)<F(4,3)**2
assert F(64,21)>3
exp3_lower=sum(F(3)**j/factorial(j) for j in range(5))
assert exp3_lower==F(131,8)>16
middle_removed_upper=F(3,5)*(F(3,16)+F(4,3))
assert middle_removed_upper==F(73,80)<1
# Early density ratio=(6/5)*exp(-3*(x-8/3)^2).
center=a*s/(s-v);curvature=(s-v)/(2*s*v)
assert center==F(8,3) and curvature==3
radius=F(1,5);ratio_lower=F(6,5)*(1-curvature*radius**2)
assert ratio_lower==F(132,125)>1
# Smoothing restores positivity at t=75/128 before the middle source arrives.
restore=F(3,8)/(1-F(3,5)**2)
assert restore==F(75,128)
mean_time=[sum(F(c0)*t for t,var,mean,c0,c1 in terms),sum(F(c1)*t for t,var,mean,c0,c1 in terms)]
second_time=[sum(F(c0)*t*t for t,var,mean,c0,c1 in terms),sum(F(c1)*t*t for t,var,mean,c0,c1 in terms)]
assert mean_time==[1,0] and second_time==[1,F(1,2)]
# Verify the Gaussian-ratio completed square as exact rational polynomial coefficients.
ratio_checks=[]
for rr,vv in [(F(1,2),F(1,8)),(F(1),F(1,8)),(F(1),F(5,8)),(F(17,32),F(5,32))]:
 K=(rr-vv)/(2*rr*vv);xx=a*rr/(rr-vv);cc=a*a/(2*(rr-vv))
 lhs=[F(1,2)/rr-F(1,2)/vv,a/vv,-a*a/(2*vv)]
 rhs=[-K,2*K*xx,cc-K*xx*xx]
 assert lhs==rhs
 ratio_checks.append([str(rr),str(vv)])
cert={'epsilon':'(3/5)*exp(-16/3)','times':['1/2','1','3/2'],
 'spatial_marginal':'N(0,1) exactly','full_caloric_transform':'1 exactly',
 'middle_density_lower_ratio':'7/80','negative_time':'1/2',
 'negative_interval':[str(center-radius),str(center+radius)],
 'early_source_ratio_lower':str(ratio_lower),'negative_Q_upper_ratio':'-7/125',
 'positivity_restored_time':str(restore),'E_time':'1','Var_time':'epsilon/2',
 'exp3_lower':str(exp3_lower),'middle_removed_upper':str(middle_removed_upper),
 'gaussian_ratio_identity_checks':ratio_checks,'scope':'positive three-time measure with Gaussian spatial marginal; NOT a deterministic graph'}
out={'certificate':cert,'runtime':{'python':platform.python_version(),'pid':os.getpid()},
 'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
Path(sys.argv[1]).write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
