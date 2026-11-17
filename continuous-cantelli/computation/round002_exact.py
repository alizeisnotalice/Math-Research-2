"""Exact rational certificate for a non-stopping caloric representing measure.
This measure is NOT a Cantelli counterexample: no graph and non-Gaussian spatial marginal.
No third-party packages, random numbers, or numerical integration.
"""
from fractions import Fraction as F
from collections import defaultdict
import json, platform, sys, os, hashlib
from pathlib import Path
s,v,w=F(1,2),F(1,8),F(2,3)
u=1-s+v
terms=[(s,v,w),(F(1),F(1),F(1)),(F(1),u,-w)]
exponents=defaultdict(F)
terminal=defaultdict(F)
for time,variance,weight in terms:
    exponents[(variance-time)/2]+=weight
    terminal[variance+1-time]+=weight
exponents={k:z for k,z in exponents.items() if z}
terminal={k:z for k,z in terminal.items() if z}
assert exponents=={F(0):F(1)}
assert terminal=={F(1):F(1)}
# Late measure p_1-w*p_u is positive: its subtracted ratio is maximized at 0.
late_ratio_squared=w*w/u
assert 0<u<1 and late_ratio_squared<1
# Early remainder p_s-w*p_v is negative around 0.
early_ratio_squared=w*w*s/v
assert early_ratio_squared==F(16,9)
# On |x|<=1/4, (w*p_v)/p_s = (4/3)*exp(-3*x*x)>=13/12.
early_ratio_lower=F(4,3)*(1-F(3,16))
assert early_ratio_lower==F(13,12)
assert early_ratio_lower>1
moments={
 'mass':sum(c for t,z,c in terms),
 'E_time':sum(c*t for t,z,c in terms),
 'E_space_squared':sum(c*z for t,z,c in terms),
 'E_space_fourth':sum(c*3*z*z for t,z,c in terms),
}
assert moments['mass']==1
assert moments['E_time']==moments['E_space_squared']==F(2,3)
fourth_gap=moments['E_space_fourth']-3*moments['E_space_squared']**2
assert fourth_gap==F(11,12)
# Hermite caloric moment equations to order 40 from exact exponential-rate cancellation.
for n in range(1,21):
    assert sum(c*(z-t)**n for t,z,c in terms)==0
certificate=dict(parameters={k:str(z) for k,z in dict(s=s,v=v,w=w,u=u).items()},
 transform_rates={str(k):str(z) for k,z in exponents.items()},
 terminal_variances={str(k):str(z) for k,z in terminal.items()},
 late_max_ratio_squared=str(late_ratio_squared),early_center_ratio_squared=str(early_ratio_squared),
 early_interval=['-1/4','1/4'],early_ratio_lower=str(early_ratio_lower),
 moments={k:str(z) for k,z in moments.items()},non_gaussian_fourth_gap=str(fourth_gap),
 caloric_moment_orders_checked=list(range(2,41,2)),
 scope='general positive time-space measure; NOT a deterministic graph; NOT Gaussian spatial marginal')
result=dict(certificate=certificate,runtime=dict(python=platform.python_version(),machine=platform.machine(),
 node=platform.node(),pid=os.getpid(),executable=sys.executable),
 script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
if len(sys.argv)>1:Path(sys.argv[1]).write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,indent=2))
