"""Independent SymPy checks for the analytic near lemma/shared hinge, no sampling."""
from pathlib import Path
import sympy as S
import json,sys,time,os,platform,hashlib
start=time.monotonic();ip=Path(sys.argv[1]);op=Path(sys.argv[2]);cfg=json.loads(ip.read_text());rat=lambda key:S.Rational(cfg[key])
m=rat('m');d=rat('reference_d');e=rat('reference_crop_e');target=rat('near_target');D=rat('far_d_min');delta=rat('distance_delta')
assert m==S.Rational(1,4)and d==S.Rational(71,200)and e==S.Rational(1,1000)
s=S.symbols('s',positive=True);x=S.symbols('x',nonnegative=True);v=S.Rational(1,4);u=(s+v*v/s)/2;z=(s-v*v/s)/2
assert S.factor(u*u-z*z)==v*v
C0=m/2+d*v-v*e+e*e;L=d+v+e
full1=(2*u-d)*(u-e)+(d-u-e)*(v-e)
assert S.factor(full1-(C0+2*z*z-L*u))==0
contact=S.factor(full1.subs(s,d));assert contact==S.Rational(691638777,10082000000)<target
G2=S.factor(full1-(s-d)*(z-e));expected2=S.Rational(22737,125000)-S.Rational(249,2000)*s-S.Rational(961,32000)/s+1/(512*s*s)
G3=S.factor((u-e)**2-(s-d)*(z-e));expected3=-s*s/4+S.Rational(71,400)*s+S.Rational(31073,500000)-S.Rational(357,32000)/s+1/(1024*s*s)
assert S.factor(G2-expected2)==0 and S.factor(G3-expected3)==0
# Formula matching at u=d-e: subtract the two uncropped formulas, factor,
# and use u=d-e directly rather than sampling a mass value.
U,Z=S.symbols('U Z',nonnegative=True)
assert S.factor(((2*U-d)*(U-e)+(d-U-e)*(v-e)-(U-e)**2).subs(U,d-e))==0
P=S.cancel((S.Rational(7,100)-G2)*s*s);assert S.denom(P)==1 or S.denom(P).is_Integer
P0=P.subs(s,d);P1=S.diff(P,s).subs(s,d);P2=S.diff(P,s,2).subs(s,d)
assert (P0,P1,P2)==(S.Rational(14101223,80000000000),-S.Rational(937829,400000000),S.Rational(41393,1000000))
assert S.factor(P.subs(s,d+x)-(P0+P1*x+P2*x*x/2+S.Rational(249,2000)*x**3))==0
lower=S.factor(P0-P1*P1/(2*P2));assert lower==S.Rational(291048892263,2649152000000000)>S.Rational(1,10000)
q=S.Rational(3,5);derivative_bound=(d-q)/2+(d+2*e)*v*v/(2*q*q)
assert derivative_bound==-S.Rational(1757,19200)<0 and u.subs(s,q)==S.Rational(169,480)<d-e
high_bound=G3.subs(s,q);assert high_bound==S.Rational(9038149,144000000)<target
assert (v-e)**2==S.Rational(62001,1000000)<target
B0=D*D-m;kappa=S.Rational(1,2)-D/4
hinges=[min(rat('total_outside_cap'),k*rat('corner_mass_cap'))-k*B0 for k in range(5)]
assert hinges==[0,S.Rational(111,400),S.Rational(111,200),S.Rational(233,400),S.Rational(11,100)]
assert max(hinges)==S.Rational(233,400)and hinges.index(max(hinges))==3
joint=S.factor(kappa*max(hinges));assert joint==S.Rational(5359,32000)<S.Rational(21,125)
assert S.Rational(7,16)**2>S.Rational(3,16)and S.Rational(3,4)**2>S.Rational(1,2)and S.Rational(19,16)<S.Rational(6,5)
crop=2*S.Rational(1,2)*e-e*e;near=target+crop+S.Rational(6,5)*delta;assert near==S.Rational(74499,1000000)<S.Rational(3,40)
assert 3*425**2-736**2==179>0 and S.Rational(3,8)-S.Rational(17,80)*S.Rational(736,425)==S.Rational(7,1000)
assert 4*S.Rational(3,40)+4*S.Rational(7,1000)+S.Rational(21,125)==S.Rational(62,125)
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
out={'status':'PASS','identities':{'full_union_first_branch':str(S.factor(full1)),'G2':str(S.apart(G2,s)),'G3':str(S.apart(G3,s)),'coverage_change_matches':True,'P_taylor_exact':True},'constants':{'zero_near_area':str((v-e)**2),'contact':str(contact),'taylor_lower':str(lower),'high_branch_derivative_upper':str(derivative_bound),'high_branch_cost_upper':str(high_bound),'active_count_hinge_bounds':list(map(str,hinges)),'hinge_max_active_count':3,'joint_sum_upper':str(joint),'robust_near_upper_using_exact_crop':str(near),'aggregate_relaxed_upper':'62/125'},'input_sha256':sha(ip),'script_sha256':sha(__file__),'scope':cfg['scope']+' Analytic branch domains and Borel bundle set inclusions require separate mathematical review, not these identities.','runtime':{'pid':os.getpid(),'python':platform.python_version(),'sympy':S.__version__,'seconds':time.monotonic()-start}}
op.write_text(json.dumps(out,indent=2)+'\n');print(json.dumps({'status':'PASS','constants':out['constants'],'runtime':out['runtime']}))
