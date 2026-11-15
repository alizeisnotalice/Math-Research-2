from fractions import Fraction as F
import math,json,sys
N=128
b=F(1);S=F(0)
for n in range(N):
 S+=b*b/F((4*n+1)**2)
 b*=F(2*n+1,2*n+2)
assert S>1
I=math.gamma(.25)*math.gamma(.5)/(4*math.gamma(.75))
rows=[]
for q in [.1,.25,.5,.75,.9,.99]:
 b=1.;H=0.
 for n in range(10000):
  term=b*q**(4*n+1)/(4*n+1);H+=term
  b*=(2*n+1)/(2*n+2)
  if abs(term)<1e-17: break
 t=q*q;B=H*H*(-math.log(q))
 B1=.5*math.sqrt(t*(1-t)/(1+t))
 B2=(math.sqrt(2)-1)/2
 assert B<B1<=B2
 rows.append(dict(q=q,H=H,R=I/H,denominator=B,first_upper=B1,uniform_upper=B2))
old=F(3590865312933,10**12)
assert F(18,5)>old
assert F(7,5)**2<2
out=dict(python=sys.version,scope="independent scalar sanity only; analytic proof carries continuum",N_terms=N,S_partial_float=float(S),S_partial_gt_one=True,I_float=I,tau_over_R2_float=float(S)/(4*I*I),universal_exact_lower="pi*(1+sqrt(2))/2",universal_lower_float=math.pi*(1+math.sqrt(2))/2,rational_lower="18/5",rational_gap_over_old=str(F(18,5)-old),samples=rows)
print(json.dumps(out,indent=2))
