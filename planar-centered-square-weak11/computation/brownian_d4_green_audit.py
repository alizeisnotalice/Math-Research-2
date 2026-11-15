from fractions import Fraction as Q
import math,sys,json
def add(*polys):
 out=[Q(0)]*max(map(len,polys))
 for p in polys:
  for i,x in enumerate(p):out[i]+=x
 return out
def mul(a,b):
 out=[Q(0)]*(len(a)+len(b)-1)
 for i,x in enumerate(a):
  for j,y in enumerate(b):out[i+j]+=x*y
 return out
p1=[-Q(1,3),Q(0),Q(1)];p2=[-Q(18,25),Q(1)]
sos=add(mul(p1,p1),[Q(5,3)*x for x in mul(p2,p2)],[Q(28,1125)])
assert sos==[Q(1),-Q(12,5),Q(1),Q(0),Q(1)]
rows=[]
for q in [Q(1,10),Q(1,4),Q(1,2),Q(3,4),Q(9,10),Q(99,100)]:
 for a in [[Q(1)],[Q(1),-Q(1,20),Q(1,100)],[Q(1,2),Q(1,3),-Q(1,5),Q(1,7)]]:
  value=sum(x*q**(4*n+1) for n,x in enumerate(a));norm=sum(x*x for x in a)
  assert norm*q*q/(1-q**8)-value*value>0
 qf=float(q);t=qf*qf
 original=math.pi*(1-qf**8)/(4*qf*qf*(-math.log(qf)))
 simplified=math.pi*(1+t+t*t+t*t*t)/(2*math.sqrt(t))
 assert original>simplified>6*math.pi/5
 rows.append(dict(q=str(q),necessary_cost_float=original,log_relaxed_float=simplified))
old=Q(3590865312933,10**12)
assert Q(18,5)>old
print(json.dumps(dict(python=sys.version,scope="exact rational SOS and finite Cauchy tests; continuum proof analytic",SOS_coefficients=list(map(str,sos)),SOS_positive_remainder="28/1125",Cauchy_exact_tests=18,universal_lower="6*pi/5",universal_lower_float=6*math.pi/5,rational_gap_over_old=str(Q(18,5)-old),samples=rows),indent=2))
