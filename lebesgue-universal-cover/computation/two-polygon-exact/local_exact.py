"""Bound each certified dual objective derivative over an entire angular cell."""
from pathlib import Path
from fractions import Fraction as F
from flint import arb,ctx
import json,math
ctx.prec=160;W=Path(__file__).resolve().parent;scan=json.loads((W/'fine-scan.json').read_text());cert=json.loads((W/'fine-certificate.json').read_text());H=[tuple(map(F,p)) for p in scan['hullP']];Q=[tuple(map(F,p)) for p in scan['hullQ']];N=scan['N'];assert N==1000000
A=lambda x:arb(str(x))
class J:
 def __init__(self,v,d=0):self.v=v if isinstance(v,arb) else A(v);self.d=d if isinstance(d,arb) else A(d)
 def __add__(a,b):
  if not isinstance(b,J):b=J(b)
  return J(a.v+b.v,a.d+b.d)
 __radd__=__add__
 def __neg__(a):return J(-a.v,-a.d)
 def __sub__(a,b):return a+(-b if isinstance(b,J) else J(-b))
 def __rsub__(a,b):return J(b)-a
 def __mul__(a,b):
  if not isinstance(b,J):b=J(b)
  return J(a.v*b.v,a.d*b.v+a.v*b.d)
 __rmul__=__mul__
 def __truediv__(a,b):
  if not isinstance(b,J):b=J(b)
  return J(a.v/b.v,(a.d*b.v-a.v*b.d)/(b.v*b.v))
dot=lambda u,v:u[0]*v[0]+u[1]*v[1];cross=lambda u,v:u[0]*v[1]-u[1]*v[0]
def ns(p):return [(b[1]-a[1],a[0]-b[0]) for a,b in zip(p,p[1:]+p[:1])]
hn,qn=ns(H),ns(Q);M=len(H);L=len(Q);h=[tuple(J(x) for x in p) for p in H];q=[tuple(J(x) for x in p) for p in Q];hnj=[tuple(J(x) for x in n) for n in hn];pi=arb.pi();half=pi/A(N);upper=arb(cert['feasible_upper']);err=arb(cert['angle_error']);target=upper-arb('1e-9');results=[];failed=[]

# The three-contact root is an exact isolated root of the rational quadratic.
tri=json.loads((W/'triple.json').read_text());u=arb(tri['tan_half_root_interval']);root=2*u.atan();assert root>pi*A(F(2*344314,10000000)) and root<pi*A(F(2*344315,10000000))
# Verify no extra zero hinge and complementarity at the exact three-contact pose.
t0=[arb(x) for x in tri['translation_intervals']];cc=(1-u*u)/(1+u*u);ss=2*u/(1+u*u);rotate=lambda x:(cc*A(x[0])-ss*A(x[1]),ss*A(x[0])+cc*A(x[1]));h0=[tuple(map(A,x)) for x in H];q0=[rotate(x) for x in Q];n0=[tuple(map(A,x)) for x in hn]+[rotate(x) for x in qn]
def da(x,y):return x[0]*y[0]+x[1]*y[1]
def mx(z):return max(x.lower() for x in z).union(max(x.upper() for x in z))
ones={8,9,10,15,16,17,20,27,28,29,33};contacts={0,7,26};margin=[]
for k,nn in enumerate(n0):
 z=mx([da(nn,x) for x in q0])-mx([da(nn,x) for x in h0])+da(nn,t0)
 if k in contacts:assert z.contains(0)
 elif k in ones:assert z>0,(k,z);margin.append(z.lower())
 else:assert z<0,(k,z);margin.append((-z).lower())
records=json.loads((W/'micro-scan.json').read_text())['results'];proof=[]
for name,rec,edge in [('left',next(x for x in records if x['i']==344313),pi*A(F(2*344314,10000000))),('right',next(x for x in records if x['i']==344315),pi*A(F(2*344315,10000000)))]:
 theta=edge.union(root);c=J(theta.cos(),-theta.sin());s=J(theta.sin(),theta.cos());rot=lambda p:(c*p[0]-s*p[1],s*p[0]+c*p[1]);qa=[rot(p) for p in q];n=hnj+[rot(tuple(J(x) for x in v)) for v in qn]
 lam=[F(str(round(x,12))) for x in rec['l']];lam=[F(0) if abs(x)<F(1,10**9) else F(1) if abs(x-1)<F(1,10**9) else x for x in lam];piv=[k for k,x in enumerate(lam) if 0<x<1];assert piv==([7,26] if name=='left' else [0,7]);assert {i for i,x in enumerate(lam) if x==1}==ones
 weights=[J(x) for x in lam];wh=[sum((lam[k]*hn[k][j] for k in range(M)),F(0)) for j in (0,1)];wq=[sum((lam[M+k]*qn[k][j] for k in range(L)),F(0)) for j in (0,1)];a,b=piv;den=cross(n[a],n[b]);res=(J(wh[0])+c*J(wq[0])-s*J(wq[1]),J(wh[1])+s*J(wq[0])+c*J(wq[1]));weights[a]-=cross(res,n[b])/den;weights[b]-=cross(n[a],res)/den
 assert weights[a].v>0 and weights[a].v<1 and weights[b].v>0 and weights[b].v<1
 derivative=arb(0);minsupport=arb(1)
 for k in sorted(ones|set(piv)):
  nn=n[k]
  if k<M:
   vals=[dot(nn,qq) for qq in qa];j=max(range(L),key=lambda j:float(vals[j].v.mid()));diff=[(vals[j]-vals[l]).v for l in range(L) if l!=j];assert all(x>0 for x in diff),(name,k,'Q support');minsupport=min(minsupport,min(x.lower() for x in diff));hh=dot(nn,h[k]);hq=vals[j]
  else:
   hq=dot(nn,qa[k-M]);vals=[dot(nn,hh) for hh in h];j=max(range(M),key=lambda j:float(vals[j].v.mid()));diff=[(vals[j]-vals[l]).v for l in range(M) if l!=j];assert all(x>0 for x in diff),(name,k,'P support');minsupport=min(minsupport,min(x.lower() for x in diff));hh=vals[j]
  derivative+=(weights[k]*(hq-hh)).d/2
 assert (derivative<0 if name=='left' else derivative>0),(name,derivative)
 proof.append({'side':name,'angle':str(theta),'pivots':piv,'derivative':str(derivative),'weight_intervals':[str(weights[k].v) for k in piv],'min_support_gap':str(minsupport)})
out={'root_angle':str(root),'hinge_nonzero_margin_lower':str(min(margin)),'contacts':sorted(contacts),'positive_hinges':sorted(ones),'proof':proof,'interval_cell':344314,'exact_contact_algebra':'triple.py','precision_bits':ctx.prec};(W/'local-exact.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
