"""Bound each certified dual objective derivative over an entire angular cell."""
from pathlib import Path
from fractions import Fraction as F
from flint import arb,ctx
import json,math
ctx.prec=160;W=Path(__file__).resolve().parent;scan=json.loads((W/'micro-scan.json').read_text());cert=json.loads((W/'micro-certificate.json').read_text());H=[tuple(map(F,p)) for p in scan['hullP']];Q=[tuple(map(F,p)) for p in scan['hullQ']];N=scan['N'];assert N==10000000
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
for z,lo in zip(scan['results'],cert['cell_lowers']):
 i=z['i'];center=pi*A(F(2*i+1,N));theta=center+arb(0,half);c=J(theta.cos(),-theta.sin());s=J(theta.sin(),theta.cos());rot=lambda p:(c*p[0]-s*p[1],s*p[0]+c*p[1]);qa=[rot(p) for p in q];n=hnj+[rot(tuple(J(x) for x in v)) for v in qn];lam=[F(str(round(x,12))) for x in z['l']];lam=[F(0) if abs(x)<F(1,10**9) else F(1) if abs(x-1)<F(1,10**9) else x for x in lam];w=[J(x) for x in lam];wh=[sum((lam[k]*hn[k][j] for k in range(M)),F(0)) for j in (0,1)];wq=[sum((lam[M+k]*qn[k][j] for k in range(L)),F(0)) for j in (0,1)]
 if any(wh+wq):
  inds=[k for k,x in enumerate(lam) if 0<x<1];assert len(inds)==2,(i,inds);pairs=[(abs(float(cross(n[a],n[b]).v.mid())),a,b) for a in inds for b in inds if a<b];
  if not pairs:failed.append([i,'pivots']);continue
  _,a,b=max(pairs);d=cross(n[a],n[b]);res=(J(wh[0])+c*J(wq[0])-s*J(wq[1]),J(wh[1])+s*J(wq[0])+c*J(wq[1]));w[a]=w[a]-cross(res,n[b])/d;w[b]=w[b]-cross(n[a],res)/d
  if not (w[a].v>0 and w[a].v<1 and w[b].v>0 and w[b].v<1):failed.append([i,'range']);continue
 deriv=arb(0); q_choice=[]
 for k,nn in enumerate(n):
  if k<M:
   # Any selected Q vertex is a valid lower support globally.
   j=max(range(L),key=lambda j:float(dot(nn,qa[j]).v.mid()));q_choice.append(j);v=dot(nn,qa[j])-dot(nn,h[k]);deriv+=(w[k]*v).d
  else:
   # Q edge support is its own first endpoint, even at parallel-edge angles.
   q_choice.append(k-M);hq=dot(nn,qa[k-M]);ps=[dot(nn,x) for x in h];best=max(x.v.lower() for x in ps);choices=[j for j,x in enumerate(ps) if x.v.upper()>=best];ds=[(w[k]*(hq-ps[j])).d for j in choices];deriv+=min(x.lower() for x in ds).union(max(x.upper() for x in ds))
 # Original midpoint bound used maximal Q support. Replace it with the exact
 # chosen vertices via a nonnegative rigorous support-gap correction.
 cm,sm=center.cos(),center.sin();rm=lambda x:(cm*A(x[0])-sm*A(x[1]),sm*A(x[0])+cm*A(x[1]));qm=[rm(x) for x in Q];nm=[tuple(map(A,x)) for x in hn]+[rm(x) for x in qn];gap=arb(0)
 for k,nn in enumerate(nm):
  allq=[nn[0]*qq[0]+nn[1]*qq[1] for qq in qm];gap+=max(A(0),(max(x.upper() for x in allq)-allq[q_choice[k]].lower()).upper())
 dmax=(abs(deriv)/2).upper();lower=(arb(lo).lower()-gap/2-dmax*half).lower();results.append({'i':i,'derivative_abs_upper':str(dmax),'lower':str(lower)})
assert not failed,failed[:20]
print('selected',len(results),'best derivative max',max((arb(r['derivative_abs_upper']).upper() for r in results)),'min lower',min(arb(r['lower']).lower() for r in results),flush=True)
(W/'micro-derivative-cells.json').write_text(json.dumps({'N':N,'results':results,'failed':failed},indent=2)+'\n')
