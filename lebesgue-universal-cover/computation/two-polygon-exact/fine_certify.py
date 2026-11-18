"""Full-angle balanced dual certificates for two fixed rational polygons."""
import sys,json,math,hashlib,time
from pathlib import Path
from fractions import Fraction as F
from flint import arb,ctx
ctx.prec=160
W=Path(__file__).resolve().parent;data=json.loads((W/'fine-scan.json').read_text());H=[tuple(map(F,p)) for p in data['hullP']];Q=[tuple(map(F,p)) for p in data['hullQ']];N=data['N'];records=data['results'];selected=json.loads((W/'selected.json').read_text());assert [r['i'] for r in records]==[100*i+j for i in selected for j in range(100)] and N==1000000
A=lambda x: arb(str(x)); dot=lambda u,v:u[0]*v[0]+u[1]*v[1]; cross=lambda u,v:u[0]*v[1]-u[1]*v[0]
def normals(p):return [(b[1]-a[1],a[0]-b[0]) for a,b in zip(p,p[1:]+p[:1])]
hn,qn=normals(H),normals(Q);M=len(hn);L=len(qn);assert len(H)==30 and len(Q)==11
inp=json.loads((W/'inputs.json').read_text()); pp=inp['P']; Pin=[(F(x,pp['denominator']),F(y,pp['denominator'])) for x,y in pp['points']]; Qin=[tuple(map(F,x)) for x in inp['Q']['points']]; assert set(H)==set(Pin) and set(Q)==set(Qin)
assert set(H)=={(x,-y) for x,y in H}
area=sum(cross(a,b) for a,b in zip(H,H[1:]+H[:1]))/2
for p in (H,Q):
 for a,b in zip(p,p[1:]+p[:1]):
  n=(b[1]-a[1],a[0]-b[0]);assert all(dot(n,x)<=dot(n,a) for x in p)
# Exact inclusion of radius .46 disk in H and .58 disk around the origin containing Q.
assert all((h:=dot(n,x))>0 and h*h>=F(46,100)**2*dot(n,n) for n,x in zip(hn,H))
assert all(dot(q,q)<=F(58,100)**2 for q in Q)
hA=[tuple(map(A,x)) for x in H];qA=[tuple(map(A,x)) for x in Q];hnA=[tuple(map(A,x)) for x in hn];qh=[dot(n,x) for n,x in zip(hn,H)];pi=arb.pi();lowest=None;bestupper=None;bestidx=None;failed=[];cell_lowers=[];best_scan_idx=min(records,key=lambda r:r['v'])['i'];t0=time.time()
for rec in records:
 i=rec['i'];theta=pi*A(F(2*i+1,N));c=theta.cos();s=theta.sin()
 rot=lambda x:(c*x[0]-s*x[1],s*x[0]+c*x[1])
 qa=[rot(q) for q in qA];ns=hnA+[rot(tuple(map(A,n))) for n in qn]
 # Any particular Q vertex provides a lower support bound; all H vertices give an upper bound.
 vals=[]
 for k,n in enumerate(ns):
  lower=max((dot(n,q).lower() for q in qa))
  upper=A(qh[k]).upper() if k<M else max((dot(n,h).upper() for h in hA))
  vals.append((lower-upper).lower())
 lam=[F(str(round(x,12))) for x in rec['l']]
 lam=[F(0) if abs(x)<F(1,10**9) else F(1) if abs(x-1)<F(1,10**9) else x for x in lam]
 assert all(0<=x<=1 for x in lam)
 wh=[sum((lam[k]*hn[k][j] for k in range(M)),F(0)) for j in (0,1)];wq=[sum((lam[M+k]*qn[k][j] for k in range(L)),F(0)) for j in (0,1)]
 lb=[A(x) for x in lam]
 if any(wh+wq):
  inds=[k for k,x in enumerate(lam) if 0<x<1];best=None
  for a in inds:
   for b in inds:
    if b<=a:continue
    det=cross(ns[a],ns[b]);score=abs(float(det.mid()))
    if best is None or score>best[0]:best=(score,a,b,det)
  if not best or not best[3].abs_lower()>0:failed.append((i,'pivots'));continue
  _,a,b,d=best;res=(A(wh[0])+c*A(wq[0])-s*A(wq[1]),A(wh[1])+s*A(wq[0])+c*A(wq[1]));lb[a]-=cross(res,ns[b])/d;lb[b]-=cross(ns[a],res)/d
  if not (lb[a]>0 and lb[a]<1 and lb[b]>0 and lb[b]<1):failed.append((i,'weights'));continue
 low=(A(area)+sum((x*y for x,y in zip(lb,vals)),A(0))/2).lower()
 cell_lowers.append(str(low))
 if lowest is None or low<lowest:lowest,bestidx=low,i
 # Fixed floating-point translation rounded to a rational gives an actual feasible pose.
 if i==best_scan_idx:
  t=[A(F(str(round(v,12)))) for v in rec['t']];term=[]
  for k,n in enumerate(ns):
   hq=max((dot(n,q).upper() for q in qa));hh=A(qh[k]).lower() if k<M else max((dot(n,h).lower() for h in hA));v=hq-hh+dot(n,t);term.append(max(A(0),v.upper()))
  bestupper=(A(area)+sum(term,A(0))/2).upper()
 if i%500==0:print(i,'lowest',lowest,'failed',len(failed),'seconds',time.time()-t0,flush=True)
assert not failed,failed[:20]
assert bestupper<A('.761')
# K contains disk radius .46.  For its area <= .761, perimeter <= 2*.761/.46.
# Rotation displacement <= .58*pi/N.  Steiner then controls area at nearest midpoint.
delta=A('.58')*pi/A(N);error=(A(2)*A('.761')/A('.46')*delta+pi*delta*delta).upper();lower=(lowest-error).lower()
out={'cell_lowers':cell_lowers,'N':N,'angle_domain':'[0,2pi]','hands':'P30 reflection invariance reduces Q to one hand','precision_bits':ctx.prec,'midpoint_lower':str(lowest),'midpoint_lower_index':bestidx,'feasible_upper':str(bestupper),'angle_error':str(error),'global_lower':str(lower),'failed':failed,'P_inradius':'0.46','Q_radius':'0.58','areaH':str(area),'source_scan_sha256':hashlib.sha256((W/'fine-scan.json').read_bytes()).hexdigest(),'elapsed_seconds':time.time()-t0};(W/'fine-certificate.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps({k:v for k,v in out.items() if k!='cell_lowers'},indent=2),flush=True)
