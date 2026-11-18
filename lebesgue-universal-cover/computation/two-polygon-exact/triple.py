"""Three P-edge support contacts: exact quadratic, Arb root and actual pose."""
import json,math
from pathlib import Path
from fractions import Fraction as F
from flint import arb,ctx
ctx.prec=192;p=Path(__file__).resolve().parent;d=json.loads((p/'micro-scan.json').read_text());H=[tuple(map(F,x)) for x in d['hullP']];Q=[tuple(map(F,x)) for x in d['hullQ']];A=lambda x:arb(str(x));dot=lambda u,v:u[0]*v[0]+u[1]*v[1];cross=lambda u,v:u[0]*v[1]-u[1]*v[0];n=[(b[1]-a[1],a[0]-b[0]) for a,b in zip(H,H[1:]+H[:1])];ids=[0,7,26];th=.216339;co,si=math.cos(th),math.sin(th)
def rotfloat(v):return (co*float(v[0])-si*float(v[1]),si*float(v[0])+co*float(v[1]))
js=[max(range(len(Q)),key=lambda j:float(dot(n[i],rotfloat(Q[j])))) for i in ids];ns=[n[i] for i in ids];ww=[cross(ns[1],ns[2]),cross(ns[2],ns[0]),cross(ns[0],ns[1])];alpha=beta=gamma=F(0)
for k in range(3):
 nx,ny=ns[k];qx,qy=Q[js[k]];alpha+=ww[k]*(nx*qx+ny*qy);beta+=ww[k]*(-nx*qy+ny*qx);gamma-=ww[k]*dot(ns[k],H[ids[k]])
# alpha cos(theta)+beta sin(theta)+gamma=0; u=tan(theta/2).
a=gamma-alpha;b=2*beta;c=gamma+alpha;disc=A(b*b-4*a*c);assert disc>0;roots=[(-A(b)+sg*disc.sqrt())/(2*A(a)) for sg in (-1,1)];u,=[x for x in roots if x>A('.108') and x<A('.109')];den=1+u*u;ct=(1-u*u)/den;st=2*u/den;rot=lambda x:(ct*A(x[0])-st*A(x[1]),st*A(x[0])+ct*A(x[1]));h=[tuple(map(A,x)) for x in H];q=[rot(x) for x in Q];nA=[tuple(map(A,x)) for x in n];nsA=[nA[i] for i in ids];bs=[];support_margin=[]
for k in range(3):
 v=dot(nsA[k],q[js[k]]);diff=[v-dot(nsA[k],other) for j,other in enumerate(q) if j!=js[k]];assert all(x>0 for x in diff);support_margin.append(str(min(x.lower() for x in diff)));bs.append(v-dot(nsA[k],h[ids[k]]))
det=cross(nsA[0],nsA[1]);assert det.abs_lower()>0;t=((-bs[0]*nsA[1][1]+bs[1]*nsA[0][1])/det,(-nsA[0][0]*bs[1]+nsA[1][0]*bs[0])/det);assert (bs[2]+dot(nsA[2],t)).contains(0)
alln=nA+[rot((b[1]-a[1],a[0]-b[0])) for a,b in zip(Q,Q[1:]+Q[:1])];area=sum(F(x*y1-y*x1,2) for (x,y),(x1,y1) in zip(H,H[1:]+H[:1]));terms=[]
for nn in alln:
 hq=max(dot(nn,x).upper() for x in q);hh=max(dot(nn,x).lower() for x in h);terms.append(max(A(0),(hq-hh+dot(nn,t)).upper()))
up=(A(area)+sum(terms,A(0))/2).upper();out={'edge_indices':ids,'Q_support_vertex_indices':js,'tan_half_root_interval':str(u),'quadratic':[str(a),str(b),str(c)],'translation_intervals':[str(x) for x in t],'support_margin_lower':support_margin,'area_upper':str(up),'claim':'actual 3-contact pose, not proof of global optimality'};(p/'triple.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
