"""Numerical dual selection followed by exact all-angle certificate."""
from pathlib import Path
import json,math,sys,hashlib
import numpy as np
from scipy.optimize import linprog
from fractions import Fraction as F
from rational_trig import PI,scale,add,mul,inv,trig
r=json.loads(Path(__file__).with_name('004-remote-certificate.json').read_text());V=np.array(r['hull_grid_points'],dtype=np.int64);S=r['coordinate_denominator'];eps=F(r['epsilon'])
edge=np.roll(V,-1,axis=0)-V;N=np.column_stack((edge[:,1],-edge[:,0])).astype(float);N/=np.linalg.norm(N,axis=1)[:,None];N=np.rint(N*10**6).astype(np.int64)
# Conservative rational supporting halfplanes for H+eps B.
b=[F(int((V@q).max()),S)+eps*(abs(int(q[0]))+abs(int(q[1]))) for q in N]
rows=[];worst=None
n=11;M=180;R=1/(2*math.cos(math.pi/(2*n)))
for j in range(M):
 angle=2*math.pi*j/(n*M);a=angle+2*math.pi*np.arange(n)/n
 pts=R*np.column_stack((np.cos(a),np.sin(a)));vals=(N.astype(float)/10**6)@pts.T;ids=vals.argmax(axis=1);h=vals.max(axis=1)
 z=linprog([0,0,1],A_ub=np.column_stack((N/10**6,-np.ones(len(N)))),b_ub=np.array([float(x)/10**6 for x in b])-h,bounds=[(None,None)]*3,method='highs');assert z.success
 active=np.argsort(-np.abs(z.ineqlin.marginals))[:3].tolist();q=[tuple(map(int,N[i])) for i in active]
 det=lambda u,v:u[0]*v[1]-u[1]*v[0]
 w=[det(q[1],q[2]),det(q[2],q[0]),det(q[0],q[1])]
 if sum(w)<0:w=[-x for x in w]
 assert all(x>=0 for x in w) and sum(w)>0
 weights=[F(x,sum(w)) for x in w]
 assert all(sum(weights[i]*q[i][k] for i in range(3))==0 for k in (0,1))
 rad=scale(inv(trig(scale(PI,F(1,2*n)),'cos')),F(1,2))
 terms=[];vertices=[]
 for i,normal in zip(active,q):
  vertex=int(ids[i]);vertices.append(vertex)
  # Full closed angle interval, including wrapping near theta=0.
  c=F(2*j,n*M)+F(2*vertex,n);d=F(1,n*M)
  angle_interval=mul(PI,(c-d,c+d))
  x=mul(rad,trig(angle_interval,'cos'));y=mul(rad,trig(angle_interval,'sin'))
  support_lower=add(scale(x,F(normal[0])),scale(y,F(normal[1])))[0]
  terms.append(b[i]-support_lower)
 gap=sum(wi*t for wi,t in zip(weights,terms));assert gap<0,(j,float(gap))
 normalized=gap/10**6;worst=normalized if worst is None else max(worst,normalized)
 rows.append({'cell':j,'normals':q,'weights':[str(x) for x in weights],'vertices':vertices,'b':[str(b[i]) for i in active],'gap_upper':str(gap)})
print(json.dumps({'n':n,'M':M,'worst_gap_normalized':str(worst),'worst_gap_display':float(worst),'cover_source_sha256':hashlib.sha256(Path(__file__).with_name('004-remote-certificate.json').read_bytes()).hexdigest(),'python':sys.version,'rows':rows,'status':'exact inequalities passed','scope':'regular diameter-one 11-gon cannot fit in fixed H+eps B; not universal lower bound'},indent=2))
