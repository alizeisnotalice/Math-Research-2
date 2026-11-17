"""Exact all-angle exclusion for h=1/2+cos(3t)/80-19cos(5t)/1200.
LP is used only to select witnesses. All acceptance tests use Fraction.
"""
from pathlib import Path
from fractions import Fraction as F
import json,math,sys,hashlib
import numpy as np
from scipy.optimize import linprog
from rational_trig import PI,scale,add,neg,mul,trig
A=F(1,80);B=F(19,1200);SAMPLES=64;M=2048
assert 8*A+24*B==F(12,25)
def material(k):
 def t(harmonic,kind):return trig(scale(PI,F(2*k*harmonic,SAMPLES)%2),kind)
 h=add((F(1,2),F(1,2)),add(scale(t(3,'cos'),A),scale(t(5,'cos'),-B)))
 hp=add(scale(t(3,'sin'),-3*A),scale(t(5,'sin'),5*B))
 x=add(mul(h,t(1,'cos')),neg(mul(hp,t(1,'sin'))));y=add(mul(h,t(1,'sin')),mul(hp,t(1,'cos')))
 return x,y
P=[material(k) for k in range(SAMPLES)];Pf=np.array([[float(sum(x)/2),float(sum(y)/2)] for x,y in P])
cover_path=Path(__file__).with_name('010-seven-body-certificate.json');raw=cover_path.read_bytes();r=json.loads(raw);V=np.array(r['hull_grid_points'],dtype=np.int64);S=r['coordinate_denominator'];eps=F(r['epsilon'])
e=np.roll(V,-1,axis=0)-V;N=np.column_stack((e[:,1],-e[:,0])).astype(float);N/=np.linalg.norm(N,axis=1)[:,None];N=np.rint(N*10**6).astype(np.int64)
b=[F(int((V@q).max()),S)+eps*(abs(int(q[0]))+abs(int(q[1]))) for q in N]
rows=[];worst=None
for j in range(M):
 theta=2*math.pi*j/M;cs,sn=math.cos(theta),math.sin(theta);Q=Pf@np.array([[cs,sn],[-sn,cs]])
 values=(N/10**6)@Q.T;ids=values.argmax(axis=1);h=values.max(axis=1)
 z=linprog([0,0,1],A_ub=np.column_stack((N/10**6,-np.ones(len(N)))),b_ub=np.array([float(v)/10**6 for v in b])-h,bounds=[(None,None)]*3,method='highs');assert z.success
 active=np.argsort(-np.abs(z.ineqlin.marginals))[:3].tolist();q=[tuple(map(int,N[i])) for i in active]
 det=lambda u,v:u[0]*v[1]-u[1]*v[0]
 w=[det(q[1],q[2]),det(q[2],q[0]),det(q[0],q[1])]
 if sum(w)<0:w=[-x for x in w]
 assert min(w)>=0 and sum(w)>0;w=[F(x,sum(w)) for x in w]
 assert all(sum(w[i]*q[i][k] for i in range(3))==0 for k in (0,1))
 interval=mul(PI,(F(2*j-1,M),F(2*j+1,M)));co=trig(interval,'cos');si=trig(interval,'sin');gap=F(0);vertices=[]
 for i,normal,weight in zip(active,q,w):
  k=int(ids[i]);vertices.append(k);x,y=P[k]
  rx=add(mul(x,co),neg(mul(y,si)));ry=add(mul(x,si),mul(y,co))
  low=add(scale(rx,F(normal[0])),scale(ry,F(normal[1])))[0]
  gap+=weight*(b[i]-low)
 assert gap<0,(j,float(gap/10**6))
 worst=gap/10**6 if worst is None else max(worst,gap/10**6)
 rows.append({'cell':j,'normals':q,'weights':[str(x) for x in w],'vertices':vertices,'b':[str(b[i]) for i in active],'gap_upper':str(gap)})
print(json.dumps({'shape':'h(t)=1/2+cos(3t)/80-19cos(5t)/1200','M':M,'samples':SAMPLES,'cover_source_sha256':hashlib.sha256(raw).hexdigest(),'worst_gap_normalized':str(worst),'worst_gap_display':float(worst),'rows':rows,'python':sys.version,'status':'passed','scope':'fixed seven-body cover not universal; no general area lower bound'},indent=2))
