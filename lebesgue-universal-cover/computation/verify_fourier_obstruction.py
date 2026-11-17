"""Pure Fraction full-angle verification of the fixed Fourier witness."""
import json,hashlib,sys
from pathlib import Path
from fractions import Fraction as F
from rational_trig import PI,scale,add,neg,mul,trig
base=Path(__file__).parent;raw=(base/'010-seven-body-certificate.json').read_bytes();H=json.loads(raw);C=json.loads((base/'011-fourier-certificate.json').read_text())
assert C['cover_source_sha256']==hashlib.sha256(raw).hexdigest()
assert C['shape']=='h(t)=1/2+cos(3t)/80-19cos(5t)/1200'
M=C['M'];count=C['samples'];assert M==2048 and count==64
assert sorted(x['cell'] for x in C['rows'])==list(range(M))
points=[]
for k in range(count):
 angles={j:scale(PI,F(2*k*j,count)%2) for j in (1,3,5)}
 co={j:trig(angles[j],'cos') for j in angles};si={j:trig(angles[j],'sin') for j in angles}
 h=add((F(1,2),F(1,2)),add(scale(co[3],F(1,80)),scale(co[5],-F(19,1200))))
 hp=add(scale(si[3],-F(3,80)),scale(si[5],F(19,240)))
 points.append((add(mul(h,co[1]),neg(mul(hp,si[1]))),add(mul(h,si[1]),mul(hp,co[1]))))
eps=F(H['epsilon']);S=H['coordinate_denominator'];worst=None
for row in C['rows']:
 j=row['cell'];normals=row['normals'];weights=list(map(F,row['weights']));vertices=row['vertices']
 assert len(normals)==len(weights)==len(vertices)==3 and min(weights)>=0 and sum(weights)==1
 assert all(sum(w*q[k] for w,q in zip(weights,normals))==0 for k in (0,1))
 interval=mul(PI,(F(2*j-1,M),F(2*j+1,M)));co=trig(interval,'cos');si=trig(interval,'sin');gap=F(0)
 for q,w,k,stored_b in zip(normals,weights,vertices,row['b']):
  assert len(q)==2 and all(type(v)==int for v in q) and 0<=k<count
  b=F(max(q[0]*v[0]+q[1]*v[1] for v in H['hull_grid_points']),S)+eps*(abs(q[0])+abs(q[1]));assert b==F(stored_b)
  x,y=points[k];rx=add(mul(x,co),neg(mul(y,si)));ry=add(mul(x,si),mul(y,co))
  low=add(scale(rx,F(q[0])),scale(ry,F(q[1])))[0];gap+=w*(b-low)
 assert gap==F(row['gap_upper']) and gap<0
 worst=gap/10**6 if worst is None else max(worst,gap/10**6)
assert worst==F(C['worst_gap_normalized'])
print(json.dumps({'python':sys.version,'verified_cells':M,'samples':count,'worst_gap_normalized':str(worst),'status':'passed','scope':'fixed seven-body cover cannot contain this Fourier body under any Euclidean isometry'},indent=2))
