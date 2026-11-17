"""Standard-library verifier; no optimizer output is trusted."""
import json,hashlib,sys
from pathlib import Path
from fractions import Fraction as F
from rational_trig import PI,scale,add,mul,inv,trig
base=Path(__file__).parent
raw=(base/'004-remote-certificate.json').read_bytes();cover=json.loads(raw);c=json.loads((base/'005-obstruction-certificate.json').read_text())
assert hashlib.sha256(raw).hexdigest()==c['cover_source_sha256']
assert c['n']==11 and c['M']==180 and sorted(r['cell'] for r in c['rows'])==list(range(180))
rad=scale(inv(trig(scale(PI,F(1,22)),'cos')),F(1,2));eps=F(cover['epsilon']);S=cover['coordinate_denominator'];worst=None
for row in c['rows']:
 j=row['cell'];q=row['normals'];w=list(map(F,row['weights']));v=row['vertices']
 assert len(q)==len(w)==len(v)==3 and all(x>=0 for x in w) and sum(w)==1
 assert all(sum(w[i]*q[i][k] for i in range(3))==0 for k in (0,1))
 gap=F(0)
 for normal,weight,vertex,bs in zip(q,w,v,row['b']):
  assert 0<=vertex<11 and len(normal)==2 and all(type(x)==int for x in normal)
  b=F(max(normal[0]*p[0]+normal[1]*p[1] for p in cover['hull_grid_points']),S)+eps*(abs(normal[0])+abs(normal[1]))
  assert b==F(bs)
  a=F(2*j,11*180)+F(2*vertex,11);d=F(1,11*180)
  interval=mul(PI,(a-d,a+d));x=mul(rad,trig(interval,'cos'));y=mul(rad,trig(interval,'sin'))
  dot=add(scale(x,F(normal[0])),scale(y,F(normal[1])))[0]
  gap+=weight*(b-dot)
 assert gap==F(row['gap_upper']) and gap<0
 worst=gap/10**6 if worst is None else max(worst,gap/10**6)
assert worst==F(c['worst_gap_normalized'])
print(json.dumps({'python':sys.version,'verified_cells':180,'worst_gap_normalized':str(worst),'all_strictly_negative':True,'status':'passed','scope':'all translations, rotations and reflections of regular diameter-one 11-gon excluded from fixed cover'},indent=2))
