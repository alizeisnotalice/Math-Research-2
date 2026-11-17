"""Standard-library verifier; no optimizer output is trusted."""
import json,hashlib,sys
from pathlib import Path
from fractions import Fraction as F
from rational_trig import PI,scale,add,mul,inv,trig
base=Path(__file__).parent
raw=Path(sys.argv[1]).read_bytes();cover=json.loads(raw);c=json.loads(Path(sys.argv[2]).read_text());n=c['n'];M=c['M']
assert hashlib.sha256(raw).hexdigest()==c['cover_source_sha256']
assert n>=3 and n%2==1 and M>0 and sorted(r['cell'] for r in c['rows'])==list(range(M))
rad=scale(inv(trig(scale(PI,F(1,2*n)),'cos')),F(1,2));eps=F(cover['epsilon']);S=cover['coordinate_denominator'];worst=None
for row in c['rows']:
 j=row['cell'];q=row['normals'];w=list(map(F,row['weights']));v=row['vertices']
 assert len(q)==len(w)==len(v)==3 and all(x>=0 for x in w) and sum(w)==1
 assert all(sum(w[i]*q[i][k] for i in range(3))==0 for k in (0,1))
 gap=F(0)
 for normal,weight,vertex,bs in zip(q,w,v,row['b']):
  assert 0<=vertex<n and len(normal)==2 and all(type(x)==int for x in normal)
  b=F(max(normal[0]*p[0]+normal[1]*p[1] for p in cover['hull_grid_points']),S)+eps*(abs(normal[0])+abs(normal[1]))
  assert b==F(bs)
  a=F(2*j,n*M)+F(2*vertex,n);d=F(1,n*M)
  interval=mul(PI,(a-d,a+d));x=mul(rad,trig(interval,'cos'));y=mul(rad,trig(interval,'sin'))
  dot=add(scale(x,F(normal[0])),scale(y,F(normal[1])))[0]
  gap+=weight*(b-dot)
 assert gap==F(row['gap_upper']) and gap<0
 worst=gap/10**6 if worst is None else max(worst,gap/10**6)
assert worst==F(c['worst_gap_normalized'])
print(json.dumps({'python':sys.version,'verified_cells':M,'n':n,'worst_gap_normalized':str(worst),'all_strictly_negative':True,'status':'passed','scope':'all translations, rotations and reflections of regular diameter-one n-gon excluded from fixed cover'},indent=2))
