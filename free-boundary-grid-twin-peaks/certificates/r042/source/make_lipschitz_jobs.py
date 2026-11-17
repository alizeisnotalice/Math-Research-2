"""Generate all 178 n5 Lipschitz profiles and 4134 exact symmetry orbits."""
import itertools,json,sys
from pathlib import Path
n=5;h=[''.join(map(str,t)) for t in itertools.product(range(1,n),repeat=n) if all(abs(x-y)<=1 for x,y in zip(t,t[1:]))]
def ts(s):
 c=''.join(str(n-int(x)) for x in s);return s,s[::-1],c,c[::-1]
groups={}
for i,a in enumerate(h):
 for b in h[i:]:
  key=min(tuple(sorted(t)) for t in zip(ts(a),ts(b)));z=groups.setdefault(key,dict(curves=list(key),diagonal_multiplicity=0,offdiagonal_multiplicity=0,members=[]));z['diagonal_multiplicity' if a==b else 'offdiagonal_multiplicity']+=1;z['members'].append([a,b])
jobs=[dict(id=i,**v) for i,(k,v) in enumerate(sorted(groups.items()))]
assert len(h)==178 and len(jobs)==4134
out=Path(sys.argv[1]);assert not out.exists()
out.write_text(json.dumps(dict(n=n,curve_count=len(h),curves=h,jobs=jobs,symmetry='simultaneous row reversal and left-right reflection; unordered pairs'),indent=2)+'\n')
print('curves=178 orbit_jobs=4134 diagonal=178 offdiagonal=15753')
