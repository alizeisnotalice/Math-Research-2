"""Certify R lies outside the C unit disk, and orientation of tangent path."""
from pathlib import Path
from fractions import Fraction as F
import runpy,io,contextlib,json,sys,hashlib
with contextlib.redirect_stdout(io.StringIO()):
 d=runpy.run_path(str(Path(__file__).with_name('certify_ah_path.py')))
I=d['I'];V=d['V'];add=d['add'];sub=d['sub'];mul=d['mul'];dot=d['dot'];norm2=d['norm2'];upper=d['upper_cross']
C,G,E,Nm,Y,Z,T=(d[k] for k in ('C','G','E','Nm','Y','Z','T'))
eL=d['dL'];lg=sub(G,C);ng=sub(E,Nm);length=dot(lg,eL)
assert length.lo>0
stack=[(F(0),F(1),F(0),F(1),0)];leaves=[];nodes=0
while stack:
 a,b,c,e,depth=stack.pop();nodes+=1
 u=I(a,b);v=I(c,e);L=add(C,mul(u,lg));N=add(Nm,mul(v,ng))
 X=upper(L,N);ell=u*length
 factor=2*dot(sub(X,L),eL)+ell
 if factor.lo>0:
  leaves.append({'u':[str(a),str(b)],'v':[str(c),str(e)],'factor_lower':str(factor.lo)});continue
 assert depth<22,'unresolved box; no proof produced'
 # Split by rational upper bound on physical span of its center segment.
 if (b-a)*F(15,1000)>(e-c)*F(26,100):
  m=(a+b)/2;stack.extend([(a,m,c,e,depth+1),(m,b,c,e,depth+1)])
 else:
  m=(c+e)/2;stack.extend([(a,b,c,m,depth+1),(a,b,m,e,depth+1)])
assert sum((F(v['u'][1])-F(v['u'][0]))*(F(v['v'][1])-F(v['v'][0])) for v in leaves)==1
# A DFS partition of the initial rectangle; no leaf is silently dropped.
def cross(a,b):return a[0]*b[1]-a[1]*b[0]
u=sub(Y,C);v=sub(T,C);z=sub(Z,C)
for a,b in ((u,v),(v,z),(u,z)):
 assert cross(a,b).lo>0 and dot(a,b).lo>0
out={'python':sys.version,'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'parent_script_sha256':hashlib.sha256(Path(__file__).with_name('certify_ah_path.py').read_bytes()).hexdigest(),'sigma_degrees_exact':'1.5494','visited_nodes':nodes,'leaf_count':len(leaves),'minimum_factor_lower':str(min(F(x['factor_lower']) for x in leaves)),'leaves':leaves,'passed':True,'scope':'fixed-angle R outside C disk plus tangent-path orientation; shortestness follows the analytic potential argument, not a new universal-cover theorem'}
print(json.dumps(out,indent=2))
