"""Independent partition/inequality replay for the fixed-angle shortest path."""
from pathlib import Path
from fractions import Fraction as F
import runpy,io,contextlib,json,sys
base=Path(__file__).parent
with contextlib.redirect_stdout(io.StringIO()):d=runpy.run_path(str(base/'certify_ah_path.py'))
r=json.loads((base/'019-ah-shortest-certificate.json').read_text());leaves=r['leaves']
I=d['I'];add=d['add'];sub=d['sub'];mul=d['mul'];dot=d['dot'];upper=d['upper_cross']
C,G,E,Nm=(d[k] for k in ('C','G','E','Nm'));lg=sub(G,C);ng=sub(E,Nm);eL=d['dL'];length=dot(lg,eL)
rects=[];total=F(0)
for leaf in leaves:
 a,b=map(F,leaf['u']);c,e=map(F,leaf['v']);assert 0<=a<b<=1 and 0<=c<e<=1
 total+=(b-a)*(e-c);rects.append((a,b,c,e))
 u=I(a,b);L=add(C,mul(u,lg));N=add(Nm,mul(I(c,e),ng));X=upper(L,N)
 factor=2*dot(sub(X,L),eL)+u*length
 assert factor.lo>=F(leaf['factor_lower'])>0
assert total==1
for i,(a,b,c,e) in enumerate(rects):
 for f,g,h,j in rects[:i]:assert min(b,g)<=max(a,f) or min(e,j)<=max(c,h)
def cross(a,b):return a[0]*b[1]-a[1]*b[0]
u=sub(d['Y'],C);v=sub(d['T'],C);z=sub(d['Z'],C)
for a,b in ((u,v),(v,z),(u,z)):assert cross(a,b).lo>0 and dot(a,b).lo>0
print(json.dumps({'python':sys.version,'partition_rectangles':len(rects),'total_area_exact':str(total),'overlap_pairs_checked':len(rects)*(len(rects)-1)//2,'all_leaf_inequalities_recomputed':True,'path_and_orientation_recomputed':True,'passed':True},indent=2))
