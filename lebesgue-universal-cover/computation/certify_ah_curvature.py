"""Certified positive turning of the AH replacement boundary at fixed angle."""
from pathlib import Path
from fractions import Fraction as F
import contextlib,io,runpy,json,sys,hashlib
base=Path(__file__).parent
saved=sys.argv;sys.argv=[str(base/'certify_ah_area.py'),'64']
with contextlib.redirect_stdout(io.StringIO()):d=runpy.run_path(str(base/'certify_ah_area.py'))
sys.argv=saved
J=d['J'];I=d['I'];Z=d['Z'];T=d['T'];bd=d['d'];C=bd['C'];G=bd['G'];E=bd['E'];Y=bd['Y']
def at(w):
 q=J(w);q.c[1]=I(1)
 return d['paired'](d['add'](Z,d['mul'](q,d['sub'](T,Z))))
def cross(a,b):return a[0]*b[1]-a[1]*b[0]
stack=[(F(0),F(1),0)];leaves=[]
while stack:
 a,b,depth=stack.pop()
 try:
  p=at(I(a,b));v=(p[0].c[1],p[1].c[1]);acc=(2*p[0].c[2],2*p[1].c[2]);k=cross(v,acc)
  certified=k.lo>0 and v[0].hi<0 and v[1].hi<0
 except AssertionError:
  certified=False  # Coarse interval may fail a radicand enclosure; keep both children.
 if certified:
  leaves.append({'interval':[str(a),str(b)],'curvature_numerator_lower':str(k.lo),'dx_upper':str(v[0].hi),'dy_upper':str(v[1].hi)});continue
 assert depth<16,'unresolved interval'
 m=(a+b)/2;stack.extend([(a,m,depth+1),(m,b,depth+1)])
intervals=sorted((F(r['interval'][0]),F(r['interval'][1])) for r in leaves)
assert intervals[0][0]==0 and intervals[-1][1]==1
assert all(b==c for (a,b),(c,e) in zip(intervals,intervals[1:]))
# Convex corner at Y: incoming E circle, outgoing C circle.
incoming=bd['perp'](bd['sub'](Y,E));outgoing=bd['perp'](bd['sub'](Y,C));corner=cross(incoming,outgoing)
assert corner.lo>0 and bd['dot'](incoming,outgoing).lo>0
end=at(I(1));Tp=(end[0].c[0],end[1].c[0]);velocity=(end[0].c[1],end[1].c[1])
orientation=bd['dot'](velocity,bd['perp'](bd['sub'](Tp,G)));assert orientation.lo>0
print(json.dumps({'python':sys.version,'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'sigma_degrees_exact':'1.5494','leaf_count':len(leaves),'minimum_curvature_numerator_lower':str(min(F(x['curvature_numerator_lower']) for x in leaves)),'Y_corner_cross_lower':str(corner.lo),'paired_end_forward_tangent_dot_lower':str(orientation.lo),'leaves':leaves,'passed':True,'scope':'local AH replacement boundary has convex turning; this does not prove all remaining boundary pieces or universal coverage'},indent=2))
