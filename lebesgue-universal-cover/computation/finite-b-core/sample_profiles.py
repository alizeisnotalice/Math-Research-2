"""Sample profile and mirror support points at every direction with Arb error checks."""
import argparse,json,sys,hashlib
from fractions import Fraction as F
from pathlib import Path
from flint import arb,ctx,__version__ as flint_version

SCALE=10**11
SOURCES=('tapered-forced-support-profiles.json','nonmonotone-curvature-profiles.json','tapered-double-deep-profile.json','ah-tapered-mirror-counterexample-input.json')
def hull(pts):
 def one(seq):
  out=[]
  for p in seq:
   while len(out)>1 and (out[-1][0]-out[-2][0])*(p[1]-out[-1][1])-(out[-1][1]-out[-2][1])*(p[0]-out[-1][0])<=0:out.pop()
   out.append(p)
  return out
 pts=sorted(set(pts))
 return pts if len(pts)<3 else one(pts)[:-1]+one(reversed(pts))[:-1]
def exact(q):return arb(q.numerator)/q.denominator

def main():
 ap=argparse.ArgumentParser();ap.add_argument('comp',type=Path);ap.add_argument('poses',type=Path);ap.add_argument('output',type=Path);ap.add_argument('--N',type=int,default=16384);args=ap.parse_args();assert args.N>=4;ctx.prec=160;comp=args.comp.resolve();sys.path.insert(0,str(comp))
 from certify_nonmonotone_ah_core import profile_jet,mirror_jet
 from certify_tapered_forced_support import Profile
 from certify_ah_tapered_mirror_counterexample import Support,Polygon
 rows=json.loads(args.poses.read_bytes())['poses'];selected=[r for r in rows if r['source'] in SOURCES];assert len(selected)==33 and len({r['index'] for r in selected})==33
 datasets={s:json.loads((comp/s).read_bytes()) for s in SOURCES};base=Support(datasets[SOURCES[-1]],arb.pi());mirror={'Support':base,'Polygon':Polygon(base)};normals=[]
 for k in range(args.N):
  sn,co=(2*arb.pi()*k/args.N).sin_cos();normals.append((co,sn))
 def point_mirror(body,q):
  z=(q-F(1,6))%2
  if z>1:z-=2
  sign=1 if z>=0 else -1;z=abs(z);width=1
  if z>F(1,2):z=1-z;width=-1
  theta=exact(z)*arb.pi();sn,co=theta.sin_cos();breaks=body.angles if isinstance(body,Support) else body.breaks;vals=[]
  for lo,hi,(a,b,d) in zip(breaks,breaks[1:],body.rows):
   if not(theta<lo or theta>hi):vals.append((arb(1)/2+width*(a*co+b*sn+d),sign*(b*co-a*sn)))
  assert vals,(q,z);h,v=vals[0]
  for a,b in vals[1:]:h=h.union(a);v=v.union(b)
  return h,v
 # The defining recurrences give identical value and derivative at each knot; assert numerically in Arb.
 for name,body in mirror.items():
  breaks=body.angles if isinstance(body,Support) else body.breaks
  for j,q in enumerate(breaks[1:-1]):
   sn,co=q.sin_cos();a,b,d=body.rows[j];c,e,f=body.rows[j+1];assert ((a-c)*co+(b-e)*sn+d-f).contains(0);assert ((b-e)*co-(a-c)*sn).contains(0)
  a,b,d=body.rows[0];assert b.contains(0);a,b,d=body.rows[-1];assert (b+d).contains(0)
 allgroups={'profiles':[],'mirror':[]};body_summary=[];maximum=arb(0);boundary=0
 for row in selected:
  name=row['source'];is_mirror=name==SOURCES[-1];obj=mirror[row['body_type']] if is_mirror else Profile(datasets[name][row['body_index']]);turn=F(row['turn'],180);sgn=row['reflection'];assert sgn in (-1,1)
  points=[];largest=arb(0)
  for k,(nx,ny) in enumerate(normals):
   q=(sgn*F(2*k,args.N)-turn)%2
   if is_mirror:h,dh=point_mirror(obj,q)
   else:
    y,v=obj.at(q.numerator,q.denominator);h=arb(1)/2+y;dh=v
    if (q*obj.grid).denominator==1:boundary+=1
   if k%127==1 and q not in (F(0),F(1),F(1,6),F(2,3),F(7,6),F(5,3)):
    py,px=(exact(q)*arb.pi()).sin_cos();legacy=(mirror_jet(obj,(px,py)) if is_mirror else profile_jet(obj,(px,py)))
    if legacy is not None:assert h.overlaps(legacy[0]) and dh.overlaps(legacy[1]),(row['index'],k,str(h),str(dh),legacy)
   x=h*nx-sgn*dh*ny;y=h*ny+sgn*dh*nx
   ix=round(float((x*SCALE).mid()));iy=round(float((y*SCALE).mid()))
   ex=abs(x-arb(ix)/SCALE);ey=abs(y-arb(iy)/SCALE);err=max(ex.upper(),ey.upper());assert err<arb('0.00000000002'),(row['index'],k,str(err));largest=max(largest,err);points.append((ix,iy))
  key='mirror' if is_mirror else 'profiles';allgroups[key]+=points;body_summary.append({'index':row['index'],'source':name,'body_index':row['body_index'],'turn':row['turn'],'reflection':sgn,'count':len(points),'max_coordinate_error_upper':str(largest),'integer_hull_vertices':len(hull(points))});maximum=max(maximum,largest);print('body',row['index'],'error',largest,flush=True)
 group_hulls={key:hull(val) for key,val in allgroups.items()};combined=hull(group_hulls['profiles']+group_hulls['mirror']);doc={'complete':True,'N':args.N,'scale':SCALE,'sample_count':len(selected)*args.N,'coordinate_error_upper':str(maximum),'coordinate_error_asserted':'1/50000000000','profile_exact_knot_samples':boundary,'pose_convention':json.loads(args.poses.read_bytes())['convention'],'input_sha256':{s:hashlib.sha256((comp/s).read_bytes()).hexdigest() for s in SOURCES},'dependency_sha256':{s:hashlib.sha256((comp/s).read_bytes()).hexdigest() for s in ('certify_tapered_forced_support.py','certify_ah_tapered_mirror_counterexample.py','certify_nonmonotone_ah_core.py')},'poses_sha256':hashlib.sha256(args.poses.read_bytes()).hexdigest(),'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'python':sys.version.split()[0],'flint':flint_version,'bodies':body_summary,'group_vertex_counts':{k:len(v) for k,v in group_hulls.items()},'combined_hull':combined};args.output.write_text(json.dumps(doc,separators=(',',':'))+'\n');print('COMPLETE bodies',len(selected),'samples',doc['sample_count'],'hull',len(combined),'maxerr',maximum,flush=True)
if __name__=='__main__':main()
