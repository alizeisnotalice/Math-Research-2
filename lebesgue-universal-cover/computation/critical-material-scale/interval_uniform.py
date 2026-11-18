"""Certify every scale-sensitive condition uniformly on a closed rational interval."""
from pathlib import Path
from fractions import Fraction as F
from functools import reduce
from itertools import combinations
from hashlib import sha256
import sys,json,tarfile
import numpy as np
from flint import arb,ctx
ctx.prec=256
D=Path(__file__).resolve().parent;comp=Path(sys.argv[1]).resolve();treepath=Path(sys.argv[2]).resolve();out=Path(sys.argv[3]);sys.path[:0]=[str(D),str(comp)]
from certify_polygon_083731 import ar,lower_radius
from certify_tomar_six_high import SourceTree
from select_tomar_083731_scores import radius_rows
from tomar_d3_mask import CUT
from cap_final import sqrt_upper
from displacements_final import POSITIVES
L=F('.99999900555134042347742929191955');U=F('.999999005551341');S=ar(L).union(ar(U));G=F(120689264254125501,2**57);R=F('.703');r=F(1,2);wall=F('.4742688878');M=F('.8372950062210252434563');m=F(57663128211,124250000000);eps=1-S;sr=S/2;pi=arb.pi();read=lambda n:json.loads((comp/n).read_bytes());tree=SourceTree(treepath)
with tarfile.open(comp/'six-input-083745/volumes/native45-receipts.tar.gz') as t: scope=json.load(t.extractfile('inputs/ordinary/scope.json'))
keys=[tuple(x) for x in scope['leaves']+scope['deferred']];assert len(keys)==len(set(keys))==3273 and len(scope['leaves'])==3071
mx=max(sum(max(abs(a),abs(b))**2 for a,b in tree.leaf(*key)[:2]) for key in keys);assert mx==F(14275988056445,18014398509481984) and mx<F('.02816')**2
seed=read('tomar-083745-future-seed.json');radii=radius_rows(read('five-body-radial-affine.json'),F('.83745'));assert seed['radii']==radii and seed['target_integer']==G*2**57
root=tree.root;assert root['ns']==[3,5,7,9] and root['split0']==[2]*11 and all(F(root['lo'][j])==0 for j in (2,5,8));assert pi/15<ar(CUT)<ar(F(root['hi'][2])) and all(ar(F(root['hi'][j]))>2*pi/n for j,n in ((5,7),(8,9)))
exterior=[]
for i,n in enumerate((3,5,7,9)):
 j=3*i;d=min(F(root['hi'][j]),F(root['hi'][j+1]),-F(root['lo'][j]),-F(root['lo'][j+1]));inner=S*lower_radius(n,pi);margin=(pi*(sr**2+inner**2)/2+ar(d)*(sr+inner)-ar(G)).lower();assert margin>0;assert F(radii[i]['integer_area_lower'],2**57)>G and F(radii[i]['integer_radius'],2**44)<min(d,F('.069'));exterior.append(str(margin))
radial=read('five-body-radial-affine.json');packet=comp/'five-body-radial-material.npz';assert sha256(packet.read_bytes()).hexdigest()==radial['table_sha']
with np.load(packet,allow_pickle=False) as z:assert list(z['counts'][:2])==[96,96];pts=z['points'][96:192]
inner=[(ar(F(int(x),2**28)),ar(F(int(y),2**28))) for x,y in pts];radius=1/(2*(pi/6).cos());corners=[]
for j in range(3):sn,cs=(2*pi*j/3).sin_cos();corners.append((radius*cs,radius*sn))
assert all((x-vx)**2+(y-vy)**2<1 for x,y in inner for vx,vy in corners)
cloud=read('finite-point-leaf-input.json');dirs=cloud['directions'];assert len(dirs)==105;direct=clip=0;min_clip=None
for x,y in dirs:
 length=arb(x*x+y*y).sqrt();nx,ny=arb(x)/length,arb(y)/length;low=reduce(arb.max,(nx*px+ny*py for px,py in inner));high=1+reduce(arb.min,(nx*vx+ny*vy for vx,vy in corners))
 if low>=ar(wall):direct+=1
 else:
  gap=(ar(r)-eps*ar(wall)-ar(F('.02816'))-S*high).lower();assert gap>0;min_clip=gap if min_clip is None else min(min_clip,gap);clip+=1
assert direct==64 and clip==41 and all(lower_radius(n,pi)>ar(wall) for n in (5,7,9))
disk=((pi-(sr/ar(R)).acos())*sr**2+sr*(ar(R)**2-sr**2).sqrt()-ar(G)).lower();assert disk>0 and ar(R*R)<2*sr**2
delta=ar(F('.6'))*eps;pure=(ar(F('.837455')-G)-ar(F(22,7))*(2*ar(R)*delta+delta**2)).lower();assert pure>0
assert read('five-global-minimum-single-worker-replay.json')['complete'] and arb(read('five-global-minimum-single-worker-replay.json')['minimum_interval'])>ar(M)
# Worst-endpoint rational cap bound, uniform in the scale interval.
vl=(1+U)*r*r+(1-U)*m*m;budget=4*(wall-m)*L*L;deficit=G-L*L*M;assert deficit>0
covered=[]
def cover(lo,hi,depth):
 assert depth<=60
 bound=(1+U)*(hi-r)+max(F(0),(vl-2*m*hi)/(hi+r))
 if lo>r:
  height=min(hi-r,sqrt_upper(deficit*sqrt_upper(hi*hi-(L*r)**2)/(L*r)))
  correction=max(F(0),(vl-2*m*lo)/(lo*lo-r*r));bound=min(bound,height*(1+U+correction))
 if bound<budget:covered.append((lo,hi,budget-bound))
 else:
  mid=(lo+hi)/2;cover(lo,mid,depth+1);cover(mid,hi,depth+1)
end=R+m*(1-L);cover(r,end,0);assert covered[0][0]==r and covered[-1][1]==end and all(x[1]==y[0] for x,y in zip(covered,covered[1:]))
# Full deformed-point diameter under any epsilon in [1-U,1-L].
P=cloud['points'];den=cloud['denominator'];Z=[]
for x,y in P:
 j=P.index([x,abs(y)]);zx,zy=POSITIVES[j];Z.append((zx,zy if y>=0 else -zy))
assert set(POSITIVES)=={j for j,(x,y) in enumerate(P) if y>=0};assert all((a,0)==(a,b) for (_,y),(a,b) in zip(P,Z) if y==0)
gaps=[]
for i,j in combinations(range(30),2):
 dif=[]
 for k in range(2):
  a=F(P[i][k]-P[j][k],den);b=F(Z[i][k]-Z[j][k],10);dif.append(max(abs(a-(1-L)*b),abs(a-(1-U)*b)))
 gap=L*L-sum(x*x for x in dif);assert gap>0;gaps.append((gap,i,j))
assert len(gaps)==435
proof_path=Path(sys.argv[4]) if len(sys.argv)>4 else D/'displacements-final/proof.json'
proof=json.loads(proof_path.read_text());assert proof['complete'] and proof['angular_rows']==2493 and proof['total_angle_cells']==384 and F(proof['old_rate'])==m and arb(proof['weighted_max_upper'])<ar(F('.419')) and F('.419')<m
res={'lower':str(L),'upper':str(U),'native_leaves':len(keys),'root_exterior_margins':exterior,'directions':[direct,clip],'minimum_clipped_margin':str(min_clip),'disk_margin':str(disk),'pure_margin':str(pure),'cap_slabs':len(covered),'cap_minimum_margin':str(min(x[2] for x in covered)),'cap_end':str(end),'point_pairs':len(gaps),'point_minimum_margin':str(min(gaps)[0]),'point_minimum_pair':list(min(gaps)[1:]),'weighted_rows':2493,'weighted_bound':proof['weighted_max_upper'],'complete':True}
out.write_text(json.dumps(res,indent=2)+'\n');print(json.dumps(res,indent=2))
