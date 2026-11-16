"""One exact full-grid fixture separating origin-triple and joint-prefix activity."""
from fractions import Fraction as F
from pathlib import Path
from itertools import combinations_with_replacement
import hashlib,json,os,platform,sys,time
start=time.monotonic();n=5;e=F(1,25);h=F(6,25)
old=[(i*h,j*h)for i in range(5)for j in range(5)]
outside=[(-F(4,25),F(1,5)),(F(1,5),-F(4,25))]
points=old+outside;masses=[e]*25+[F(17,200)]*2
axes=[sorted({p[j]for p in points})for j in (0,1)]
rows={}
for a,b in combinations_with_replacement(axes[0],2):
    for c,d in combinations_with_replacement(axes[1],2):
        mask=sum(1<<i for i,(x,y)in enumerate(points)if a<=x<=b and c<=y<=d)
        if not mask or mask in rows:continue
        ids=[i for i in range(len(points))if mask>>i&1]
        span=max(max(points[i][j]for i in ids)-min(points[i][j]for i in ids)for j in (0,1))
        mass=sum((masses[i]for i in ids),F(0))
        slack=max(F(1),span*span)-mass;assert slack>=0
        rows[mask]=(mass,span)
eligible={mask:v for mask,v in rows.items()if v[0]>v[1]*v[1]}
L=max(m for m,d in eligible.values());winners=[mask for mask,(m,d)in eligible.items()if m==L]
assert L==1 and winners==[(1<<25)-1]
component=set(range(25));changed=True
while changed:
    changed=False
    for mask in eligible:
        ids={i for i in range(27)if mask>>i&1}
        if ids&component and not ids<=component:component|=ids;changed=True
assert len(component)==27
span=max(axes[j][-1]-axes[j][0]for j in (0,1));total=sum(masses,F(0))
nonterminal=span*span-F(35,33)*total;assert nonterminal==F(371,27500)>0
def label(ids):
    aa=min(points[i][0]for i in ids);bb=max(points[i][0]for i in ids)
    cc=min(points[i][1]for i in ids);dd=max(points[i][1]for i in ids)
    delta=max((bb-aa)/2,(dd-cc)/2,bb+e,dd+e,-F(1,2)-aa,-F(1,2)-cc,F(0))
    mass=sum((masses[i]for i in ids),F(0));gap=mass-4*delta*delta
    return {'indices':ids,'bbox':list(map(str,[aa,bb,cc,dd])),'mass':str(mass),'delta_Gamma':str(delta),'activity_gap':str(gap),'active':gap>0}
origin=label([0,25,26]);prefix=label([0,1,5,6,25,26])
assert origin['active']is False and F(origin['activity_gap'])==-F(51,2500)
assert prefix['active']is True and F(prefix['activity_gap'])==F(41,2500)
x=(-F(41,1000),-F(41,1000));r=F(57,200)
assert all(-F(1,2)<v<-e for v in x)
captured=[i for i,p in enumerate(points)if max(abs(p[j]-x[j])for j in (0,1))<r]
assert captured==prefix['indices']
fullgap=sum((masses[i]for i in captured),F(0))-4*r*r;assert fullgap==F(51,10000)>0
residual=[]
for mask in range(1,4):
    ids=[25+j for j in range(2)if mask>>j&1]
    distance=max(abs(points[i][j]-x[j])for i in ids for j in (0,1))
    mass=sum((masses[i]for i in ids),F(0));gap=4*distance*distance-mass;assert gap>0
    residual.append({'indices':ids,'minimum_witness_radius_at_x':str(distance),'nonhalo_gap':str(gap)})
out={'status':'PASS','n':5,'positions':[[str(v)for v in p]for p in points],'masses':list(map(str,masses)),
 'literal_bbox_closures':len(rows),'eligible_closures':len(eligible),'all_square_caps':True,'unique_maximum_mask':winners[0],
 'maximum_eligible_mass':'1','connected':True,'nonterminal_gap':str(nonterminal),
 'origin_triple':origin,'actual_2x2_prefix_joint_label':prefix,'strict_full_witness':{'center':list(map(str,x)),'radius':str(r),'full_margin':str(fullgap)},
 'complete_residual_label_exclusions_at_center':residual,
 'scope':'One directed finite admissible unique-maximum/connected/nonterminal fixture: origin triple is null but a larger old-prefix joint label is active. No area/payment violation claimed; the center is already in the old singleton halo.',
 'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
 'runtime':{'python':platform.python_version(),'pid':os.getpid(),'seconds':time.monotonic()-start}}
Path(sys.argv[1]).write_text(json.dumps(out,indent=2)+'\n');print(json.dumps({k:v for k,v in out.items()if k not in ('positions','masses')}))
