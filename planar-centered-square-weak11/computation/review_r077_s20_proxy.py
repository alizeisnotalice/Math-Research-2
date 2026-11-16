"""One independent raw-S20 capacity audit and exact synthetic-proxy area."""
from pathlib import Path
from fractions import Fraction as F
from functools import cmp_to_key
from collections import defaultdict,Counter
import hashlib,json,os,platform,sys,time
import sympy as S
start=time.monotonic();cp=Path(sys.argv[1]);ap=Path(sys.argv[2]);cr=cp.read_bytes();ar=ap.read_bytes()
case=next(c for c in json.loads(cr)['cases']if c['id']=='n20-corrected');actual=json.loads(ar)
n=20;scale=n*n;e=F(1,scale);h=F(n+1,scale);p=F(1,2)-e
points=[tuple(map(F,z))for z in actual['positions']];masses=list(map(F,actual['masses']))
old=[(i*h,j*h)for i in range(n)for j in range(n)]
proxy_points=[tuple(map(F,z))for z in case['positions']];proxy_mass=list(map(F,case['masses']))
assert points[:400]==old and masses[:400]==[e]*400 and actual['S_indices']==list(range(400))
assert points[400:]==proxy_points[1:]==[(-F(71,200),F(0)),(F(0),-F(171,200))]
assert masses[400:]==proxy_mass[1:]==[F(1,4),F(1,2)] and F(actual['threshold'])==1
local_old=[i for i,z in enumerate(old)if all(-1<v<p for v in z)]
mP=len(local_old)*e;assert mP==proxy_mass[0]==F(1,4) and len(local_old)==100
P=[tuple(int(v*scale)for v in z)for z in points];assert all(F(v,scale)==q for z,qq in zip(P,points)for v,q in zip(z,qq))
units=[int(m*scale)for m in masses];assert list(map(F,units))==[m*scale for m in masses]
groups=[];axis_masks=[]
for axis in (0,1):
    g=defaultdict(int)
    for i,z in enumerate(P):g[z[axis]]|=1<<i
    g=sorted(g.items());groups.append(g);rr=[]
    for j in range(len(g)):
        mask=0
        for coord,part in g[j:]:mask|=part;rr.append(mask)
    axis_masks.append(rr)
closures={a&b for a in axis_masks[0]for b in axis_masks[1]}-{0}
oldmask=(1<<400)-1;D=scale*scale;class_min={1:None,2:None,3:None};tags=Counter();eligible=0;maxnum=0;winners=[];jointmargin=None
for mask in closures:
    ranges=[[v for v,g in gg if mask&g]for gg in groups]
    span=max(v[-1]-v[0]for v in ranges);tag=mask>>400;count=(mask&oldmask).bit_count()
    massnum=count+units[400]*bool(tag&1)+units[401]*bool(tag&2)
    cap=max(D,span*span);assert massnum*scale<=cap
    tags[tag]+=1
    if tag:
        rhs=F(cap-count*scale,D)
        if class_min[tag]is None or rhs<class_min[tag]:class_min[tag]=rhs
    if massnum*scale>span*span:
        eligible+=1
        if massnum>maxnum:maxnum=massnum;winners=[mask]
        elif massnum==maxnum:winners.append(mask)
    if mask==(1|(1<<400)|(1<<401)):
        jointmargin=F(massnum*scale-span*span,D)
assert maxnum==400 and winners==[oldmask] and jointmargin>0
# Full old block is eligible; this eligible joint label joins both outside
# atoms to its origin, proving a single eligible connected component.
span=max(max(z[j]for z in points)-min(z[j]for z in points)for j in (0,1))
nonterminal=span*span-F(35,33)*sum(masses,F(0));assert nonterminal>0
a1=-points[400][0];a2=-points[401][1];b1,b2=masses[400:]
v1=[max(F(1),(a1+(k-1)*h)**2)-F(k,n)for k in range(1,n+1)]
v2=[max(F(1),(a2+(l-1)*h)**2)-F(l,n)for l in range(1,n+1)]
u=[max(F(1),max(a1+(k-1)*h,a2+(l-1)*h)**2)-F(k*l,n*n)for k in range(1,n+1)for l in range(1,n+1)]
assert all(b1<=v for v in v1)and all(b2<=v for v in v2)and all(b1+b2<=v for v in u)
assert (min(v1),min(v2),min(u))==tuple(class_min[t]for t in (1,2,3))
activity=[]
for k in range(1,11):
    for l in range(1,11):
        loX=-a1;hiX=(k-1)*h;loY=-a2;hiY=(l-1)*h
        delta=max((hiX-loX)/2,(hiY-loY)/2,hiX+e,hiY+e,-F(1,2)-loX,-F(1,2)-loY,F(0))
        if F(k*l,400)+b1+b2>4*delta*delta:activity.append([k,l])
Q=(-a1,1-a1,-a2,1-a2)
inside=lambda z:Q[0]<=z[0]<=Q[1]and Q[2]<=z[1]<=Q[3]
captured=[i for i,z in enumerate(points)if inside(z)];original=sum((masses[i]for i in captured),F(0))
collapsed=original-sum((masses[i]for i in local_old if inside(points[i])),F(0))+mP
assert original==F(339,400)and collapsed==F(409,400)>1

# Separate synthetic mixed-only proxy: no 400-source halo is built.
pp=[tuple(S.Rational(v.numerator,v.denominator)for v in z)for z in proxy_points]
pm=[S.Rational(v.numerator,v.denominator)for v in proxy_mass]
def cmp(a,b):
    z=S.simplify(a-b)
    return 0 if z==0 else(-1 if bool(z<0)else 1)
def lt(a,b):return cmp(a,b)<0
def mn(a,b):return a if cmp(a,b)<=0 else b
def mx(a,b):return b if lt(a,b)else a
def merge(iv):
    result=[]
    for a,b in sorted(iv,key=cmp_to_key(lambda i,j:cmp(i[0],j[0]))):
        if not result or lt(result[-1][1],a):result.append([a,b])
        else:result[-1][1]=mx(result[-1][1],b)
    return result
ee=S.Rational(1,400);full=[];res=[];catalog=[]
for mask in (2,3,4,5,6,7):
    ids=[i for i in range(3)if mask>>i&1];mass=sum(pm[i]for i in ids);r=S.sqrt(mass)/2
    if mask&1:r=mn(r,S.Rational(1,2))
    rr=(mx(max(pp[i][0]for i in ids)-r,-S.Rational(1,2)),mn(min(pp[i][0]for i in ids)+r,-ee),mx(max(pp[i][1]for i in ids)-r,-S.Rational(1,2)),mn(min(pp[i][1]for i in ids)+r,-ee))
    nonempty=lt(rr[0],rr[1])and lt(rr[2],rr[3]);catalog.append({'mask':mask,'radius':str(r),'clipped':list(map(str,rr))if nonempty else None})
    if nonempty:(full if mask&1 else res).append(tuple(map(S.simplify,rr)))
yy=sorted(set(v for rr in full+res for v in rr[2:]),key=cmp_to_key(cmp));area=S.Integer(0);slabs=[]
for y0,y1 in zip(yy,yy[1:]):
    y=(y0+y1)/2;ff=merge([(a,b)for a,b,c,d in full if lt(c,y)and lt(y,d)]);rr=merge([(a,b)for a,b,c,d in res if lt(c,y)and lt(y,d)])
    length=sum((b-a for a,b in ff),S.Integer(0))
    for a,b in ff:
        for c,d in rr:
            lo=mx(a,c);hi=mn(b,d)
            if lt(lo,hi):length-=hi-lo
    part=S.simplify((y1-y0)*length);assert cmp(part,0)>=0;area+=part;slabs.append({'y0':str(y0),'y1':str(y1),'mixed_only_area':str(part)})
area=S.simplify(area);formula=43*S.sqrt(3)/200-101*S.sqrt(2)/1600-S.Rational(24641,160000)
assert S.simplify(area-formula)==0
assert F(433,250)**2<3 and F(283,200)**2>2
lower=F(43,200)*F(433,250)-F(101,1600)*F(283,200)-F(24641,160000)
assert lower==F(206483,1600000)>F(1,8)
out={'status':'PASS','original_full_grid':{'n':20,'old_atoms':400,'total_atoms':402,'grid_spacing':str(h),'internal_atom_mass':str(e),'actual_bbox_closures':len(closures),'closure_counts_by_outside_mask':dict(tags),'eligible_closures':eligible,'maximum_eligible_mass':'1','unique_maximum_is_full_S20':True,'connected':True,'outside_joint_label_connecting_margin':str(jointmargin),'nonterminal_gap':str(nonterminal),'raw_full_bbox_bounds_V1_V2_U':{str(t):str(v)for t,v in class_min.items()},'independent_axis_prefix_rows':440,'local_internal_atom_count':100,'local_internal_mass':'1/4'},
 'actual_old_prefix_joint_activity':{'prefix_count':100,'active_count':len(activity),'active_prefixes':activity,'scope':'activity only, no area'},
 'synthetic_proxy':{'mixed_masks':[3,5,7],'complete_residual_masks':[2,4,6],'old_singleton_added':False,'full_radius_cutoff':'1/2','Gamma_cut':'1/400','receiver_catalog':catalog,'horizontal_slabs':slabs,'exact_area':str(area),'decimal_display':float(area),'simple_strict_lower':str(lower),'target':'1/8','lower_excess':str(lower-F(1,8))},
 'collapse_while_retaining_other_old_atoms':{'unit_square':list(map(str,Q)),'original_mass':str(original),'synthetic_after_collapse_mass':str(collapsed),'violates_cap':True},
 'scope':'Original full S20 input is qualified; synthetic mixed-only clipped/truncated proxy exceeds1/8. No actual S20 halo/Delta area, C2 or main corner counterexample is claimed. Mixed-only omits the separate old-singleton union term but does not subtract mixed overlap with old halo.',
 'proxy_input_sha256':hashlib.sha256(cr).hexdigest(),'actual402_input_sha256':hashlib.sha256(ar).hexdigest(),'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
 'runtime':{'python':platform.python_version(),'sympy':S.__version__,'pid':os.getpid(),'seconds':time.monotonic()-start}}
Path(sys.argv[3]).write_text(json.dumps(out,indent=2)+'\n');print(json.dumps({k:v for k,v in out.items()if k not in ('synthetic_proxy',)}|{'proxy_exact':str(area),'proxy_decimal':float(area),'proxy_lower':str(lower)}))
