#!/usr/bin/env python3
"""Independent all-bbox single-root maximum certificate for the finite axis catalog."""
from fractions import Fraction as F
from pathlib import Path
from math import lcm,isqrt
import argparse,hashlib,json,os,sys,time
assert os.path.ismount('/Volumes/data')
t0=time.perf_counter();here=Path(__file__).resolve()
ap=argparse.ArgumentParser();ap.add_argument('--catalog',type=Path,default=here.parent/'inputs/reflection_axis_common_144.json');ap.add_argument('--seam',type=Path);args=ap.parse_args();catalog=args.catalog
blob=catalog.read_bytes();assert hashlib.sha256(blob).hexdigest()=='cbe911db74b373f11cb2bec5e47ec3effdad9a3dd7af1d8693ac2d6681d868cc'
z=json.loads(blob);pts=[tuple(map(F,t)) for t in z['pts']];assert len(pts)==144 and F(z['lam'])==1
xm=max(x for x,y,m in pts);target_m=max(m for x,y,m in pts if x==xm);q=target_m/4
xs=sorted({x for x,y,m in pts});ys=sorted({y for x,y,m in pts});nx,ny=len(xs),len(ys)
W=lcm(*(m.denominator for x,y,m in pts));xindex={x:i for i,x in enumerate(xs)};yindex={y:i for i,y in enumerate(ys)}
prefix=[[0]*(ny+1) for i in range(nx+1)]
for x,y,m in pts:prefix[xindex[x]+1][yindex[y]+1]+=int(m*W)
for i in range(1,nx+1):
    row=0
    for j in range(1,ny+1):row+=prefix[i][j];prefix[i][j]=prefix[i-1][j]+row
xpairs=[(i,k,a,b,(b-a)**2) for i,a in enumerate(xs) for k,b in enumerate(xs[i:],i)]
ypairs=[(j,l,c,d,(d-c)**2) for j,c in enumerate(ys) for l,d in enumerate(ys[j:],j)]
counts=dict(boxes=0,zero_mass=0,empty_receivers=0,nonempty_receivers=0,p_le_q=0,L_le_zero=0,squared_comparison=0,strict_upper=0,equality=0)
ties=[];min_square_margin=None;min_mass_slack=None
for i,k,a,b,xspan2 in xpairs:
  for j,l,c,d,yspan2 in ypairs:
    counts['boxes']+=1
    M=prefix[k+1][l+1]-prefix[i][l+1]-prefix[k+1][j]+prefix[i][j]
    if not M:counts['zero_mass']+=1;continue
    span2=max(xspan2,yspan2)
    if span2.numerator*W>=M*span2.denominator:counts['empty_receivers']+=1;continue
    counts['nonempty_receivers']+=1
    p=F(M,4*W);gap=xm-a;assert gap>=0
    mass_slack=q+gap*gap-p;assert mass_slack>=0
    if mass_slack>0:min_mass_slack=mass_slack if min_mass_slack is None else min(min_mass_slack,mass_slack)
    equal=False
    if p<=q:
        counts['p_le_q']+=1;equal=(p==q and gap==0)
    else:
        L=p-q-gap*gap
        if L<=0:counts['L_le_zero']+=1
        else:
            counts['squared_comparison']+=1
            margin=4*gap*gap*q-L*L
            assert margin>=0,(a,b,c,d,p,q,gap,L)
            equal=(margin==0)
            if not equal:min_square_margin=margin if min_square_margin is None else min(min_square_margin,margin)
    if equal:
        counts['equality']+=1;ties.append(dict(bbox=list(map(str,(a,b,c,d))),mass=str(F(M,W))))
    else:counts['strict_upper']+=1
assert counts['boxes']==nx*(nx+1)//2*ny*(ny+1)//2
assert counts['nonempty_receivers']==5066 and counts['equality']==4
assert all(t['bbox'][0]==t['bbox'][1]==str(xm) and t['bbox'][2]==t['bbox'][3] and F(t['mass'])==target_m for t in ties)
assert isqrt(q.numerator)**2!=q.numerator or isqrt(q.denominator)**2!=q.denominator
ops=args.seam;opsblob=None;lo=hi=None
if ops is not None:
    opsblob=ops.read_bytes();o=json.loads(opsblob)
    assert F(o['max_lower_witness']['mass'])==target_m and F(o['max_lower_witness']['bbox'][0])==xm
    lo,hi=F(o['c_edge_lower']),F(o['c_edge_upper'])
    assert 0<lo-xm<hi-xm and (lo-xm)**2<=q<(hi-xm)**2
out=dict(status='PASS',python=sys.version,arithmetic='Fraction; exact rational squaring of a+sqrt(p)<=xmax+sqrt(q)',random_seed=None,floating_tolerance=None,catalog_sha256=hashlib.sha256(blob).hexdigest(),catalog_atoms=len(pts),catalog_x=nx,catalog_y=ny,xmax=str(xm),target_singleton_mass=str(target_m),target_radius_squared=str(q),exact_C_formula='xmax + sqrt(target_radius_squared)',C_is_irrational=True,counts=counts,equality_witnesses=ties,minimum_positive_squared_margin=None if min_square_margin is None else str(min_square_margin),all_nonempty_bbox_mass_bound="M_B <= target_singleton_mass + 4*(xmax-a)^2",minimum_positive_mass_bound_slack=str(min_mass_slack),operations_target_source=str(ops) if ops else None,operations_target_sha256=hashlib.sha256(opsblob).hexdigest() if opsblob else None,operations_interval_independently_contains_exact_C=True if ops else None,c_edge_lower=str(lo) if lo is not None else None,c_edge_upper=str(hi) if hi is not None else None,script_sha256=hashlib.sha256(here.read_bytes()).hexdigest(),elapsed_seconds=time.perf_counter()-t0,scope='Complete 63000 coordinate bboxes of the already certified common 144-atom owner catalog; no dyadic root pruning and no author seam implementation reused')
path=here.parent/'results/reflection-axis-single-root/exact-certificate.json';path.parent.mkdir(parents=True,exist_ok=True);path.write_text(json.dumps(out,ensure_ascii=False,indent=2)+'\n')
print(json.dumps({k:out[k] for k in ['status','catalog_x','catalog_y','xmax','target_singleton_mass','counts','elapsed_seconds']},ensure_ascii=False,indent=2))
