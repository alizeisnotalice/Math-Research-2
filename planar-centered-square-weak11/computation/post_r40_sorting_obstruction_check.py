#!/usr/bin/env python3
from fractions import Fraction as F
from itertools import combinations,product
from pathlib import Path
import hashlib,json,platform,sys,time
ROOT=Path(__file__).resolve().parents[1]
start=time.perf_counter();a=F(99,100);eps=F(1,10000)
points=tuple((x,y) for x,y in product((-a,F(0),a),repeat=2) if (x,y)!=(0,0))
M=1+8*eps;delta=1-a
assert a*a<1+eps and 8*eps<a*a and 1+8*eps<4*a*a
rays=((1,0),(1,1),(0,1),(-1,1),(-1,0),(-1,-1),(0,-1),(1,-1))
sectors=((2,1),(1,2),(-1,2),(-2,1),(-2,-1),(-1,-2),(1,-2),(2,-1))
counts=[]
for u,v in rays+sectors:
    positive=sum(u*x+v*y>0 for x,y in points)
    zero=sum(u*x+v*y==0 for x,y in points)
    assert positive>=3 and zero in (0,2)
    counts.append({'functional':[u,v],'strict_positive_light_points':positive,'zero_light_points':zero})
for p,q in combinations(((F(0),F(0)),)+points,2):
    overlap=max(F(0),1-abs(p[0]-q[0]))*max(F(0),1-abs(p[1]-q[1]))
    assert overlap<=delta
quad_count=0
for group in combinations(points,4):
    span=max(max(p[0] for p in group)-min(p[0] for p in group),max(p[1] for p in group)-min(p[1] for p in group))
    assert span==2*a;quad_count+=1
assert quad_count==70
assert (2*a-F(1,100))**2>1+3*eps
projection_lower=4-6*delta;free_exact=1+11*eps
assert projection_lower==F(197,50)
assert projection_lower/M==F(4925,1251)
assert free_exact==F(10011,10000) and free_exact/M==F(3337,3336)
assert 2-delta>free_exact
# Corner cluster: p mass1 at0, q mass1 at1, r mass14 at4.
assert F(1)<1+1
assert F(4)**2>=14+1
assert F(4)**2==14+1+1
assert F(3)**2<14+1
assert F(9)**2>32  # 9>4sqrt2, hence sqrt15>2+sqrt2
assert F(15)-2==13>2
# 100-center cluster, individual-mass ascending failure and good free order.
center_diameter=F(99,100000);perimeter_mass=F(1,50);cluster_total=F(1)
assert F(1,100)<perimeter_mass
assert center_diameter**2<1
assert (a+center_diameter)**2<1+perimeter_mass
assert center_diameter<F(1,100)
M_cluster=1+8*perimeter_mass
mass_sort_lower=(2*a+1)**2
assert M_cluster==F(29,25) and mass_sort_lower==F(22201,2500)
assert mass_sort_lower/M_cluster==F(22201,2900)
assert 8*perimeter_mass<a*a and 4*a*a>M_cluster
assert F(53,50)<F(103,100)**2
free_cluster_upper=F(4,25)+F(53,50)+F(1,1000)*F(103,100)
assert free_cluster_upper==F(122103,100000)
# Six-cycle: no deletion leaves a component of size <=4.
h=F(49,100)
hexagon=((h,F(0)),(h/2,h),(-h/2,h),(-h,F(0)),(-h/2,-h),(h/2,-h))
adj={j:set() for j in range(6)}
for j,k in combinations(range(6),2):
    span=max(abs(hexagon[j][0]-hexagon[k][0]),abs(hexagon[j][1]-hexagon[k][1]))
    if span*span<F(1,3): adj[j].add(k);adj[k].add(j)
assert {tuple(sorted((j,k))) for j in adj for k in adj[j]}=={tuple(sorted((j,(j+1)%6))) for j in range(6)}
deletions=[]
for removed in range(6):
    vertices=set(range(6))-{removed};seen={min(vertices)}
    while True:
        grown=seen|{k for j in seen for k in adj[j] if k in vertices}
        if grown==seen:break
        seen=grown
    assert seen==vertices
    pts=[hexagon[j] for j in vertices]
    span=max(max(p[0] for p in pts)-min(p[0] for p in pts),max(p[1] for p in pts)-min(p[1] for p in pts))
    assert span==2*h and span*span>F(5,6)
    deletions.append({'removed':removed,'connected_old_component_size':len(seen),'whole_component_eligible':False,'bbox_span_squared':str(span*span)})
assert (2*h)**2<1
proper_rank_checks=0
for size in (3,4,5):
    for subset in combinations(hexagon,size):
        span=max(max(p[0] for p in subset)-min(p[0] for p in subset),max(p[1] for p in subset)-min(p[1] for p in subset))
        assert span*span>=F(size,6)
        proper_rank_checks+=1
assert proper_rank_checks==41
proof=ROOT/'proofs/post-r40-sorting-obstructions.md'
output={'status':'PASS','python':sys.version,'platform':platform.platform(),
        'arithmetic':'fractions.Fraction exact; no random inputs or floating decisions',
        'heavy_mass':'1','light_mass':str(eps),'a':str(a),'total_mass':str(M),
        'projection_area_lower':str(projection_lower),'projection_ratio_lower':str(projection_lower/M),
        'free_proxy_area_exact':str(free_exact),'free_proxy_ratio_exact':str(free_exact/M),
        'projection_sign_cases':counts,'four_light_bbox_checks':quad_count,
        'mass_sort_cluster':{'center_count':100,'center_mass':'1/100','perimeter_mass':'1/50','total_mass':str(M_cluster),'bad_order_area_lower':str(mass_sort_lower),'bad_order_ratio_lower':str(mass_sort_lower/M_cluster),'good_free_order_area_upper_strict':str(free_cluster_upper)},
        'six_cycle_deletions':deletions,'six_cycle_proper_rank_two_checks':proper_rank_checks,'six_cycle_rank_root_proxy_upper':'8/3',
        'corner_interleaving_proxy_outside_exact':'13','corner_own_mass':'1',
        'elapsed_seconds':time.perf_counter()-start,
        'source_sha256':{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in (proof,Path(__file__).resolve())},
        'scope':'Structural sorting-class obstruction and a bad-interleaving example; not a free-order minproxy>3 or C2 lower bound.'}
result=ROOT/'computation/results/post-r40-sorting/exact-obstructions.json'
result.parent.mkdir(parents=True,exist_ok=True);result.write_text(json.dumps(output,indent=2)+'\n')
print(json.dumps(output,indent=2))
