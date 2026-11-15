from fractions import Fraction as F
from itertools import combinations
from pathlib import Path
import json,hashlib,platform

def span(points):
    if not points:return F(0)
    return max(max(p[k] for p in points)-min(p[k] for p in points) for k in (0,1))
def subsets(n):
    for k in range(n+1):
        yield from combinations(range(n),k)
def components(points,masses):
    parent=list(range(len(points)));edges=[]
    def find(i):
        while parent[i]!=i:i=parent[i]
        return i
    for S in subsets(len(points)):
        if not S:continue
        if span([points[i] for i in S])**2<sum((masses[i] for i in S),F(0)):
            edges.append(S)
            for i in S[1:]:parent[find(i)]=find(S[0])
    groups={}
    for i in range(len(points)):groups.setdefault(find(i),[]).append(i)
    return list(groups.values()),edges
# Corner cluster and failure of arbitrary internal-source accounting.
p=(F(0),F(0));old=[(F(1),F(0)),(F(4),F(0))];masses=[F(1),F(14)];b=F(1)
eligible=[]
for S in subsets(2):
    if span([p]+[old[i] for i in S])**2<sum((masses[i] for i in S),b):eligible.append(S)
assert eligible==[(),(0,)]
assert span([p,old[1]])**2==b+sum(masses)
cc,edges=components(old,masses);assert cc==[[0,1]]
# Q at q has area15 and contains Q at p of area2: sqrt(15)>2+sqrt(2).
assert F(9,4)**2>2 # 6+4sqrt(2)<15.
corner=dict(points=[[str(z) for z in q] for q in old],masses=[str(z) for z in masses],eligible_with_p=[list(S) for S in eligible],old_components=cc,internal_proxy_outer_loss='13')
# General overlapping eligible labels are not union-closed.
p2=(F(9,10),F(0));old2=[(F(0),F(0)),(F(47,50),F(0)),(F(46,25),F(0))]
for S in [(0,1),(1,2)]:assert span([p2]+[old2[i] for i in S])**2<F(3,2)
assert span([p2]+old2)**2>2
# The 8-light + heavy separator example, all old nonempty subsets.
e=F(1,10000);d=F(99,100);light=[(x,y) for x in (-d,F(0),d) for y in (-d,F(0),d) if (x,y)!=(F(0),F(0))]
cc,edges=components(light,[e]*8);assert len(cc)==8 and len(edges)==8
maxold=F(0);eligible_count=0
for S in subsets(8):
    m=len(S)*e
    if span([p]+[light[i] for i in S])**2<1+m:maxold=max(maxold,m);eligible_count+=1
assert maxold==3*e
assert F(1001,1000)**2>1+3*e and d>(F(1001,1000)+F(1,100))/2
assert 1+11*e==F(10011,10000)
assert (F(197,50)/(1+8*e))==F(4925,1251)
for u,v in [(1,0),(0,1),(1,1),(1,-1),(2,1)]:assert sum(u*x+v*y>0 for x,y in light)>=3
# Mass-dominant leaf: M16, p0 of mass15, eight old unit-perimeter sources totaling1.
oldmass=F(1);totalmass=F(16);D=F(1)
assert D+F(1,2)<=F(2) # D+sqrt(R)/2 <= sqrt(M)/2.
assert (2*D)**2<totalmass
# Uniform nesting leaf: 8 unit masses over bbox diameter1/1000.
nesting_D=F(1,1000);nesting_M=F(8);minmass=F(1)
assert 16*nesting_D**2*nesting_M<minmass**2
leaf_checks=dict(mass_dominant=dict(total='16',last_mass='15',old_mass='1',D='1',old_proxy_outer_radius_bound='3/2',final_radius='2'),uniform_nesting=dict(total='8',min_mass='1',bbox_diameter='1/1000',squared_test=str(16*nesting_D**2*nesting_M)))
# Six-point eligible pair cycle; every deletion has one noneligible 5-point component.
a=F(49,100);cycle=[(a,F(0)),(a/2,a),(-a/2,a),(-a,F(0)),(-a/2,-a),(a/2,-a)]
pair_edges=[]
for i,j in combinations(range(6),2):
    valid=span([cycle[i],cycle[j]])**2<F(1,3)
    expected=(j-i in (1,5))
    assert valid==expected
    if valid:pair_edges.append([i,j])
deletions=[]
for gone in range(6):
    pts=[cycle[i] for i in range(6) if i!=gone]
    comps,_=components(pts,[F(1,6)]*5)
    assert len(comps)==1 and len(comps[0])==5
    assert span(pts)**2==F(2401,2500)>F(5,6)
    deletions.append(dict(deleted=gone,component_sizes=[5],whole_component_eligible=False))
assert span(cycle)**2<F(1)
proper_eligible=[]
for S in subsets(6):
    if not S or len(S)==6:continue
    if span([cycle[i] for i in S])**2<F(len(S),6):proper_eligible.append(S)
assert max(map(len,proper_eligible))==2 and len(proper_eligible)==12
assert F(1)+2*F(5,6)==F(8,3)
cycle_check=dict(a=str(a),pair_edges=pair_edges,deletions=deletions,whole_six_eligible=True,proper_eligible_rank=2,proper_eligible_count=len(proper_eligible),rank_root_proxy_bound='8/3')
# Arbitrarily long connected rank2 class: exact finite exemplar N8.
chain=[(F(11*i,10),F(0)) for i in range(8)]
chain_comps,chain_edges=components(chain,[F(1)]*8)
assert len(chain_comps)==1 and len(chain_edges)==15 and max(map(len,chain_edges))==2
rank_check=dict(chain_N=8,chain_components=chain_comps,chain_eligible_rank=2,chain_eligible_count=15,general_rank_k_new_mass_payment='T_j<=k*m_j')
# One exact quarter-strip container example, no large geometry enumeration.
strip_points=[(F(0),F(0)),(F(49,50),F(0)),(F(49,50),F(49,50))]
strip_masses=[F(1,4),F(1,4),F(1,2)];cut=F(1,4)
R=sum((m for q,m in zip(strip_points,strip_masses) if q[0]<cut),F(0))
assert R==F(1,4) and cut+F(1,4)==F(1,2)
container=(F(-1,4),F(37,25),F(-1,2),F(37,25))
old=[];proxy_caps=[];proxy_rects=[]
for i,q in enumerate(strip_points):
    cap=strip_masses[i]
    for S in subsets(len(old)):
        idx=[old[z] for z in S]+[i]
        m=sum((strip_masses[z] for z in idx),F(0))
        if span([strip_points[z] for z in idx])**2<m:cap=max(cap,m)
    proxy_caps.append(cap);old.append(i)
    assert cap in (F(1,4),F(1))
    radius=F(1,4) if cap==F(1,4) else F(1,2)
    rect=(q[0]-radius,q[0]+radius,q[1]-radius,q[1]+radius)
    assert container[0]<=rect[0]<=rect[1]<=container[1] and container[2]<=rect[2]<=rect[3]<=container[3]
    proxy_rects.append(rect)
assert proxy_caps==[F(1,4),F(1,4),F(1)]
for a,b in combinations(proxy_rects,2):assert min(a[1],b[1])<=max(a[0],b[0]) or min(a[3],b[3])<=max(a[2],b[2])
strip_bound=(container[1]-container[0])*(container[3]-container[2])
assert strip_bound==F(17127,5000)<F(7,2)
strip_check=dict(M='1',a='1/4',R=str(R),W='49/50',H='49/50',proxy_caps=[str(z) for z in proxy_caps],container=[str(z) for z in container],container_area=str(strip_bound),actual_disjoint_proxy_area='3/2')
root=Path(__file__).resolve().parents[1]
source_names=['proofs/sol-effective-envelope.md','proofs/sol-effective-ordering-draft.md','proofs/general-dual-interface.md','proofs/mass-capacity-equality-attempt.md']
print(json.dumps(dict(python=platform.python_version(),arithmetic='exact Fraction, deterministic single process; no seed or tolerance',corner_example=corner,noncorner_union_failure=True,light_components=cc,light_eligible_edges=len(edges),heavy_last_max_old=str(maxold),heavy_last_proxy_area=str(1+11*e),fixed_projection_ratio_lower=str(F(4925,1251)),eligible_light_subsets_with_heavy=eligible_count,leaf_checks=leaf_checks,six_cycle=cycle_check,rank_checks=rank_check,quarter_strip=strip_check,source_sha256={s:hashlib.sha256((root/s).read_bytes()).hexdigest() for s in source_names},proof_sha256=hashlib.sha256((root/'proofs/sol-effective-separator-recursion.md').read_bytes()).hexdigest()),indent=2))
