from fractions import Fraction as F
from itertools import combinations,product
from pathlib import Path
import json,hashlib,platform,sys
B=F(1994727599733432151936139,604462909807314587353088)
K=F(7,2)
assert B<F(33,10) and B*F(35,33)<K

def area(R):return (R[1]-R[0])*(R[3]-R[2])
def intersect(R,S):
    T=(max(R[0],S[0]),min(R[1],S[1]),max(R[2],S[2]),min(R[3],S[3]))
    return T if T[0]<T[1] and T[2]<T[3] else None

def union_area(rects):
    rects=[R for R in rects if R is not None];xs=sorted({x for R in rects for x in R[:2]});out=F(0)
    for a,b in zip(xs,xs[1:]):
        x=(a+b)/2;ints=sorted((R[2],R[3]) for R in rects if R[0]<x<R[1])
        if not ints:continue
        l,u=ints[0];height=F(0)
        for ll,uu in ints[1:]:
            if ll>u:height+=u-l;l,u=ll,uu
            else:u=max(u,uu)
        out+=(b-a)*(height+u-l)
    return out

def bbox(atoms):return (min(p[0] for p,m in atoms),max(p[0] for p,m in atoms),min(p[1] for p,m in atoms),max(p[1] for p,m in atoms))
def diam2(atoms):
    R=bbox(atoms);return max(R[1]-R[0],R[3]-R[2])**2

def mass(atoms):return sum((m for p,m in atoms),F(0))
q=3;eta=F(1,100);sides=[]
for c in (-F(11,10),F(11,10)):
    sides.append([((c-eta/2+eta*F(u,q-1),-eta/2+eta*F(v,q-1)),F(1,q*q)) for u,v in product(range(q),repeat=2)])
pivot=((F(0),F(0)),F(1));atoms=[pivot]+sides[0]+sides[1]
assert len(atoms)==19 and mass(atoms)==3 and diam2(atoms)==F(48841,10000)>3
assert all(mass(C)==1 and diam2(C)==eta*eta<1 for C in sides)
assert all(diam2([pivot]+C)<2 for C in sides)
# Any crossing source set has span >=219/100, squared>total root mass3.
crossmin=min(p[0] for p,m in sides[1])-max(p[0] for p,m in sides[0])
assert crossmin==F(219,100) and crossmin*crossmin>mass(atoms)
# Enumerate every complete source bbox, not receiver/radius grid.
xs=sorted({p[0] for p,m in atoms});ys=sorted({p[1] for p,m in atoms})
cap=F(1);boxes=0;pivot_boxes=0;eligible_pivot_boxes=0
for ia,ib,ic,id in product(range(len(xs)),range(len(xs)),range(len(ys)),range(len(ys))):
    if ia>ib or ic>id:continue
    boxes+=1
    if not(xs[ia]<=0<=xs[ib] and ys[ic]<=0<=ys[id]):continue
    pivot_boxes+=1
    S=[(p,m) for p,m in atoms if xs[ia]<=p[0]<=xs[ib] and ys[ic]<=p[1]<=ys[id]]
    if diam2(S)<mass(S):eligible_pivot_boxes+=1;cap=max(cap,mass(S))
assert cap==2
# Safe rational outer pivot rectangle; not a false claim sqrt(2) is rational.
assert F(3,2)**2>cap
P=(-F(3,4),F(3,4),-F(3,4),F(3,4));cores=[]
for C in sides:
    R=bbox(C);G=(R[1]-F(1,2),R[0]+F(1,2),R[3]-F(1,2),R[2]+F(1,2));cores.append(G)
    assert area(G)==F(9801,10000) and area(intersect(P,G))==F(2871,20000)
    # q3 singleton halos have halfside1/6 and are contained in full-leaf core.
    assert F(1,3)**2==C[0][1]
    for p,m in C:
        assert G[0]<=p[0]-F(1,6) and p[0]+F(1,6)<=G[1]
        assert G[2]<=p[1]-F(1,6) and p[1]+F(1,6)<=G[3]
assert intersect(*cores) is None
Parea=union_area([P]);union=union_area([P]+cores);sumcores=sum(map(area,cores),F(0))
omega=Parea+sumcores-union
assert Parea==F(9,4) and omega==F(2871,10000)
score=union+2*B-sumcores
conservative=F(85629,10000)
assert score==Parea+2*B-omega<conservative<K*mass(atoms)
assert B*diam2(atoms)>K*mass(atoms)
# Explicit global tree budget (8), with all leaves eligible and total leaf mass2.
lhs=Parea-K;rhs=(K-B)*2+omega
assert score-K*3==lhs-rhs and lhs<rhs
root=Path(__file__).resolve().parents[1];proof=root/'proofs/r058-bbox-terminal-masked-tree.md'
out=dict(status='PASS',python=platform.python_version(),executable=sys.executable,arithmetic='exact Fraction; no float/randomness/tolerance/MATLAB',atoms=len(atoms),side_grid_q=q,total_mass='3',all_complete_source_boxes=boxes,pivot_containing_source_boxes=pivot_boxes,eligible_pivot_boxes=eligible_pivot_boxes,root_bbox_diameter_squared=str(diam2(atoms)),pivot_cap=str(cap),safe_pivot_halfside='3/4',pivot_union_area=str(Parea),leaf_core_areas=[str(area(G)) for G in cores],all_cover_core_union_area=str(union),overlap_credit=str(omega),certified_score=str(score),strict_rational_score_bound=str(conservative),target_total_area='21/2',budget_lhs=str(lhs),budget_rhs=str(rhs),script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),proof_sha256=hashlib.sha256(proof.read_bytes()).hexdigest())
print(json.dumps(out,indent=2))
