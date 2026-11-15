from fractions import Fraction as F
from itertools import combinations, product
from pathlib import Path
import hashlib, json, platform, sys

def rect(lx,ux,ly,uy):
    return (lx,ux,ly,uy) if lx<ux and ly<uy else None

def area(R):return F(0) if R is None else (R[1]-R[0])*(R[3]-R[2])

def inside(p,R):return R is not None and R[0]<p[0]<R[1] and R[2]<p[1]<R[3]

def union_area(rects):
    rects=[R for R in rects if R is not None]
    xs=sorted({x for R in rects for x in R[:2]});out=F(0)
    for a,b in zip(xs,xs[1:]):
        xm=(a+b)/2;ints=sorted((R[2],R[3]) for R in rects if R[0]<xm<R[1])
        if not ints:continue
        lo,hi=ints[0];height=F(0)
        for l,u in ints[1:]:
            if l>hi:height+=hi-lo;lo,hi=l,u
            else:hi=max(hi,u)
        height+=hi-lo;out+=(b-a)*height
    return out

def subsets(n):return [J for h in range(1,n+1) for J in combinations(range(n),h)]

def bounds(boxes,J):return max(boxes[i][0] for i in J),min(boxes[i][1] for i in J),max(boxes[i][2] for i in J),min(boxes[i][3] for i in J)

def saturate(boxes,J):
    lx,ux,ly,uy=bounds(boxes,J)
    return tuple(i for i,B in enumerate(boxes) if B[0]<=lx and B[1]>=ux and B[2]<=ly and B[3]>=uy)

boxes=[(F(0),F(1,10),F(0),F(1,10)),(F(1,4),F(7,20),F(0),F(1,10))]
masses=[F(9,25),F(16,25)]
centers=[((B[0]+B[1])/2,(B[2]+B[3])/2) for B in boxes]
deltas=[max(abs(v-c[k]) for k in (0,1) for v in B[2*k:2*k+2]) for B,c in zip(boxes,centers)]
radii={(0,):F(3,10),(1,):F(2,5),(0,1):F(1,2)}
allJ=subsets(2);family=sorted({saturate(boxes,J) for J in allJ})
Ds=[];Es=[];Fs=[];sum_error=F(0)
for J in allJ:
    K=saturate(boxes,J);assert set(J)<=set(K) and bounds(boxes,J)==bounds(boxes,K) and saturate(boxes,K)==K
    mass=sum((masses[i] for i in J),F(0));r=radii[J];assert 4*r*r==mass
    lx,ux,ly,uy=bounds(boxes,J);D=rect(lx-r,ux+r,ly-r,uy+r);Ds.append(D)
    cx=[centers[i][0] for i in J];cy=[centers[i][1] for i in J];delta=max(deltas[i] for i in J)
    E=rect(max(cx)-r,min(cx)+r,max(cy)-r,min(cy)+r);Es.append(E)
    G=rect(max(cx)-r-delta,min(cx)+r+delta,max(cy)-r-delta,min(cy)+r+delta);Fs.append(G)
    if D is not None:
        assert G[0]<=D[0] and D[1]<=G[1] and G[2]<=D[2] and D[3]<=G[3]
    err=area(G)-area(E);assert 0<=err<=8*delta*r+4*delta*delta;sum_error+=err
    if E is not None:assert D[0]<=E[0] and E[1]<=D[1] and D[2]<=E[2] and E[3]<=D[3]
A_D=union_area(Ds);A_C=union_area(Es);gap=A_D-A_C;assert 0<=gap<=sum_error
# Eight atoms, each block's four corners, independent exact strict-capture tests.
atoms=[((x,y),masses[i]/4) for i,B in enumerate(boxes) for x in B[:2] for y in B[2:]]
checks=0;strict_witnesses=0
for x,y,k in product(range(-10,21),range(-10,21),range(1,13)):
    p=(F(x,20),F(y,20));r=F(k,20);caught=[(a,m) for a,m in atoms if max(abs(a[0]-p[0]),abs(a[1]-p[1]))<=r]
    mass=sum((m for a,m in caught),F(0));checks+=1
    if 4*r*r<mass:
        strict_witnesses+=1;assert any(inside(p,D) for D in Ds)
assert strict_witnesses>0
# All saturation subsets for overlapping boxes; mass ownership remains labels.
overlap=[(F(0),F(2),F(0),F(2)),(F(1),F(3),F(1),F(3)),(F(1,2),F(5,2),F(1,2),F(5,2))]
saturation_checks=0
for J in subsets(3):
    K=saturate(overlap,J);assert set(J)<=set(K) and bounds(overlap,J)==bounds(overlap,K) and saturate(overlap,K)==K;saturation_checks+=1
assert saturate(overlap,(0,1))==(0,1,2)
# Threshold branch: original contracted receiver empty in x, expanded one positive.
r=F(1,2);delta=F(1,20);E=rect(F(1)-r,r,-r,r);G=rect(F(1)-r-delta,r+delta,-r-delta,r+delta)
assert E is None and area(G)==F(11,100) and area(G)<=8*delta*r+4*delta*delta
root=Path(__file__).resolve().parents[1]
proof=root/'proofs/r053-quantitative-bbox-contraction.md'
out=dict(status='PASS',python=platform.python_version(),executable=sys.executable,arithmetic='Fraction exact; no random seed, tolerance, float, or MATLAB',blocks=2,atoms=8,total_mass='1',representative_radius=str(max(deltas)),block_subsets=len(allJ),saturated_states=len(family),saturation_overlap_subsets=saturation_checks,receiver_point_radius_checks=checks,strict_witnesses=strict_witnesses,outer_envelope_area=str(A_D),contracted_halo_area=str(A_C),exact_union_excess=str(gap),sum_individual_geometric_errors=str(sum_error),threshold_empty_then_positive_area=str(area(G)),script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),proof_sha256=hashlib.sha256(proof.read_bytes()).hexdigest())
print(json.dumps(out,indent=2))
