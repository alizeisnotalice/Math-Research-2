#!/usr/bin/env python3
"""R69 seven-source net-mixed obstruction: exact subsets and rational patch.

No full-halo numerical evaluation. Rectangle bounds are analytic outer covers
of the internal and residual halos, and a fixed strict witness covers the
two inner rectangles. Their rational difference is a certified pure patch.
"""
import argparse
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path
import platform
import time


def intersect(a, b):
    q = (max(a[0],b[0]),min(a[1],b[1]),max(a[2],b[2]),min(a[3],b[3]))
    return q if q[0] < q[1] and q[2] < q[3] else None


def area(q):
    return (q[1]-q[0])*(q[3]-q[2]) if q else F(0)


def square(p, r):
    return (p[0]-r,p[0]+r,p[1]-r,p[1]+r)


def difference_area(inners, outers):
    xs = sorted({q[k] for q in inners+outers for k in (0,1)})
    ys = sorted({q[k] for q in inners+outers for k in (2,3)})
    total = F(0)
    def inside(x,y,q):
        return q[0] < x < q[1] and q[2] < y < q[3]
    for xa,xb in zip(xs,xs[1:]):
        x=(xa+xb)/2
        for ya,yb in zip(ys,ys[1:]):
            y=(ya+yb)/2
            if any(inside(x,y,q) for q in inners) and not any(inside(x,y,q) for q in outers):
                total+=(xb-xa)*(yb-ya)
    return total


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--input',type=Path,required=True);ap.add_argument('--output',type=Path,required=True)
    args=ap.parse_args();start=time.perf_counter();inp=json.loads(args.input.read_text())
    points=[tuple(map(F,row['point'])) for row in inp['sources']]
    masses=[F(row['mass']) for row in inp['sources']]
    n=len(points);S=inp['S_mask'];assert n==7 and S==31 and all(m>0 for m in masses)
    allrows=[];eligible=[];closed=[]
    for mask in range(1,1<<n):
        ids=[i for i in range(n) if mask>>i&1]
        mass=sum((masses[i] for i in ids),F(0))
        lo=[min(points[i][t] for i in ids) for t in (0,1)]
        hi=[max(points[i][t] for i in ids) for t in (0,1)]
        widths=[hi[t]-lo[t] for t in (0,1)]
        closure=sum(1<<i for i,p in enumerate(points) if all(lo[t]<=p[t]<=hi[t] for t in (0,1)))
        row={'mask':mask,'mass':str(mass),'widths':list(map(str,widths)),
             'eligible':max(widths)**2<mass,'closed':closure==mask}
        allrows.append(row)
        if row['eligible']:eligible.append(row)
        if row['closed']:closed.append(row)
    L=max(F(r['mass']) for r in eligible)
    maxima=[r['mask'] for r in eligible if F(r['mass'])==L]
    assert L==1 and maxima==[S]
    fullmass=sum(masses,F(0));span=max(F(v) for v in allrows[-1]['widths'])
    wide_margin=span*span-F(35,33)*fullmass;assert wide_margin>0
    parent=list(range(n))
    def find(i):
        while parent[i]!=i:
            parent[i]=parent[parent[i]];i=parent[i]
        return i
    for row in eligible:
        ids=[i for i in range(n) if row['mask']>>i&1]
        for i in ids[1:]:parent[find(i)]=find(ids[0])
    assert len({find(i) for i in range(n)})==1
    # A separate exact Prim calculation at L=1.
    vertices=[i for i in range(n) if S>>i&1];taken={vertices[0]};tree=[];treecost=F(0)
    while len(taken)<len(vertices):
        candidates=[]
        for i in taken:
            for j in vertices:
                if j not in taken:
                    dx,dy=(abs(points[i][t]-points[j][t]) for t in (0,1))
                    candidates.append((dx+dy-dx*dy,i,j))
        cost,i,j=min(candidates);taken.add(j);tree.append([i,j]);treecost+=cost
    assert treecost==F(27,10)>F(5,2)
    tiny=sum(masses[2:5],F(0));oneheavy=masses[0]+tiny
    rI,rV,rR=F(177,500),F(23,1000),F(353,1000)
    assert tiny==F(1,500) and oneheavy==F(501,1000)
    assert 4*rI*rI>oneheavy and 4*rV*rV>tiny and 4*rR*rR>masses[5]
    # All internal labels contain two heavy sources, one, or none.
    Icover=[square(points[i],rI) for i in (0,1)]
    Icover.append((F(0),F(1,2),F(-1,2),F(1,2)))
    Icover.extend(square(points[i],rV) for i in (2,3,4))
    Rcover=[square(points[i],rR) for i in (5,6)]
    assert all(intersect(a,b) is None for a in Icover for b in Rcover)
    assert (points[6][0]-points[5][0])**2>masses[5]+masses[6]
    r=F(499,1000)
    J=[(F(-499,1000),F(-301,1000),-r,r),
       (F(801,1000),F(999,1000),-r,r)]
    for rect,ids in zip(J,[(0,5),(1,6)]):
        assert sum((masses[i] for i in ids),F(0))>4*r*r
        for i in ids:
            for t in (0,1):
                assert max(abs(rect[2*t]-points[i][t]),abs(rect[2*t+1]-points[i][t]))<=r
    paid=Icover+Rcover
    patch=difference_area(J,paid)
    assert patch==F(12231,50000) and patch>F(1,5)
    # Exact evaluation of the new heavy-core sufficient criterion, not a
    # separate proof of that general lemma (see its proof attachment).
    radius=F(43,50);M=masses[0]+masses[1]
    assert 4*radius*radius<=3*M
    assert all(min(max(abs(points[q][t]-points[v][t]) for t in (0,1)) for v in (2,3,4))**2>=tiny for q in (5,6))
    C=(F(-1,2),F(7,5),F(-1,2),F(7,5))
    K=(F(1,2)-radius,radius,-radius,radius)
    Z=[square(points[i],F(1,2)) for i in (0,1)]
    Kinside=intersect(C,K)
    ZK=area(Kinside)-difference_area([Kinside],Z)
    gamma=area(C)-area(Kinside)+ZK
    bound=gamma+9*tiny
    assert bound==F(1993,625)<F(7,2)
    result={'schema':'r069-net-mixed-obstruction-v1','status':'PASS',
            'scope':'127 exact subsets and a rational positive-area patch; no complete halo area evaluation. Refutes net-mixed 1/5 increment, not the 3.5 deletion/Wprime selection.',
            'python':platform.python_version(),'platform':platform.platform(),
            'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
            'input_sha256':hashlib.sha256(args.input.read_bytes()).hexdigest(),
            'source_count':n,'subset_count':len(allrows),'closed_count':len(closed),
            'eligible_count':len(eligible),'eligible_closed_count':sum(r['eligible'] for r in closed),
            'maximum_mass':str(L),'maximum_masks':maxima,'connected':True,
            'nonterminal_margin':str(wide_margin),'mst_cost':str(treecost),'mst_edges':tree,
            'internal_residual_disjoint':True,'pure_mixed_patch_area':str(patch),
            'net_mixed_lower':str(patch),'excess_over_one_fifth':str(patch-F(1,5)),
            'fixed_witness_radius':str(r),'witness_mass':str(masses[0]+masses[5]),
            'witness_cost':str(4*r*r),
            'I_outer_rectangles':[[str(v) for v in q] for q in Icover],
            'R_outer_rectangles':[[str(v) for v in q] for q in Rcover],
            'full_halo_inner_rectangles':[[str(v) for v in q] for q in J],
            'heavy_core_criterion':{'H_mask':3,'M':str(M),'E':str(tiny),'beta':'0','r':str(radius),'Gamma':str(gamma),'bound':str(bound),'paid':True},
            'all_subsets':allrows,'wall_seconds':time.perf_counter()-start}
    args.output.parent.mkdir(parents=True,exist_ok=True);args.output.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:result[k] for k in ['status','subset_count','closed_count','eligible_closed_count','pure_mixed_patch_area','excess_over_one_fifth','heavy_core_criterion','wall_seconds']}))


if __name__=='__main__':main()
