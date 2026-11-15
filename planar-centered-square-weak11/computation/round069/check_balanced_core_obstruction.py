#!/usr/bin/env python3
"""255 exact subsets, 15 outward-radius core lower bounds, four whole patches."""
from fractions import Fraction as F
from pathlib import Path
import hashlib
import itertools
import json
import os
import sys
import time


def box(pts):
    return (min(x for x,y in pts),max(x for x,y in pts),min(y for x,y in pts),max(y for x,y in pts))


def inside(r,x,y):
    return r[0]<x<r[1] and r[2]<y<r[3]


def measure(rectangles,predicate):
    xs=sorted({r[k] for r in rectangles for k in (0,1)})
    ys=sorted({r[k] for r in rectangles for k in (2,3)})
    total=F(0)
    for a,b in zip(xs,xs[1:]):
        for c,d in zip(ys,ys[1:]):
            if predicate((a+b)/2,(c+d)/2):total+=(b-a)*(d-c)
    return total


def square(p,r):
    x,y=p
    return (x-r,x+r,y-r,y+r)


def rotate(p):
    x,y=p
    return (F(9,10)-y,x)


def rotate_rect(r):
    a,b,c,d=r
    return (F(9,10)-d,F(9,10)-c,a,b)


def main():
    if not __debug__:raise RuntimeError('Assertions must be enabled')
    assert os.path.ismount('/Volumes/data')
    started=time.perf_counter()
    pts=[(F(0),F(0)),(F(9,10),F(0)),(F(9,10),F(9,10)),(F(0),F(9,10))]
    q=(F(-3,4),F(-11,50))
    for j in range(4):pts.append(q);q=rotate(q)
    masses=[F(1,4)]*4+[F(749,1000)]*4
    S=15;total=sum(masses);rows=[]
    for mask in range(1,256):
        ids=[i for i in range(8) if mask>>i&1];a,b,c,d=box([pts[i] for i in ids])
        m=sum(masses[i] for i in ids);sp=max(b-a,d-c)**2
        rows.append({'mask':mask,'mass':str(m),'span_squared':str(sp),'eligible':m>sp})
    maximum=max(F(z['mass']) for z in rows if z['eligible'])
    maxima=[z['mask'] for z in rows if z['eligible'] and F(z['mass'])==maximum]
    assert maximum==1 and maxima==[S]
    assert rows[S-1]['eligible']
    for i in range(4):assert rows[((1<<i)|(1<<(i+4)))-1]['eligible']
    nonterminal=F(12,5)**2-F(35,33)*total
    assert nonterminal>0
    # Internal halo: only four singletons and the full S are eligible.
    internal_eligible=[z['mask'] for z in rows if z['mask']&~S==0 and z['eligible']]
    residual_eligible=[z['mask'] for z in rows if z['mask']&S==0 and z['eligible']]
    assert internal_eligible==[1,2,4,8,15] and residual_eligible==[16,32,64,128]
    I=[square(p,F(1,4)) for p in pts[:4]]+[(F(2,5),F(1,2),F(2,5),F(1,2))]
    rout=F(433,1000);assert 4*rout*rout>F(749,1000)
    Rplus=[square(p,rout) for p in pts[4:]]
    Icredit=measure(I+Rplus,lambda x,y:any(inside(t,x,y) for t in I) and any(inside(t,x,y) for t in Rplus))
    assert Icredit==0
    radius=F(999,2000);assert 4*radius**2<F(999,1000)
    J=[(-radius,radius-F(3,4),-radius,radius-F(11,50))]
    for i in range(3):J.append(rotate_rect(J[-1]))
    for i,r in enumerate(J):
        for j in (i,i+4):
            assert max(abs(pts[j][0]-r[0]),abs(pts[j][0]-r[1]),abs(pts[j][1]-r[2]),abs(pts[j][1]-r[3]))<=radius
    assert measure(J+I,lambda x,y:any(inside(t,x,y) for t in J) and any(inside(t,x,y) for t in I))==0
    patch=measure(J+Rplus,lambda x,y:any(inside(t,x,y) for t in J) and not any(inside(t,x,y) for t in Rplus))
    assert patch==F(255759,1000000) and patch>F(1,5)

    C=(F(-1,2),F(7,5),F(-1,2),F(7,5))
    # r_up >= sqrt(3*m(H))/2. Gamma is nonincreasing in r, so
    # computing Gamma(r_up) gives a LOWER bound at the optimal allowed r.
    radius_uppers={1:F(434,1000),2:F(613,1000),3:F(3,4),4:F(867,1000)}
    cores=[]
    for H in range(1,16):
        ids=[i for i in range(4) if H>>i&1];M=F(len(ids),4);r=radius_uppers[len(ids)]
        assert r*r>=3*M/4
        a,b,c,d=box([pts[i] for i in ids]);K=(b-r,a+r,d-r,c+r)
        Z=[square(pts[i],F(1,2)) for i in ids]
        gamma=measure([C,K]+Z,lambda x,y:inside(C,x,y) and (not inside(K,x,y) or any(inside(z,x,y) for z in Z)))
        assert gamma>=F(1419,400)>F(7,2)
        cores.append({'H_mask':H,'mass':str(M),'outward_radius':str(r),'radius_square_margin':str(r*r-3*M/4),'Gamma_at_outward_radius':str(gamma),'Phi_lower_bound':str(gamma)})
    assert min(F(v['Phi_lower_bound']) for v in cores)==F(1419,400)
    # The side path has minimum tree cost: all six edges cost at least9/10.
    costs=[abs(pts[i][0]-pts[j][0])+abs(pts[i][1]-pts[j][1])-abs((pts[i][0]-pts[j][0])*(pts[i][1]-pts[j][1])) for i,j in itertools.combinations(range(4),2)]
    assert min(costs)==F(9,10) and 3*min(costs)==F(27,10)>F(5,2)
    root=Path(__file__).resolve().parent
    spec={'id':'r069-balanced-four-corner-all-core-counterexample','points':[list(map(str,p)) for p in pts],'initial_masses':list(map(str,masses)),'S_mask':S,'threshold':'1','geometry':'axis-parallel centered common radius','patch_J_open':[[str(v) for v in r] for r in J],'fixed_half_side':str(radius),'residual_outer_half_side':str(rout),'core_upper_radii':{str(k):str(v) for k,v in radius_uppers.items()}}
    ip=root/'r069-balanced-core-obstruction-input.json';ip.write_text(json.dumps(spec,indent=2)+'\n')
    result={'schema':'r069-balanced-all-core-and-net-obstruction-v1','status':'PASS','python':sys.version,'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'input_sha256':hashlib.sha256(ip.read_bytes()).hexdigest(),'source_count':8,'complete_subset_count':255,'subsets':rows,'maximum_mass':'1','all_maximum_masks':maxima,'connected_witnesses':'S plus four adjacent source/residual pairs','full_mass':str(total),'nonterminal_margin':str(nonterminal),'MST_minimum':'27/10','internal_halo_labels':internal_eligible,'residual_halo_labels':residual_eligible,'internal_residual_outer_credit':'0','patch_area_net_lower':str(patch),'excess_above_one_fifth':str(patch-F(1,5)),'core_count':len(cores),'all_cores':cores,'minimum_Gamma_lower':'1419/400','strict_margin_to_7over2':'19/400','wall_seconds':time.perf_counter()-started,'scope':'All source subsets, unique maximum, all15 cores and rational whole-patch certificate refute the conditional finite-core/net-one-fifth alternative. No general3.5/C2 refutation or complete full halo area is asserted.'}
    (root/'r069-balanced-core-obstruction-arithmetic.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({'status':'PASS','net_lower':str(patch),'all_core_Gamma_lower':'1419/400','seconds':result['wall_seconds']}))


if __name__=='__main__':main()
