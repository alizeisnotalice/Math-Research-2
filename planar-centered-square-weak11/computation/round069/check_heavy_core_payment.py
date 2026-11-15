#!/usr/bin/env python3
"""Exact fixed Γ/beta/domain checks for R69; no full halo/W/Wprime areas."""
from fractions import Fraction as F
from pathlib import Path
from itertools import combinations
import argparse,hashlib,json,os,sys,time
import sympy as sp


def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def bbox(p):
    return (min(x for x,y in p),max(x for x,y in p),min(y for x,y in p),max(y for x,y in p))


def span2(p):
    a,b,c,d=bbox(p)
    return max(b-a,d-c)**2


def area(r):
    return F(0) if r is None else max(F(0),r[1]-r[0])*max(F(0),r[3]-r[2])


def clip(a,b):
    if a is None or b is None:
        return None
    q=(max(a[0],b[0]),min(a[1],b[1]),max(a[2],b[2]),min(a[3],b[3]))
    return q if q[0]<=q[1] and q[2]<=q[3] else None


def contains(r,x,y):
    return r is not None and r[0]<x<r[1] and r[2]<y<r[3]


def union_area(rects):
    rects=[r for r in rects if area(r)>0]
    if not rects:
        return F(0)
    xx=sorted({r[i] for r in rects for i in (0,1)})
    yy=sorted({r[i] for r in rects for i in (2,3)})
    return sum(((xr-xl)*(yr-yl) for xl,xr in zip(xx,xx[1:]) for yl,yr in zip(yy,yy[1:]) if any(contains(r,(xl+xr)/2,(yl+yr)/2) for r in rects)),F(0))


def direct_gamma(C,K,core_boxes):
    # Separate predicate-cell integration: C AND (not K OR core support).
    rects=[C]+([K] if K is not None else [])+core_boxes
    xx=sorted({r[i] for r in rects for i in (0,1)})
    yy=sorted({r[i] for r in rects for i in (2,3)})
    result=F(0)
    for xl,xr in zip(xx,xx[1:]):
        for yl,yr in zip(yy,yy[1:]):
            x=(xl+xr)/2;y=(yl+yr)/2
            if contains(C,x,y) and (not contains(K,x,y) or any(contains(r,x,y) for r in core_boxes)):
                result+=(xr-xl)*(yr-yl)
    return result


def root_gt(m,q):
    # Exact comparison sqrt(m)>q, with m>0.
    return q<=0 or m>q*q


def receiver_meets_closed_K(points,m,K):
    if K is None:
        return False
    l,r,b,t=bbox(points)
    return all(root_gt(m,z) for z in [2*(r-K[1]),2*(K[0]-l),2*(t-K[3]),2*(K[2]-b)])


def mst(points):
    edges=sorted((abs(points[i][0]-points[j][0])+abs(points[i][1]-points[j][1])-abs((points[i][0]-points[j][0])*(points[i][1]-points[j][1])),i,j) for i,j in combinations(range(len(points)),2))
    parent=list(range(len(points)))
    def root(i):
        while parent[i]!=i:
            parent[i]=parent[parent[i]];i=parent[i]
        return i
    selected=[]
    for cost,i,j in edges:
        a,b=root(i),root(j)
        if a!=b:
            parent[a]=b;selected.append((i,j,cost))
    assert len(selected)==len(points)-1
    return sum((q[2] for q in selected),F(0))


def analyze(case):
    p=[tuple(map(F,z)) for z in case['points']]
    mass=list(map(F,case['masses']));n=len(p);S=case['S_mask'];H=case['H_mask'];V=S&~H
    assert H and H&S==H and len(set(p))==n and min(mass)>0
    sids=[i for i in range(n) if S>>i&1];hids=[i for i in range(n) if H>>i&1];vids=[i for i in range(n) if V>>i&1]
    L=sum(mass[i] for i in sids);M=sum(mass[i] for i in hids);E=L-M;r=F(case['radius'])
    assert L==1 and E>=0 and 4*r*r<=3*M
    eligible=[]
    for mask in range(1,1<<n):
        ids=[i for i in range(n) if mask>>i&1];m=sum(mass[i] for i in ids)
        if span2([p[i] for i in ids])<m:
            eligible.append((mask,m))
    maximum=max(m for mask,m in eligible);max_masks=[mask for mask,m in eligible if m==maximum]
    assert maximum==L and max_masks==[S]
    seen=1
    while True:
        old=seen
        for mask,m in eligible:
            if seen&mask:
                seen|=mask
        if seen==old:
            break
    assert seen==(1<<n)-1
    nonterminal=span2(p)-F(35,33)*sum(mass);assert nonterminal>0
    left,right,bottom,top=bbox([p[i] for i in sids]);w=right-left;h=top-bottom
    assert w>F(3,4) and h>F(3,4) and max(w,h)<1
    C=(left-F(1,2),right+F(1,2),bottom-F(1,2),top+F(1,2))
    kl=max(p[i][0]-r for i in hids);kr=min(p[i][0]+r for i in hids)
    kb=max(p[i][1]-r for i in hids);kt=min(p[i][1]+r for i in hids)
    K=(kl,kr,kb,kt) if kl<=kr and kb<=kt else None
    core_boxes=[(p[i][0]-F(1,2),p[i][0]+F(1,2),p[i][1]-F(1,2),p[i][1]+F(1,2)) for i in hids]
    assert all(clip(z,C)==z for z in core_boxes)
    cap=clip(C,K);zcap=union_area([clip(z,K) for z in core_boxes]);Gamma=area(C)-area(cap)+zcap
    assert Gamma==direct_gamma(C,K,core_boxes)==F(case['expected_Gamma'])
    near=[];profiles=[]
    for q in range(n):
        if S>>q&1:
            continue
        d2=min((max(abs(p[q][0]-p[i][0]),abs(p[q][1]-p[i][1]))**2 for i in vids),default=None)
        member=d2 is not None and d2<E
        if member:
            near.append(q)
        profiles.append(dict(source=q,distance_to_V_squared=str(d2) if d2 is not None else None,strict_near=member))
    B=sum(1<<i for i in near);beta=sum((mass[i] for i in near),F(0));nu=V|B
    assert sum(mass[i] for i in range(n) if nu>>i&1)==E+beta
    branch={'touch_core':0,'low_miss_core':0,'high_miss_core':0};rows=[]
    for mask,m in eligible:
        if not mask&S:
            continue
        ids=[i for i in range(n) if mask>>i&1]
        if mask&H:
            category='touch_core'
        elif m<=E:
            category='low_miss_core';assert mask&~nu==0
            for q in ids:
                if not S>>q&1:
                    assert any(max(abs(p[q][0]-p[v][0]),abs(p[q][1]-p[v][1]))**2<E for v in vids if mask>>v&1)
        else:
            category='high_miss_core';assert M+m>L
            assert span2([p[i] for i in hids])<L<M+m
            assert span2([p[i] for i in ids])<m<M+m
            assert span2([p[i] for i in set(hids+ids)])>=M+m
            assert not receiver_meets_closed_K([p[i] for i in ids],m,K)
        branch[category]+=1;rows.append(dict(mask=mask,mass=str(m),branch=category))
    total=Gamma+9*(E+beta);assert total==F(case['expected_total9']) and total<F(7,2)*L
    low_rank=len(vids)*E;refined=Gamma+min(9*(E+beta),low_rank)
    tree=mst([p[i] for i in sids]);assert tree==F(case['expected_MST']) and tree>F(5,2)
    corner=False
    if len(hids)==1:
        x,y=p[hids[0]];corner=x in (left,right) and y in (bottom,top)
    uniform=corner and E+beta<=L/32
    if case['id']=='r068-eight-point-heavy-corner':
        assert uniform
    return dict(id=case['id'],source_count=n,subset_count=(1<<n)-1,eligible_subsets=len(eligible),unique_maximum_masks=max_masks,S_mask=S,H_mask=H,L=str(L),M=str(M),E=str(E),r=str(r),radius_margin_3M_minus_4r2=str(3*M-4*r*r),near_outside_indices=near,beta=str(beta),near_profiles=profiles,connected=True,nonterminal_margin=str(nonterminal),S_width=str(w),S_height=str(h),C=list(map(str,C)),K=list(map(str,K)) if K is not None else None,K_contained_in_C=K is not None and clip(K,C)==K,C_area=str(area(C)),C_intersect_K_area=str(area(cap)),ZH_intersect_K_area=str(zcap),Gamma_exact=str(Gamma),low_mass_weak9=str(9*(E+beta)),total9_upper=str(total),gap_to_7over2=str(F(7,2)*L-total),low_rank_upper=str(low_rank),min_low_refined_total=str(refined),uniform_dominant_corner_corollary_applies=uniform,uniform_corner_constant='2767/800' if uniform else None,MST_cost=str(tree),MST_5over2_test_passes=False,branch_counts=branch,all_touching_eligible_labels=rows,W_or_Wprime_area_evaluated=False,scope='Exact certificate for the proved heavy-core upper bound only, not the actual Wprime value.')


def main():
    assert os.path.ismount('/Volumes/data');start=time.perf_counter()
    ap=argparse.ArgumentParser();ap.add_argument('--input',type=Path,required=True);ap.add_argument('--out',type=Path,required=True);args=ap.parse_args()
    data=json.loads(args.input.read_text());x=sp.symbols('x')
    identities=[sp.Rational(4)-sp.Rational(329,400)+sp.Rational(9,32)-sp.Rational(2767,800),sp.Rational(93,32)-sp.Rational(289,100)-sp.Rational(13,800),sp.expand((3*x)**2-9*x*x)]
    assert all(v==0 for v in identities)
    out=dict(round=69,status='PASS_FIXED_EXACT_HEAVY_CORE_GEOMETRY',python=sys.version,executable=sys.executable,sympy=sp.__version__,process_pid=os.getpid(),script_sha256=sha(Path(__file__)),input_sha256=sha(args.input),symbolic_zero_identities=list(map(str,identities)),cases=[analyze(c) for c in data['cases']],elapsed_seconds=time.perf_counter()-start,scope='Only three existing fixed fixtures: all-subset domain/branch predicates, strict beta, rational clipped Gamma geometry and conditional upper-bound arithmetic. No fullhalo/W/Wprime areas, no LP, MATLAB, new input family or search.')
    assert str(args.out.resolve()).startswith('/Volumes/data/');args.out.parent.mkdir(parents=True,exist_ok=True);args.out.write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(dict(status=out['status'],python=out['python'],sympy=out['sympy'],process_pid=out['process_pid'],cases=[{k:c[k] for k in ['id','M','E','beta','Gamma_exact','total9_upper','K_contained_in_C','branch_counts','uniform_dominant_corner_corollary_applies']} for c in out['cases']],seconds=out['elapsed_seconds']),indent=2))

if __name__=='__main__':
    main()
