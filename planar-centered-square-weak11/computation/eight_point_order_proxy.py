#!/usr/bin/env python3
"""All orderings for the 8-point perimeter insertion proxy; exact certificates."""
from fractions import Fraction as F
from itertools import permutations, product
from pathlib import Path
import csv, hashlib, json, math, os, platform, sys, time

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'computation/results/eight-point-order-proxy'
BITS=32
A=F(49,100)
POINTS=tuple((x,y) for x,y in product((-A,F(0),A),repeat=2) if (x,y)!=(0,0))
N=len(POINTS)
FULL=(1<<N)-1

def bbox_span_squared(mask):
    pts=[POINTS[j] for j in range(N) if mask>>j&1]
    span=max(max(p[0] for p in pts)-min(p[0] for p in pts),
             max(p[1] for p in pts)-min(p[1] for p in pts))
    return span*span

def build_eff():
    spans={mask:bbox_span_squared(mask) for mask in range(1,1<<N)}
    eff=[[None]*(1<<N) for _ in range(N)]
    eligibility_checks=0
    for j in range(N):
        bit=1<<j
        for old in range(1<<N):
            if old&bit: continue
            k=old.bit_count()
            eligibility_checks+=1
            if spans[old|bit]<F(k+1,N):
                eff[j][old]=k
            else:
                eff[j][old]=max((eff[j][old^(1<<i)] for i in range(N) if old>>i&1),default=0)
    return eff,eligibility_checks

def union_area(rects):
    xs=sorted({v for r in rects for v in r[:2]})
    total=0
    for lo,hi in zip(xs,xs[1:]):
        # Closed endpoint membership is irrelevant to area. No midpoint rounding.
        intervals=sorted((y0,y1) for x0,x1,y0,y1 in rects if x0<=lo and hi<=x1)
        length=0; start=None; end=None
        for y0,y1 in intervals:
            if start is None: start,end=y0,y1
            elif y0>end: length+=end-start;start,end=y0,y1
            else: end=max(end,y1)
        if start is not None: length+=end-start
        total+=(hi-lo)*length
    return total

def build_rectangles():
    denominator=100*(1<<BITS)
    floats={};inners={};outers={};bounds={}
    for num in range(1,N+1):
        # r^2=(num/8)/4=num/32. isqrt gives a rigorous dyadic bracket.
        radicand=num << (2*BITS-5)
        low=math.isqrt(radicand)
        high=low if low*low==radicand else low+1
        assert low*low<=radicand<=high*high
        bounds[str(num)]={'radius_squared':str(F(num,32)),
                          'lower_radius_numerator':low,'upper_radius_numerator':high,
                          'radius_denominator':1<<BITS}
        radius=math.sqrt(num/N)/2
        for j,(x,y) in enumerate(POINTS):
            xf,yf=float(x),float(y)
            floats[j,num]=(xf-radius,xf+radius,yf-radius,yf+radius)
            xi=int(x*denominator); yi=int(y*denominator)
            assert F(xi,denominator)==x and F(yi,denominator)==y
            rlo=100*low;rhi=100*high
            inners[j,num]=(xi-rlo,xi+rlo,yi-rlo,yi+rlo)
            outers[j,num]=(xi-rhi,xi+rhi,yi-rhi,yi+rhi)
    return denominator,floats,inners,outers,bounds

def main():
    start=time.perf_counter()
    assert os.path.ismount('/Volumes/data')
    OUT.mkdir(parents=True,exist_ok=True)
    eff,eligibility_checks=build_eff()
    D,floats,inners,outers,bounds=build_rectangles()
    D2=D*D
    sigs={};records=[];orbit_stats={}
    last_options=(('corner',0),('edge_midpoint',1))
    for name,last in last_options:
        best_float=math.inf;best_sig=None;best_order=None
        count=0
        for before in permutations(j for j in range(N) if j!=last):
            order=before+(last,);mask=0;signature=[None]*N
            for j in order:
                signature[j]=eff[j][mask]+1
                mask|=1<<j
            assert mask==FULL
            signature=tuple(signature)
            if signature not in sigs:
                inner=union_area([inners[j,num] for j,num in enumerate(signature)])
                outer=union_area([outers[j,num] for j,num in enumerate(signature)])
                area=union_area([floats[j,num] for j,num in enumerate(signature)])
                assert inner<=outer
                sigs[signature]={'id':len(sigs),'types':signature,'representative_order':order,
                                 'inner_area_numerator':inner,'outer_area_numerator':outer,
                                 'floating_area':area,'multiplicity':0}
            row=sigs[signature];row['multiplicity']+=1
            records.append((name,''.join(map(str,order)),row['id']))
            count+=1
            if row['floating_area']<best_float:
                best_float=row['floating_area'];best_sig=signature;best_order=order
        orbit_stats[name]={'last':last,'permutations':count,'best_floating_area':best_float,
                           'best_order':best_order,'best_signature_id':sigs[best_sig]['id']}
    rows=sorted(sigs.values(),key=lambda row:row['id'])
    min_inner=min(row['inner_area_numerator'] for row in rows)
    min_outer=min(row['outer_area_numerator'] for row in rows)
    optimal_candidates=[row['id'] for row in rows if row['inner_area_numerator']<=min_outer]
    strict_excess=min_inner-2*D2
    assert all(row['inner_area_numerator']>2*D2 for row in rows), 'proxy<=2 still possible; inspect recorded bounds'
    ledger=OUT/'all-order-signatures.csv'
    with ledger.open('w',newline='') as f:
        writer=csv.writer(f);writer.writerow(('last_orbit','order_vertex_indices','signature_id'));writer.writerows(records)
    output={'status':'PASS: all representative orders have strictly more than 2 proxy area',
        'definition':'At each step choose maximum cardinality old subset S with bbox(new point union S) max span squared < (|S|+1)/8; square is centered at the new point and has area (M_eff+new_mass).',
        'python':sys.version,'platform':platform.platform(),'bits':BITS,
        'arithmetic':'Fraction exact eligibility; integer isqrt radius brackets; exact integer rectangle-union area. Floating values are exploratory only.',
        'input':{'points':[[str(x),str(y)] for x,y in POINTS],'masses':['1/8']*N,'M':'1'},
        'D4_reduction':'All final vertices lie in two D4 orbits: corner 0 or edge midpoint 1. Each has all remaining 7! orders enumerated.',
        'eligibility_checks':eligibility_checks,'permutation_count':len(records),'unique_signature_count':len(rows),
        'observed_square_area_types':sorted({num for sig in sigs for num in sig}),
        'area_denominator':D2,'radius_bounds':bounds,
        'minimum_proxy_area_lower':str(F(min_inner,D2)),'minimum_proxy_area_upper':str(F(min_outer,D2)),
        'minimum_proxy_area_float_interval':[float(F(min_inner,D2)),float(F(min_outer,D2))],
        'strict_excess_over_2_lower':str(F(strict_excess,D2)),
        'possible_optimal_signature_ids':optimal_candidates,'orbit_stats':orbit_stats,'signatures':rows,
        'ledger_sha256':hashlib.sha256(ledger.read_bytes()).hexdigest(),
        'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        'elapsed_seconds':time.perf_counter()-start,
        'scope':'A strict lower bound for an insertion upper-bound proxy, not a lower bound for C2.'}
    result=OUT/'exact-eight-point-proxy.json'
    result.write_text(json.dumps(output,indent=2)+'\n')
    brief={k:v for k,v in output.items() if k not in ('signatures','radius_bounds')}
    print(json.dumps(brief,indent=2))
if __name__=='__main__': main()
