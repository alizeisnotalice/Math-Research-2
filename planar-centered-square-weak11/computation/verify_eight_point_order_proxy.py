#!/usr/bin/env python3
"""Independent integer-cell verification of the eight-point proxy certificate."""
from fractions import Fraction as F
from pathlib import Path
import csv,hashlib,json,math,platform,sys,time
ROOT=Path(__file__).resolve().parents[1]
RESULT=ROOT/'computation/results/eight-point-order-proxy/exact-eight-point-proxy.json'
LEDGER=RESULT.with_name('all-order-signatures.csv')

def cell_area(rects):
    xs=sorted({t for r in rects for t in r[:2]})
    ys=sorted({t for r in rects for t in r[2:]})
    area=0
    for x0,x1 in zip(xs,xs[1:]):
        for y0,y1 in zip(ys,ys[1:]):
            if any(a<=x0 and x1<=b and c<=y0 and y1<=d for a,b,c,d in rects):
                area+=(x1-x0)*(y1-y0)
    return area

def exact_type(order,pts):
    # Independent eligibility classification: before final insertion no full-span
    # subset is eligible, and a half-span bbox contains at most 3 perimeter points.
    old=[];out=[None]*8
    for step,j in enumerate(order):
        if step==7:
            num=8
        else:
            near=[i for i in old if max(abs(pts[i][0]-pts[j][0]),abs(pts[i][1]-pts[j][1]))<=F(49,100)]
            num=1
            if near: num=2
            if any(max(abs(pts[i][0]-pts[k][0]),abs(pts[i][1]-pts[k][1]))<=F(49,100)
                   for z,i in enumerate(near) for k in near[z+1:]):
                num=3
        out[j]=num;old.append(j)
    return tuple(out)

def main():
    start=time.perf_counter();r=json.loads(RESULT.read_text());b=r['bits']
    assert b==32
    pts=tuple(tuple(F(v) for v in p) for p in r['input']['points'])
    assert len(pts)==8 and all(m=='1/8' for m in r['input']['masses'])
    assert F(49,100)**2<F(1,4)
    assert F(98,100)**2>F(7,8) and F(98,100)**2<1
    D=100*(1<<b);D2=D*D
    assert r['area_denominator']==D2
    sigs={row['id']:row for row in r['signatures']}
    counts={sid:0 for sid in sigs};orders={'corner':set(),'edge_midpoint':set()}
    assert hashlib.sha256(LEDGER.read_bytes()).hexdigest()==r['ledger_sha256']
    with LEDGER.open() as f:
        for row in csv.DictReader(f):
            name=row['last_orbit'];order=tuple(map(int,row['order_vertex_indices']));sid=int(row['signature_id'])
            assert set(order)==set(range(8)) and len(order)==8
            assert order[-1]==(0 if name=='corner' else 1)
            assert order not in orders[name];orders[name].add(order)
            assert exact_type(order,pts)==tuple(sigs[sid]['types'])
            counts[sid]+=1
    assert all(len(values)==math.factorial(7) for values in orders.values())
    min_inner=None;min_outer=None
    for sid,row in sigs.items():
        assert counts[sid]==row['multiplicity']
        inner=[];outer=[]
        for (x,y),num in zip(pts,row['types']):
            assert num in (1,2,3,8)
            z=num << (2*b-5);low=math.isqrt(z);high=low+(low*low<z)
            assert low*low<=z<=high*high
            xi=int(x*D);yi=int(y*D)
            inner.append((xi-100*low,xi+100*low,yi-100*low,yi+100*low))
            outer.append((xi-100*high,xi+100*high,yi-100*high,yi+100*high))
        lo=cell_area(inner);hi=cell_area(outer)
        assert lo==row['inner_area_numerator'] and hi==row['outer_area_numerator']
        assert lo>2*D2
        min_inner=lo if min_inner is None else min(min_inner,lo)
        min_outer=hi if min_outer is None else min(min_outer,hi)
    assert str(F(min_inner,D2))==r['minimum_proxy_area_lower']
    assert str(F(min_outer,D2))==r['minimum_proxy_area_upper']
    # Independent closed-form scalar and overlap comparisons for winner 237.
    constant=F(19,8)-2*F(13,50)**2-4*F(1,32)-2*F(3,32)-2*F(6,25)**2
    assert constant==F(18121,10000)
    scale=1<<64
    sqrt_bounds={}
    for d in (2,6):
        z=d <<128; lo=math.isqrt(z); hi=lo+(lo*lo<z)
        sqrt_bounds[d]=(F(lo,scale),F(hi,scale))
    rlo,rhi=(q/8 for q in sqrt_bounds[2]); ulo,uhi=(q/8 for q in sqrt_bounds[6])
    a=F(49,100)
    assert a>rhi+uhi and a-rhi>F(1,4)
    assert a-uhi>F(1,2)-a and ulo>F(6,25)
    formula_lo=constant+F(3,25)*sqrt_bounds[6][0]-F(1,200)*sqrt_bounds[2][1]
    formula_hi=constant+F(3,25)*sqrt_bounds[6][1]-F(1,200)*sqrt_bounds[2][0]
    assert F(min_inner,D2)<=formula_lo<=formula_hi<=F(min_outer,D2)
    winner=sigs[237]
    assert all(row['inner_area_numerator']>winner['outer_area_numerator'] for sid,row in sigs.items() if sid!=237)
    output={'status':'PASS','python':sys.version,'platform':platform.platform(),
            'method':'independent local-bbox eligibility classification and integer rectangular-cell coverage; no production union-area helper imported',
            'orders_checked':sum(len(v) for v in orders.values()),'signatures_checked':len(sigs),
            'minimum_proxy_area_lower':str(F(min_inner,D2)),
            'minimum_proxy_area_upper':str(F(min_outer,D2)),
            'closed_form':'18121/10000 + 3*sqrt(6)/25 - sqrt(2)/200',
            'closed_form_interval':[str(formula_lo),str(formula_hi)],
            'proof_sha256':hashlib.sha256((ROOT/'proofs/eight-point-order-proxy-obstruction.md').read_bytes()).hexdigest(),
            'elapsed_seconds':time.perf_counter()-start,
            'result_sha256':hashlib.sha256(RESULT.read_bytes()).hexdigest(),
            'ledger_sha256':hashlib.sha256(LEDGER.read_bytes()).hexdigest(),
            'verifier_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
            'scope':'Finite 8-point proxy certificate only, not a C2 lower bound.'}
    target=RESULT.with_name('independent-cell-verification.json')
    target.write_text(json.dumps(output,indent=2)+'\n');print(json.dumps(output,indent=2))
if __name__=='__main__':main()
