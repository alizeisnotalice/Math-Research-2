#!/usr/bin/env python3
"""One bounded quarter-cell cumulative-prefix group-proxy audit."""
from fractions import Fraction as F
from itertools import permutations,product
from pathlib import Path
import csv,hashlib,json,math,os,platform,sys,time
ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'computation/results/quarter-cell-prefix-group-proxy'
BANDS=((F(0),F(1,4)),(F(1,4),F(3,4)),(F(3,4),F(1)))
CELLS=tuple((BANDS[i][0],BANDS[i][1],BANDS[j][0],BANDS[j][1])
            for i,j in product(range(3),repeat=2) if (i,j)!=(1,1))
BITS=32

def area(rectangles):
    xs=sorted({t for q in rectangles for t in q[:2]})
    total=0
    for lo,hi in zip(xs,xs[1:]):
        intervals=sorted((y0,y1) for x0,x1,y0,y1 in rectangles if x0<=lo and hi<=x1)
        length=0;start=None;end=None
        for a,b in intervals:
            if start is None:start,end=a,b
            elif a>end:length+=end-start;start,end=a,b
            else:end=max(end,b)
        if start is not None:length+=end-start
        total+=(hi-lo)*length
    return total

def main():
    start=time.perf_counter();assert os.path.ismount('/Volumes/data')
    OUT.mkdir(parents=True,exist_ok=True)
    float_rects={}
    for j,cell in enumerate(CELLS):
        for k in range(1,9):
            r=math.sqrt(k/8)/2
            x0,x1,y0,y1=map(float,cell)
            float_rects[j,k]=(x0-r,x1+r,y0-r,y1+r)
    records=[];orbits={}
    for name,last in (('corner',0),('edge_midpoint',1)):
        best=math.inf;best_order=None
        for before in permutations(j for j in range(8) if j!=last):
            order=before+(last,);rank=[None]*8
            for k,j in enumerate(order,1):rank[j]=k
            value=area([float_rects[j,k] for j,k in enumerate(rank)])
            row={'orbit':name,'order':order,'rank':tuple(rank),'float_area':value}
            records.append(row)
            if value<best:best=value;best_order=order
        orbits[name]={'last':last,'count':math.factorial(7),'best_float_area':best,'best_order':best_order}
    best=min(records,key=lambda row:row['float_area'])
    output={'definition':'For each cell j, expand its complete bbox by half-side sqrt(k/8)/2 where k is its insertion rank; union all eight expanded rectangles. No actual common-capture eligibility is used.',
        'input':{'bands':[[str(v) for v in b] for b in BANDS],
                 'cells':[[str(v) for v in cell] for cell in CELLS],
                 'cell_masses':['1/8']*8,'central_cell_mass':'0','total_mass':'1'},
        'python':sys.version,'platform':platform.platform(),'case_count':1,'orders_checked':len(records),
        'D4_reduction':'Two final-cell orbits: corner index0, edge midpoint index1, every remaining 7! ordering; covers all 8! orders by symmetry.',
        'orbits':orbits,'best_floating_record':best,'certificate':None,
        'scope':'Coarse cumulative-prefix mass plus quarter-grid bbox grouping only; distinct from original effective insertion U_pi and from C2.'}
    if best['float_area']>3.500001:
        D=4*(1<<BITS);D2=D*D;inner={};outer={};bounds={}
        for k in range(1,9):
            z=k << (2*BITS-5);lo=math.isqrt(z);hi=lo+(lo*lo<z)
            assert lo*lo<=z<=hi*hi
            bounds[str(k)]={'radius_squared':str(F(k,32)),'lower_numerator':lo,'upper_numerator':hi,'denominator':1<<BITS}
            for j,cell in enumerate(CELLS):
                x0,x1,y0,y1=(int(v*D) for v in cell)
                inner[j,k]=(x0-4*lo,x1+4*lo,y0-4*lo,y1+4*lo)
                outer[j,k]=(x0-4*hi,x1+4*hi,y0-4*hi,y1+4*hi)
        for row in records:
            rank=row['rank']
            row['inner_area_numerator']=area([inner[j,k] for j,k in enumerate(rank)])
            row['outer_area_numerator']=area([outer[j,k] for j,k in enumerate(rank)])
            assert row['inner_area_numerator']<=row['outer_area_numerator']
        lower=min(row['inner_area_numerator'] for row in records)
        upper=min(row['outer_area_numerator'] for row in records)
        assert 2*lower>7*D2
        output['certificate']={'status':'PASS: every representative order has area strictly greater than 7/2',
            'method':'32-bit isqrt inner/outer radii and exact integer rectangle union',
            'bits':BITS,'coordinate_denominator':D,'area_denominator':D2,'radius_bounds':bounds,
            'minimum_lower':str(F(lower,D2)),'minimum_upper':str(F(upper,D2)),
            'minimum_float_interval':[float(F(lower,D2)),float(F(upper,D2))],
            'excess_over_7_2_lower':str(F(lower,D2)-F(7,2))}
    ledger=OUT/'all-orders.csv'
    with ledger.open('w',newline='') as f:
        names=('orbit','order','rank_by_cell','floating_area','inner_area_numerator','outer_area_numerator')
        writer=csv.writer(f);writer.writerow(names)
        for row in records:
            writer.writerow((row['orbit'],''.join(map(str,row['order'])),''.join(map(str,row['rank'])),repr(row['float_area']),
                             row.get('inner_area_numerator',''),row.get('outer_area_numerator','')))
    output['ledger_sha256']=hashlib.sha256(ledger.read_bytes()).hexdigest()
    output['script_sha256']=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    output['elapsed_seconds']=time.perf_counter()-start
    result=OUT/'quarter-cell-single-input.json'
    result.write_text(json.dumps(output,indent=2)+'\n')
    print(json.dumps(output,indent=2))
if __name__=='__main__':main()
