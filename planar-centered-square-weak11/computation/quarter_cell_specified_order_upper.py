#!/usr/bin/env python3
"""A single chosen order: exact outward integer area, not optimality."""
from fractions import Fraction as F
from pathlib import Path
import hashlib,json,math,os,platform,sys,time
ROOT=Path(__file__).resolve().parents[1]
DIR=ROOT/'computation/results/quarter-cell-prefix-group-proxy'
SOURCE=DIR/'quarter-cell-single-input.json'
ORDER=(5,6,7,3,0,4,2,1)
BITS=32

def cell_area(rectangles):
    xs=sorted({v for q in rectangles for v in q[:2]})
    ys=sorted({v for q in rectangles for v in q[2:]})
    area=0
    for x0,x1 in zip(xs,xs[1:]):
        for y0,y1 in zip(ys,ys[1:]):
            if any(a<=x0 and x1<=b and c<=y0 and y1<=d for a,b,c,d in rectangles):
                area+=(x1-x0)*(y1-y0)
    return area

def main():
    start=time.perf_counter();assert os.path.ismount('/Volumes/data')
    source=json.loads(SOURCE.read_text())
    cells=[tuple(F(v) for v in q) for q in source['input']['cells']]
    assert len(cells)==8 and set(ORDER)==set(range(8))
    assert all(F(v)==F(1,8) for v in source['input']['cell_masses'])
    D=4*(1<<BITS);outer=[];rank=[None]*8;radius_brackets=[]
    for k,j in enumerate(ORDER,1):
        rank[j]=k
        z=k << (2*BITS-5);low=math.isqrt(z);high=low+(low*low<z)
        assert low*low<=z<=high*high
        x0,x1,y0,y1=(int(v*D) for v in cells[j]);h=4*high
        outer.append((x0-h,x1+h,y0-h,y1+h))
        radius_brackets.append({'rank':k,'cell':j,'radius_squared':str(F(k,32)),
            'upper_radius_numerator':high,'radius_denominator':1<<BITS})
    upper_num=cell_area(outer);upper=F(upper_num,D*D)
    assert upper<F(7,2)
    output={'status':'PASS: specified order has a rigorous proxy area upper bound below 7/2',
        'python':sys.version,'platform':platform.platform(),
        'definition':'Complete cell bbox expanded by half-side sqrt(cumulative prefix mass)/2; no effective eligibility.',
        'input':source['input'],'specified_order':ORDER,'rank_by_cell':rank,
        'method':'32-bit isqrt outward radii; exact integer rectangular-cell area coverage',
        'coordinate_denominator':D,'area_denominator':D*D,'radius_upper_bounds':radius_brackets,
        'outer_rectangles_integer':outer,'area_upper_numerator':upper_num,
        'area_upper_exact':str(upper),'area_upper_float':float(upper),
        'strict_margin_below_7_2':str(F(7,2)-upper),
        'scope':'Upper bound for this one reproducible order and mass-1 input; no exact global-order optimality or general four-side concentration theorem claimed.',
        'production_result_sha256':hashlib.sha256(SOURCE.read_bytes()).hexdigest(),
        'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        'elapsed_seconds':time.perf_counter()-start}
    out=DIR/'specified-order-outer-certificate.json'
    out.write_text(json.dumps(output,indent=2)+'\n');print(json.dumps(output,indent=2))
if __name__=='__main__':main()
