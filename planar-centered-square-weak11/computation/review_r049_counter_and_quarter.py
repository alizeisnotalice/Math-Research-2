from fractions import Fraction as F
from itertools import combinations
from pathlib import Path
import argparse,hashlib,json,platform,time

def span(pts):
    return max(max(p[k] for p in pts)-min(p[k] for p in pts) for k in (0,1)) if pts else F(0)
def union_area_vertical_sweep(rects):
    xs=sorted({q for a,b,c,d in rects for q in (a,b)});answer=0
    for lo,hi in zip(xs,xs[1:]):
        intervals=sorted((c,d) for a,b,c,d in rects if a<=lo and hi<=b)
        length=0;end=None
        for c,d in intervals:
            if end is None or c>end:length+=d-c;end=d
            elif d>end:length+=d-end;end=d
        answer+=(hi-lo)*length
    return answer

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--source-root',type=Path,required=True);args=ap.parse_args()
    root=args.source_root;start=time.perf_counter()
    paths=['proofs/post-r40-sorting-obstructions.md','proofs/quarter-cell-prefix-group-single-test.md','computation/post_r40_sorting_obstruction_check.py','computation/results/post-r40-sorting/exact-obstructions.json','computation/results/post-r40-sorting/HANDOFF-MANIFEST.json','computation/quarter_cell_specified_order_upper.py','computation/results/quarter-cell-prefix-group-proxy/quarter-cell-single-input.json','computation/results/quarter-cell-prefix-group-proxy/specified-order-outer-certificate.json','computation/results/quarter-cell-prefix-group-proxy/HANDOFF-MANIFEST.json']
    sources={s:dict(bytes=(root/s).stat().st_size,sha256=hashlib.sha256((root/s).read_bytes()).hexdigest()) for s in paths}
    assert sources[paths[0]]['sha256']=='47f0ac8ff9675a02458ef313654f7f44b0153fb7fd55e4cd1d340313f1e08beb'
    assert sources[paths[1]]['sha256']=='8faa60748b4215b8bd8662c3047e31da7ddaf106bd9801d358e53fec17f1f510'
    qdir=root/'computation/results/quarter-cell-prefix-group-proxy'
    inp=json.loads((qdir/'quarter-cell-single-input.json').read_text());cert=json.loads((qdir/'specified-order-outer-certificate.json').read_text())
    assert inp['input']==cert['input']
    assert all(F(m)==F(1,8) for m in inp['input']['cell_masses']) and F(inp['input']['total_mass'])==1
    cells=[tuple(F(v) for v in q) for q in inp['input']['cells']]
    expected=[(x0,x1,y0,y1) for x0,x1 in [(F(0),F(1,4)),(F(1,4),F(3,4)),(F(3,4),F(1))] for y0,y1 in [(F(0),F(1,4)),(F(1,4),F(3,4)),(F(3,4),F(1))] if not (x0==F(1,4) and y0==F(1,4))]
    assert cells==expected
    order=[5,6,7,3,0,4,2,1];assert cert['specified_order']==order
    den=1<<32;D=4*den;assert cert['coordinate_denominator']==D and cert['area_denominator']==D*D
    rects=[]
    for k,i in enumerate(order,1):
        r=cert['radius_upper_bounds'][k-1];H=r['upper_radius_numerator']
        assert r['rank']==k and r['cell']==i and r['radius_denominator']==den and F(r['radius_squared'])==F(k,32)
        target=k*(den*den)//32
        assert H*H>=target and (H-1)*(H-1)<target
        coords=[v*D for v in cells[i]];assert all(v.denominator==1 for v in coords)
        x0,x1,y0,y1=map(int,coords);h=4*H
        rects.append([x0-h,x1+h,y0-h,y1+h])
    assert rects==cert['outer_rectangles_integer']
    area_integer=union_area_vertical_sweep(rects)
    assert area_integer==943119706436290073136==cert['area_upper_numerator']
    upper=F(area_integer,D*D);assert upper==F(cert['area_upper_exact'])==F(58944981652268129571,18446744073709551616)<F(7,2)
    margin=F(7,2)-upper;assert margin==F(cert['strict_margin_below_7_2'])
    # Independent exact scalar/subset reread of the post-R40 mathematical examples.
    e=F(1,10000);a=F(99,100);light=[(x,y) for x in (-a,F(0),a) for y in (-a,F(0),a) if (x,y)!=(F(0),F(0))]
    four_count=0
    for S in combinations(range(8),4):assert span([light[i] for i in S])==2*a;four_count+=1
    assert F(197,50)/(1+8*e)==F(4925,1251)
    assert (1+11*e)/(1+8*e)==F(3337,3336)
    assert F(199,100)>1+11*e
    assert F(1001,1000)**2>1+3*e and F(1001,1000)+F(1,100)<2*a
    qmax=F(99,100000)
    assert qmax*qmax<1 and (a+qmax)**2<1+F(1,50) and qmax<F(1,100)
    assert F(22201,2500)/F(29,25)==F(22201,2900)
    assert F(103,100)**2>F(53,50)
    assert F(4,25)+F(53,50)+F(1,1000)*F(103,100)==F(122103,100000)
    assert F(122103,100000)/F(29,25)==F(122103,116000)
    aa=F(49,100);cycle=[(aa,F(0)),(aa/2,aa),(-aa/2,aa),(-aa,F(0)),(-aa/2,-aa),(aa/2,-aa)]
    pair_edges=[]
    for i,j in combinations(range(6),2):
        valid=span([cycle[i],cycle[j]])**2<F(1,3);assert valid==(j-i in (1,5))
        if valid:pair_edges.append([i,j])
    proper_count=0
    for k in (3,4,5):
        for S in combinations(range(6),k):assert span([cycle[i] for i in S])**2>=F(k,6);proper_count+=1
    for gone in range(6):assert span([p for i,p in enumerate(cycle) if i!=gone])**2==F(2401,2500)>F(5,6)
    assert 9*9>4*4*2 # q-proxy contains p-proxy: 9>4sqrt(2).
    assert 15-2==13 and 1+2*F(5,6)==F(8,3)
    out=dict(status='PASS',python=platform.python_version(),method='Independent integer vertical-sweep union of the one specified outward certificate; no source union helper or floating ordering search imported',specified_order=order,radius_ceiling_checks=8,outer_rectangles_checked=8,integer_area=area_integer,coordinate_denominator=D,area_upper=str(upper),strict_margin_below_7_2=str(margin),four_light_bbox_checks=four_count,cycle_pair_checks=15,cycle_proper_higher_subset_checks=proper_count,cycle_deleted_bbox_checks=6,scope='Fixed eight labelled perimeter blocks each mass1/8; covers arbitrary sources assigned once inside each closed block bbox. Not general nine-cell mass, free-order optimality, or C2.',not_rerun=['10080 floating candidate orders','production 2D endpoint-cell union','16-sector production checker','any large search'],source_files=sources,elapsed_seconds=time.perf_counter()-start,review_script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
    print(json.dumps(out,indent=2))
if __name__=='__main__':main()
