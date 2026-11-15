#!/usr/bin/env python3
"""R49: independent strip-boundary accounting and one exact block cover."""
from pathlib import Path
from fractions import Fraction as F
import argparse,hashlib,json,math,platform,sys,time

def union_area(rects):
    xs=sorted({x for r in rects for x in r[:2]});out=0
    for a,b in zip(xs,xs[1:]):
      active=sorted((c,d) for x,y,c,d in rects if x<=a and b<=y)
      length=0;top=None
      for c,d in active:
        if top is None or c>top:length+=d-c;top=d
        elif d>top:length+=d-top;top=d
      out+=(b-a)*length
    return out

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--certificate',type=Path,required=True);ap.add_argument('--out',type=Path,required=True);args=ap.parse_args();start=time.perf_counter()
    record=json.loads(args.certificate.read_text());order=record['specified_order'];assert order==[5,6,7,3,0,4,2,1]
    cells=[tuple(map(F,r)) for r in record['input']['cells']];assert len(cells)==8
    b=32;D=1<<b;rects=[]
    for k,j in enumerate(order,1):
      q=k<<(2*b-5);lo=math.isqrt(q);hi=lo+(lo*lo<q)
      assert (hi-1)*(hi-1)<q<=hi*hi
      a,z,c,d=cells[j];ints=[int(v*D) for v in cells[j]];assert all(F(v,D)==v0 for v,v0 in zip(ints,cells[j]))
      rects.append((ints[0]-hi,ints[1]+hi,ints[2]-hi,ints[3]+hi))
    area=F(union_area(rects),D*D);assert str(area)==record['area_upper_exact'];assert area<F(7,2)
    assert [[4*x for x in r] for r in rects]==record['outer_rectangles_integer']
    checks=0
    for a in [F(0),F(1,4),F(41,200),F(1,2)]:
      vals=sorted({F(0),a,F(1,2),1-a,F(1)})
      for x in vals:
       for y in vals:
        edges=[x<a,x>1-a,y<a,y>1-a]
        horizontal=int(edges[0])+int(edges[1]);vertical=int(edges[2])+int(edges[3])
        assert horizontal in [0,1] and vertical in [0,1]
        c=int(horizontal==1 and vertical==1);z=int(horizontal==0 and vertical==0)
        assert sum(edges)==1+c-z;checks+=1
    k=F(359,100);a=(4-k)/2;rho=(k-3)**2
    assert a==F(41,200) and rho==F(3481,10000)
    assert 4*rho-1==F(981,2500)
    assert 4-2*a==k and a+(k-3)/2==F(1,2)
    assert F(197,50)/F(1251,1250)==F(4925,1251)
    assert F(22201,2500)/F(29,25)==F(22201,2900)
    output=dict(status='PASS',python=sys.version,platform=platform.platform(),algorithm='Independent integer x-sweep of specified outward block rectangles; exact indicator checks including all cut boundaries; rational mass-budget identities.',specified_order=order,area_upper_exact=str(area),margin_to_7_2=str(F(7,2)-area),boundary_indicator_checks=checks,K_359_strip_width=str(a),K_359_mass_threshold=str(rho),K_359_corner_minus_center_necessary_lower=str(4*rho-1),arithmetic='integers and fractions.Fraction; no random input or tolerance',elapsed_seconds=time.perf_counter()-start,source_sha256={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in [args.certificate,Path(__file__).resolve()]},scope='One prescribed block-mass distribution and one block order; not all mass-simplex inputs or a general C2 bound. General indicator identity is proved analytically in the review.')
    args.out.parent.mkdir(parents=True,exist_ok=True);args.out.write_text(json.dumps(output,indent=2)+'\n');print(json.dumps(output,indent=2))
if __name__=='__main__':main()
