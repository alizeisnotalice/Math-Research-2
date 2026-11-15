"""Exact finite check of all-subset/bbox equivalence and one-axis affine pieces.
Two labels of masses 9,16 at (0,0),(s,0), lambda=1, s in [-1,7].
All three subset radii are rational (3/2,2,5/2), R=3 is a valid tail cap.
"""
from pathlib import Path
from fractions import Fraction as F
from itertools import combinations,product
from math import isqrt
from datetime import datetime,timezone
import argparse,hashlib,json,sys,time

def rectangles(s,saturated=False):
    points=[(F(0),F(0),F(9)),(s,F(0),F(16))];out=[]
    sources=[tuple(points[j]for j in ids)for n in(1,2)for ids in combinations(range(2),n)]
    if saturated:
        xs=sorted({p[0]for p in points});ys=[F(0)];sources=[]
        for a,b in combinations(xs,2):sources.append(tuple(p for p in points if a<=p[0]<=b))
        sources.extend(tuple(p for p in points if p[0]==x)for x in xs)
    for S in sources:
        mass=sum(p[2]for p in S);r2=mass/4;n,d=isqrt(r2.numerator),isqrt(r2.denominator);assert n*n==r2.numerator and d*d==r2.denominator;r=min(F(3),F(n,d))
        x=[p[0]for p in S];y=[p[1]for p in S];L=max(x)-r;U=min(x)+r;B=max(y)-r;T=min(y)+r
        L,U=max(F(-10),L),min(F(11),U);B,T=max(F(-4),B),min(F(4),T)
        if L<U and B<T:out.append((L,U,B,T))
    return out

def area(s,saturated=False):
    rs=rectangles(s,saturated);xs=sorted({x for r in rs for x in r[:2]});ys=sorted({y for r in rs for y in r[2:]});total=F(0)
    for (a,b),(c,d)in product(zip(xs,xs[1:]),zip(ys,ys[1:])):
        x,y=(a+b)/2,(c+d)/2
        if any(L<x<U and B<y<T for L,U,B,T in rs):total+=(b-a)*(d-c)
    return total

def events():
    es={F(-1),F(0),F(7)}
    for l,u in((F(-1),F(0)),(F(0),F(7))):
        # affine (slope,constant) for all fixed-subset x endpoints and F clipping.
        ep=[(0,F(-10)),(0,F(11)),(0,F(-3,2)),(0,F(3,2)),(1,F(-2)),(1,F(2))]
        ep += [(0,F(-5,2)),(1,F(5,2))]if u<=0 else [(1,F(-5,2)),(0,F(5,2))]
        for (a,b),(c,d)in combinations(ep,2):
            if a!=c:
                root=(d-b)/(a-c)
                if l<root<u:es.add(root)
    return sorted(es)

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--output',type=Path,required=True);a=ap.parse_args();start=time.perf_counter();ev=events();pieces=[];tests=[]
    assert F(25)<=4*F(3)**2
    for l,u in zip(ev,ev[1:]):
        p,q=(2*l+u)/3,(l+2*u)/3;beta=(area(q)-area(p))/(q-p);alpha=area(p)-beta*p
        samples=[l+(u-l)*F(k,7)for k in range(8)]
        for s in samples:
            A=area(s);assert A==alpha+beta*s and A==area(s,True);tests.append({'s':str(s),'area':str(A)})
        pieces.append({'interval':[str(l),str(u)],'alpha':str(alpha),'beta':str(beta)})
    gain=area(ev[-1])-area(ev[0]);telescoped=sum(F(z['beta'])*(F(z['interval'][1])-F(z['interval'][0]))for z in pieces);assert gain==telescoped
    # Collision at s=0 is handled by labelled subsets and saturated boxes alike.
    assert area(F(0))==25 and area(F(0),True)==25
    # The joint two-label rectangle dies at s=5, with continuous area at that wall.
    assert any(r[1]-r[0]==F(1,10)for r in rectangles(F(49,10)))and len(rectangles(F(5)))==2
    out={'status':'passed','scope':'one finite two-atom rational example only; not a proof for arbitrary inputs','model':{'weights':['9','16'],'positions':['(0,0)','(s,0)'],'parameter_interval':['-1','7'],'level':'1','R':'3'},'events':list(map(str,ev)),'pieces':pieces,'sample_checks':len(tests),'all_subset_vs_saturated_bbox_areas_equal':True,'persistent_source_y_equals':True,'source_collision_at_s0_area':'25','empty_rectangle_wall_s':'5','area_at_start':str(area(ev[0])),'area_at_end':str(area(ev[-1])),'net_area_change':str(gain),'slope_integral':str(telescoped),'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'python':sys.version,'command':[sys.executable,*sys.argv],'verified_utc':datetime.now(timezone.utc).isoformat(),'elapsed_seconds':time.perf_counter()-start}
    a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
if __name__=='__main__':main()
