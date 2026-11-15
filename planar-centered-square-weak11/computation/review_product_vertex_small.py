"""Exact two-atom halo check of axis-product partitions and vertex domination.
A finite example only. General product reduction is justified analytically in review.
"""
from pathlib import Path
from fractions import Fraction as F
from itertools import combinations,product
from math import isqrt
from datetime import datetime,timezone
import argparse,hashlib,json,sys,time

def area(x,y):
    pts=[(F(0),F(0),F(9)),(x,y,F(16))];rs=[]
    for ids in ((0,),(1,),(0,1)):
        S=[pts[j]for j in ids];mass=sum(p[2]for p in S);r2=mass/4;n,d=isqrt(r2.numerator),isqrt(r2.denominator);assert n*n==r2.numerator and d*d==r2.denominator;r=min(F(3),F(n,d))
        L=max(p[0]for p in S)-r;U=min(p[0]for p in S)+r;B=max(p[1]for p in S)-r;T=min(p[1]for p in S)+r
        if L<U and B<T:rs.append((L,U,B,T))
    xx=sorted({v for r in rs for v in r[:2]});yy=sorted({v for r in rs for v in r[2:]});ans=F(0)
    for (a,b),(c,d)in product(zip(xx,xx[1:]),zip(yy,yy[1:])):
        m,k=(a+b)/2,(c+d)/2
        if any(L<m<U and B<k<T for L,U,B,T in rs):ans+=(b-a)*(d-c)
    return ans

def full_axis_events():
    # Every candidate endpoint for every possible subset extremal label, plus F.
    e=[(F(0),F(-10)),(F(0),F(11))]
    for ids,r in (((0,),F(3,2)),((1,),F(2)),((0,1),F(5,2))):
        for j in ids:
            e.extend(((F(j),-r),(F(j),r)))
    roots={F(-1),F(0),F(7)}
    for (a,b),(c,d)in combinations(e,2):
        if a!=c:
            z=(d-b)/(a-c)
            if -1<=z<=7:roots.add(z)
    return sorted(roots)

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--output',type=Path,required=True);a=ap.parse_args();start=time.perf_counter();ev=full_axis_events();cache={}
    def A(x,y):
        key=x,y
        if key not in cache:cache[key]=area(x,y)
        return cache[key]
    cells=[];checks=0
    for (l,u),(b,t)in product(zip(ev,ev[1:]),repeat=2):
        corners=[A(l,b),A(u,b),A(l,t),A(u,t)];upper=max(corners)
        for p,q in product((F(0),F(1,3),F(1,2),F(2,3),F(1)),repeat=2):
            x,y=l+p*(u-l),b+q*(t-b);v=A(x,y);interp=(1-p)*(1-q)*corners[0]+p*(1-q)*corners[1]+(1-p)*q*corners[2]+p*q*corners[3]
            assert v==interp and v<=upper;checks+=1
        cells.append({'x':[str(l),str(u)],'y':[str(b),str(t)],'max_corner_area':str(upper)})
    # Low-dimensional product face x=0 retains all collisions on that axis.
    face_checks=0
    for b,t in zip(ev,ev[1:]):
        for q in (F(0),F(1,3),F(1,2),F(2,3),F(1)):
            y=b+q*(t-b);assert A(F(0),y)==(1-q)*A(F(0),b)+q*A(F(0),t);face_checks+=1
    assert A(F(0),F(0))==25
    # A cross-axis coupled feasible line can maximize in its relative interior.
    coupled=[(F(k,8),F(1)-F(k,8)**2)for k in range(-8,9)];assert max(coupled,key=lambda z:z[1])==(F(0),F(1))and coupled[0][1]==coupled[-1][1]==0
    out={'status':'passed','scope':'finite two-atom exact halo samples and polynomial counterexample only; general theorem is analytic','weights':['9','16'],'positions':['(0,0)','(x,y)'],'level':'1','R':'3','domain':['[-1,7]','[-1,7]'],'all_candidate_axis_events':list(map(str,ev)),'closed_product_cells':len(cells),'bilinear_interpolation_and_vertex_domination_checks':checks,'low_dimensional_x0_face_checks':face_checks,'persistent_and_collision_face_covered':True,'unique_area_evaluations':len(cache),'maximum_product_vertex_area':str(max(A(x,y)for x,y in product(ev,repeat=2))),'coupled_example':'f=1+xy, x=s, y=-s on [-1,1], max=1 at s=0; endpoints=0; not claimed to be a halo example','script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'python':sys.version,'command':[sys.executable,*sys.argv],'verified_utc':datetime.now(timezone.utc).isoformat(),'elapsed_seconds':time.perf_counter()-start}
    a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
if __name__=='__main__':main()
