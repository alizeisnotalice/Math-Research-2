"""Upper witnesses for each single deletion at the exact five-body minimizer."""
import argparse,json,sys
from fractions import Fraction as F
from math import isqrt
from pathlib import Path
from flint import arb,ctx
from template import build
from hinge_bank import q

def run(c,root):
    ctx.prec=192;sys.path.insert(0,str(c.resolve()))
    from certify_nine import append_point,hull
    import certify_five_body_contact_stationary as st
    _,actual,_,_=st.model(*(st.interval(*s) for s in st.BOX))
    p=[int((x.mid()*2**44).floor().unique_fmpz()) for x in actual]
    assert all(abs(x-q(F(v,2**44)))<q(F(2,2**44)) for x,v in zip(actual,p))
    k=2048;h,D=build(c,False,k,True,box=[[z,z] for z in p]);S=D-4
    raw=json.loads((root/'new-points.json').read_bytes())
    cloud=[(q(F(x,raw['denominator'])),q(F(y,raw['denominator']))) for x,y in raw['points']]
    inside,scale=build(c,False,192,True,box=[[int((x.lower()*2**44).floor().unique_fmpz()),int((x.upper()*2**44).ceil().unique_fmpz())] for x in actual])
    entries=json.loads((root/'deletions.json').read_bytes())
    assert sorted(row['deleted'] for row in entries)==list(range(11))
    result={}
    for row in entries:
        ang=arb.pi()*row['angle_numerator']/row['angle_denominator']
        sn,cs=ang.sin_cos();hand=row['hand'];assert hand in (1,-1)
        tx=q(F(row['tx'],row['translation_grid']));ty=q(F(row['ty'],row['translation_grid']))
        pts=list(h);placed=[]
        for j,(x,y) in enumerate(cloud):
            if j!=row['deleted']:
                X,Y=tx+cs*x-hand*sn*y,ty+sn*x+hand*cs*y
                append_point(pts,X,Y);placed.append((X,Y))
        z=hull(pts);edges=list(zip(z,z[1:]+z[:1]));twice=sum(x*Y-y*X for (x,y),(X,Y) in edges)
        per=sum(1+isqrt((X-x)**2+(Y-y)**2) for (x,y),(X,Y) in edges)
        # All host points differ from the snapped component positions by less
        # than 6/2^44: translation norm <=sqrt(2)*2/2^44 and rotation costs
        # at most .58 * 2/2^44. The disk and posed rational points are fixed.
        eps=F(2,S)+F(22,7)**2/(72*k*k)+F(6,2**44)
        hi=F(twice,2*S*S)+eps*F(per,S)+F(22,7)*eps**2
        assert hi<F('.83750'),(row['deleted'],hi)
        min_edge=min((((X-x)*(Yp*scale-y)-(Y-y)*(Xp*scale-x))/scale**2).lower()
                     for (x,y),(X,Y) in zip(inside,inside[1:]+inside[:1]) for Xp,Yp in placed)
        assert min_edge>0
        result[row['deleted']]=dict(area_upper=str(F((hi*10**12).__ceil__(),10**12)),minimum_edge_determinant=str(min_edge))
    return result

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('comp',type=Path);p.add_argument('output',type=Path)
    a=p.parse_args();assert __debug__;v=run(a.comp,Path(__file__).resolve().parent)
    print(max(F(z['area_upper']) for z in v.values()),flush=True);a.output.write_text(json.dumps(v,indent=2)+'\n')
