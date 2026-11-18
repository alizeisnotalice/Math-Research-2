"""Test the full eleven-point cloud in the hard leaf and in the exact stationary five-body host."""
import argparse,json
from fractions import Fraction as F
from pathlib import Path
from hashlib import sha256
from math import isqrt
import sys
sys.path.insert(0,str(Path(__file__).resolve().parent))
from template import build
from flint import arb,ctx


def q(x):
    x=F(x);return arb(x.numerator)/x.denominator


def posed(points,entry):
    x,y,a=map(q,entry['pose'])
    assert entry['hand'] in (-1,1)
    s,c=a.sin_cos()
    return [(x+c*q(X)-s*entry['hand']*q(Y),y+s*q(X)+c*entry['hand']*q(Y)) for X,Y in points]


def containment(h,D,points):
    minimum=None
    for i,(a,b) in enumerate(h):
        c,d=h[(i+1)%len(h)]
        row=min(((c-a)*(Y*D-b)-(d-b)*(X*D-a)).lower() for X,Y in points)
        minimum=row if minimum is None or row<minimum else minimum
    return minimum/(D*D)


def area_bound(h,D,k):
    S=D-4
    per=sum(1+isqrt((x2-x1)**2+(y2-y1)**2) for (x1,y1),(x2,y2) in zip(h,h[1:]+h[:1]))
    twice=sum(x1*y2-y1*x2 for (x1,y1),(x2,y2) in zip(h,h[1:]+h[:1]))
    e=F(2,S)+F(22,7)**2/(72*k*k)
    return F(twice,2*D*D),F(twice,2*S*S)+e*F(per,S)+F(22,7)*e*e


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for key in ('comp','points11','input','output'):p.add_argument(key,type=Path)
    a=p.parse_args();sys.path.insert(0,str(a.comp.resolve()));assert __debug__
    from certify_nine import append_point,hull
    import certify_five_body_contact_stationary as st
    ctx.prec=192
    extra=json.loads(a.input.read_bytes())
    cloud30=json.loads((a.comp/'finite-point-leaf-input.json').read_bytes())
    p30=[(F(x,cloud30['denominator']),F(y,cloud30['denominator'])) for x,y in cloud30['points']]
    raw=a.points11.read_bytes()
    assert sha256(raw).hexdigest()=='05fb76a261b46c7b939885becb32b31b3e5faca2b017dde5ff70123ff22f6a0d'
    p11=json.loads(raw)['points']
    assert len(p11)==11 and len(p30)==30
    lohost=extra['hard']['center']
    box=[[z,z] for z in lohost]
    assert len(box)==11
    k=4096
    h,D=build(a.comp,False,k,True,box=box)
    low11=containment(h,D,posed(p11,extra['hard']['eleven']))
    assert low11>0
    lo,hi=area_bound(h,D,k)
    assert hi<F('.837458260623')
    # P30 is included at its own separately specified orientation and translation.
    ps=posed(p30,extra['hard']['p30'])
    rounded=list(h)
    for x,y in ps:append_point(rounded,x,y)
    hmix=hull(rounded)
    bot,top=area_bound(hmix,D,k)
    assert top<F('.840')
    # Evaluate the full two-angle isolating box from the global uniqueness theorem.
    ctx.prec=192
    _,actual,_,_=st.model(*(st.interval(*s) for s in st.BOX))
    grid=2**44
    bounds=[[int((v.lower()*grid).floor().unique_fmpz()),int((v.upper()*grid).ceil().unique_fmpz())] for v in actual]
    inside,D=build(a.comp,False,256,True,box=bounds)
    star11=containment(inside,D,posed(p11,extra['star']['eleven']))
    star5=containment(inside,D,posed([p11[i] for i in (0,2,5,7,9)],extra['star']['five']))
    assert star5>0
    result=dict(hard_five_curves_full_eleven_area=[str(lo),str(hi)],hard_eleven_min_edge_determinant=str(low11),
                hard_all_seven_at_listed_poses_area=[str(bot),str(top)],hard_all_seven_upper_ceil_12=str(F((top*10**12).__ceil__(),10**12)),
                star_stationary_parameter_box=st.BOX,star_grid_cover=bounds,star_five_min_determinant=str(star5),
                star_eleven_min_determinant=str(star11),star_contains_all_eleven=bool(star11>0),
                hard_arc_samples=k,star_arc_samples=256, input_sha256=sha256(a.input.read_bytes()).hexdigest(),
                geometry_sha256={n:sha256((a.comp/n).read_bytes()).hexdigest() for n in ('certify_five_body_contact_stationary.py','certify_four_body_contact_stationary.py','five-global-minimum-single-worker-replay.json','finite-point-leaf-input.json')})
    a.output.write_text(json.dumps(result,indent=2)+'\n')
    print(result['star_contains_all_eleven'],result['star_five_min_determinant'],result['hard_all_seven_upper_ceil_12'])
if __name__=='__main__':main()
