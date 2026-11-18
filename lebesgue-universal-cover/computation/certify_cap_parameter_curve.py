"""Sharp scalar choice of cap depth and point shrink for the fixed R195 materials."""
from fractions import Fraction as F
from math import isqrt
import json
from pathlib import Path
from flint import ctx

from certify_refined_material_contraction import SCALE, ar, arb, read


def verify(comp):
    ctx.prec=160
    s=SCALE; eps=1-s; r=F(1,2); R=F(".703"); a=F(".431"); radius=F(".528619")
    t=F(500000000377835707,500000000000000000)
    core=F(10**12-5,1000000229098); base=F(".837295")
    old=read(comp/"refined-material-contraction.json")
    assert old["point_factor"]==str(t/(1-F(".000000757")))
    def square(T):
        return (T-s*s*base)*R/(s*r)
    assert 0<square(F(".83742"))<square(F(".837531"))<F(".019")**2
    assert (1+s)*r*r>2*a*(r+F(".019"))
    assert 1+s>(2*a*r+(1+s)*r*r+eps*a*a)/(4*r*r)
    assert 8*r*s*s-2*(r+F(".019"))>0
    assert 8*r*a*s*s-(1+s)*F(".019")*(2*r+F(".019"))-(1+s)*r*r>0
    scale=10**12
    samples=[]
    for i in range(83742,83754):
        G=F((F(i,100000)*2**57).__floor__()+1,2**57)
        exact_h=ar(square(G)).sqrt()
        def root(h):
            length=2*ar(r)+h
            B=4*ar(s*s)*length-2*(ar(r)+h)
            C=4*ar(a*s*s)*length-ar(1+s)*(h*length+ar(r*r))
            return 2*C/(B+(B*B+4*ar(eps)*C).sqrt())
        opt=root(exact_h)*ar(eps/radius)
        # Round the height up and the point shrink down; prove strict comparison rationally.
        sq=square(G)
        h=F(isqrt(sq.numerator*scale**2//sq.denominator)+1,scale)
        rho=F(int((root(ar(h))*ar(eps/radius)*scale).lower().floor().unique_fmpz()),scale)
        m=rho*radius/eps;v=(1+s)*r*r+eps*m*m
        slack=4*(a-m)*s*s-(1+s)*h-(v-2*m*(r+h))/(2*r+h)
        assert m<a and h*h>sq and slack>0 and 0<opt-ar(rho)<ar(F(2,scale))
        factor=min(core/s,t/(1-rho))
        ceiling=(ar(core/s)).min(ar(t)/(1-opt))**2*ar(G)
        area=G*factor**2
        assert ceiling-ar(area)<ar(F(".000000000004"))
        samples.append(dict(target=str(F(i,100000)),integer_goal=str(G),h=str(h),rho=str(rho),
                            supremum_point_shrink=str(opt),arithmetic_slack=str(slack),
                            conditional_floor=str(area),conditional_display=float(area),
                            scalar_area_supremum=str(ceiling)))
    return dict(curve_scale=str(s),radius=str(radius),samples=samples)


if __name__=="__main__":
    print(json.dumps(verify(Path(__file__).resolve().parent),indent=2))
