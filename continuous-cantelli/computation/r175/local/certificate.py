#!/usr/bin/env python3
"""Directed interval signs for one exact dyadic two-peak graph."""
import json
from fractions import Fraction as F
from pathlib import Path

import mpmath as mp

iv=mp.iv
iv.dps=50
N=512
S=F(5,4)
PEAKS=((F(0),F(3,32),F(3,4),F(1,8),F(7,8)),
       (F(27,32),F(5,32),F(1,2),F(1,8),F(3,4)))


def exact(x):
    return iv.mpf(x.numerator)/iv.mpf(x.denominator)


def hull(a,b):
    return iv.mpf([exact(a).a,exact(b).b])


def density(x,v):
    return iv.exp(-x**2/(2*v))/iv.sqrt(2*iv.pi*v)


def integrate(a,b,count,fn):
    step=(b-a)/count
    w=exact(step)
    out=None
    for i in range(count):
        row=[w*z for z in fn(hull(a+i*step,a+(i+1)*step))]
        out=row if out is None else [p+q for p,q in zip(out,row)]
    return out


def target(u,L):
    p=density(L*(u-iv.mpf("0.5")),exact(S))
    return p,p*u,p*u*u,p*iv.sin(iv.pi*u)


def local(t,left,width,h,t0,t1,s0,s1,L):
    u=exact(left+width*t0)+exact(width*(t1-t0))*t
    v=exact(h*s0)+exact(h*(s1-s0))*t
    p=density(L*(u-iv.mpf("0.5")),1)*exact(width*(t1-t0))
    return p,p*u,p*(u*u-v),p*iv.sin(iv.pi*u)*iv.exp(iv.pi**2*v/2)


def deficits(L):
    R=integrate(F(0),F(1),N,lambda u:target(u,L))
    P=[iv.mpf(0)]*4
    for left,w,h,lo,hi in PEAKS:
        for t0,t1,s0,s1 in ((F(0),lo,F(0),F(1)),(lo,hi,F(1),F(1)),(hi,F(1),F(1),F(0))):
            vals=integrate(F(0),F(1),N,lambda t:local(t,left,w,h,t0,t1,s0,s1,L))
            P=[p+q for p,q in zip(P,vals)]
    return [r-p for r,p in zip(R,P)]


def main():
    d0,d1,d2,ds=deficits(iv.mpf(1))
    para=d1-d2
    det=d0*d2-d1*d1
    physical_x2=d2-d1+d0/4
    assert d0.a>iv.mpf("0.25")
    assert det.a>iv.mpf("0.013")
    assert para.a>iv.mpf("0.005")
    assert physical_x2.a>iv.mpf("0.05")
    assert ds.b<iv.mpf("-0.07")

    z0,z1,z2,zs=deficits(hull(F(0),F(1,4)))
    scaled_det=z0*z2-z1*z1
    scaled_para=z1-z2
    assert z0.a>iv.mpf("0.25")
    assert scaled_det.a>iv.mpf("0.013")
    assert scaled_para.a>iv.mpf("0.001")
    assert zs.b<iv.mpf("-0.07")
    # Optional mean-one completion: B equals 1 on the hull and decays
    # to zero at |x|=3/4; D equals 1 up to |x|=2 and decays to 0 at 3.
    eB=2*(integrate(F(0),F(1,2),N,lambda x:(density(x,1),))[0]
          +integrate(F(1,2),F(3,4),N,lambda x:((3-4*x)*density(x,1),))[0])
    denom_min=2*integrate(F(3,4),F(2),N,lambda x:(density(x,1),))[0]
    eH_max=exact(F(19,128))/iv.sqrt(2*iv.pi)
    beta=iv.mpf(1)/64
    alpha_upper=(eB/4+eH_max+beta)/denom_min
    assert alpha_upper.b<iv.mpf("0.55")

    report=dict(s=str(S),a="-1/2",b="1/2",normalized_peak_parameters=[list(map(str,p)) for p in PEAKS],
                subboxes_per_segment=N,mpmath=mp.__version__,interval_digits=iv.dps,
                D0=str(d0),D1=str(d1),D2=str(d2),positive_quadratic_det=str(det),
                concave_parabola_slack=str(para),physical_x_squared_slack=str(physical_x2),
                sine_slack_RHS_minus_LHS=str(ds),mean_completion_alpha_upper=str(alpha_upper),
                uniform_scaled_width_family=dict(L="0<L<=1/4",D0=str(z0),
                    positive_quadratic_det=str(scaled_det),concave_parabola_slack=str(scaled_para),
                    sine_slack_RHS_minus_LHS=str(zs)))
    Path("certificate.json").write_text(json.dumps(report,indent=2)+"\n")
    for k,v in report.items():
        print(k,v)


if __name__=="__main__":
    main()
