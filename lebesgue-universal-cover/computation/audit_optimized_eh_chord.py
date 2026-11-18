#!/usr/bin/env python3
"""Independent Arb audit of the signed Jensen overlap test at delta=0.0000154."""
from fractions import Fraction as Q
import argparse
import hashlib
import json
from pathlib import Path
import sys

import flint
from flint import arb, ctx


ctx.prec=160
HALF=arb(1)/2
ONE=arb(1)
PI=arb.pi()
SIGMA=arb("1.5494009336023418793")*PI/180


def exact(x):
    return arb(x.numerator)/x.denominator


def unit(theta):
    s,c=theta.sin_cos()
    return c,s


def plus(x,y):
    return x[0]+y[0],x[1]+y[1]


def neg(x):
    return -x[0],-x[1]


def minus(x,y):
    return plus(x,neg(y))


def scale(s,x):
    return s*x[0],s*x[1]


def rot(x):
    return -x[1],x[0]


def dot(x,y):
    return x[0]*y[0]+x[1]*y[1]


def cross(x,y):
    return x[0]*y[1]-x[1]*y[0]


def norm(x):
    return dot(x,x).sqrt()


def line_pair(u,v,alpha=HALF,beta=HALF):
    det=cross(u,v)
    return (alpha*v[1]-beta*u[1])/det,(beta*u[0]-alpha*v[0])/det


def circle_line(center,normal,sign):
    offset=dot(center,normal)-HALF
    return plus(minus(center,scale(offset,normal)),scale(sign*(1-offset*offset).sqrt(),rot(normal)))


def circle_pair(a,b,sign):
    displacement=minus(b,a)
    length=norm(displacement)
    offset=(1-length*length/4).sqrt()
    return plus(scale(HALF,plus(a,b)),scale(sign*offset/length,rot(displacement)))


def upper_exact(a):
    mantissa,exponent=map(int,a.upper().man_exp())
    return arb((mantissa,exponent))


def setup():
    n=[unit(PI*j/3) for j in range(6)]
    m=[unit(PI/6+SIGMA+PI*j/3) for j in range(6)]
    c3=line_pair(n[0],m[5])
    b3=line_pair(n[1],m[0])
    p=unit(7*PI/6)
    u=circle_line(c3,p,-1)
    s=circle_line(line_pair(m[2],n[1],HALF,arb(0)),n[0],-1)
    t=circle_pair(s,b3,1)
    j=circle_pair(b3,c3,-1)
    assert u[1]>arb("-.2") and t[0]<0 and j[0]<0
    return n,m,c3,b3,p,u,s,t,j


def trial(delta_text, tight=False):
    n,m,c3,b3,p,u,s,t,j=setup()
    delta=arb(delta_text)
    c,sine=unit(delta)
    rad=minus(u,c3)
    U=plus(c3,minus(scale(c,rad),scale(sine,rot(rad))))
    chord=minus(t,U)
    assert chord[0]>0 and chord[1]<0
    k=chord[0]/(-chord[1])
    assert 0<k<arb(13)/200 and U[0]<u[0] and 0<U[1]-u[1]<delta
    eta=HALF-dot(p,U)
    gamma=HALF-arb(3).sqrt()*k/2
    assert eta>0 and gamma>0
    eprime=unit(7*PI/6-SIGMA)
    assert U[1]<s[1] and t[1]<s[1] and t[0]<s[0]
    assert norm(minus(U,s))>1
    assert dot(minus(t,s),minus(U,t))>0
    assert dot(eprime,U)>HALF and dot(eprime,t)>HALF
    U_rad,Jc=minus(U,c3),minus(j,c3)
    Jb,Tb=minus(j,b3),minus(t,b3)
    for v,w in (((-ONE,arb(0)),U_rad),(U_rad,rad)):
        assert cross(v,w)>0 and dot(v,w)>0
    for v,w in ((U_rad,Jc),(Jb,Tb)):
        assert cross(v,w)>0 and dot(v,w)>0
        assert all(x<0 for x in v+w)
    assert all(dot(chord,v)<0 for v in (U_rad,Jc))
    assert all(dot(chord,v)>0 for v in (Jb,Tb))
    assert cross(rot(U_rad),chord)>0 and cross(chord,rot(Tb))>0
    assert dot(rot(U_rad),chord)>0 and dot(chord,rot(Tb))>0
    assert all(dot(eprime,rot(v))>0 for v in (U_rad,Jc))
    assert all(dot(eprime,rot(v))<0 for v in (Jb,Tb))
    assert all(dot(minus(c3,s),rot(v))>0 for v in (U_rad,Jc))
    assert all(dot(minus(b3,s),rot(v))<0 for v in (Jb,Tb))

    a=c3[1]-U[1]
    H=HALF-U[0]
    c=H*k
    Y=upper_exact(eta/gamma)
    assert c>a>0 and H>2*k*Y
    coefficient=-k*k if tight else ONE
    def f(y):
        return (a*a+2*c*y+coefficient*y*y).sqrt()-a-y

    # The direct squared-distance estimate tightens the looser source
    # estimate. Both are concave; for the tighter one use f'(2Y)>0.
    if tight:
        ymax=2*Y
        assert (c+coefficient*ymax)/(a*a+2*c*ymax+coefficient*ymax*ymax).sqrt()-1 >0
    R=upper_exact(f(Y))
    rind=upper_exact(f(2*Y))
    assert rind>R>0
    m5=m[5]
    M=-m5[1]
    assert 0<M*rind<1
    A=minus(c3,m5)
    N0=scale(HALF,n[5])
    horizontal=rot(m5)
    excursion=(2*M*R-(M*R)**2).sqrt()
    xi=minus(plus(minus(A,N0),(arb(0),R)),
             minus(scale(excursion,horizontal),scale(M*R,m5)))
    length=norm(xi)
    v=scale(1/length,xi)
    vm=dot(v,m5)
    beta=dot(v,horizontal)
    lambda_signed=v[1]+M*vm
    assert vm<0 and beta>0 and lambda_signed>0
    derivative=beta*M*(1-M*rind)/(2*M*rind).sqrt()-lambda_signed
    assert derivative>0
    exact_derivative=beta*M*(1-M*rind)/(2*M*rind-(M*rind)**2).sqrt()-lambda_signed
    assert exact_derivative>derivative
    gap=length-1
    gain=(delta+cross(minus(c3,t),minus(u,U)))/2
    assert gain>0
    return {"delta":delta_text,"tight_squared_distance":tight,
            "major_axis_projection":str(vm),"tangent_projection":str(beta),
            "signed_linear_coefficient":str(lambda_signed),
            "derivative_lower":str(derivative),"endpoint_derivative":str(exact_derivative),
            "mean_height_upper":str(Y),"mean_CB_upper":str(R),
            "individual_CB_upper":str(rind),"strict_overlap_margin":str(gap),
            "gain_area":str(gain),"passes":bool(gap>0)}


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output",type=Path)
    parser.add_argument("--old-area",type=Path)
    args=parser.parse_args()
    main_trial=trial("0.0000154")
    next_trial=trial("0.0000155")
    refined_trial=trial("0.0000154059",tight=True)
    further_trial=trial("0.000015406",tight=True)
    assert main_trial["passes"] and not next_trial["passes"]
    assert refined_trial["passes"] and not further_trial["passes"]
    out={"precision":ctx.prec,"python":sys.version,"python_flint":flint.__version__,
         "scope":"Independently constructs U,T,J coordinates and certifies the signed Jensen projection calculation.",
         "trials":[main_trial,next_trial,refined_trial,further_trial],
         "script_sha256":hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    if args.old_area:
        old_json=json.loads(args.old_area.read_bytes())
        for row in (main_trial,refined_trial):
            gain=arb(row["gain_area"])
            area=exact(Q(old_json["area_interval"][1]))-gain
            row["candidate_area_upper_from_archived_baseline"]=str(area.upper())
    args.output.write_text(json.dumps(out,indent=2)+"\n")
    print(json.dumps(out,indent=2))


if __name__=="__main__":
    main()
