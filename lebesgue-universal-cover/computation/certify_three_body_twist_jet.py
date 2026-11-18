"""Second-order zero-twist certificate for the three-body contact manifold."""
from hashlib import sha256
import json
from pathlib import Path
import platform

import flint
from flint import arb, arb_series, ctx
from certify_three_body_symmetric_slice import R3, add, c36, c72, det, dot, g, interval, pi, q, r, scale, sub

ctx.cap = 3


def J(x):
    return -x[1], x[0]


def rotate(x, angle):
    return add(scale(angle.cos(), x), scale(angle.sin(), J(x)))


def normal(angle):
    return angle.cos(), angle.sin()


def tangent_normal(p, sign):
    radius2 = dot(p,p)
    return add(scale(r/radius2,p), scale(sign*(radius2-r**2).sqrt()/radius2,J(p)))


def contacts(a,y,phi):
    right = (R3+a,y)
    top3 = (a-R3/2,r+y)
    bottom3 = (a-R3/2,-r+y)
    plus = tangent_normal(right, 1)
    minus = tangent_normal(right, -1)
    Cplus = rotate((-c36,-g),phi)
    Cminus = rotate((-c36,g),phi)
    Rplus = -r-dot(Cplus,plus)
    Rminus = -r-dot(Cminus,minus)
    cross = det(plus,minus)
    t5 = ((Rplus*minus[1]-plus[1]*Rminus)/cross,
          (plus[0]*Rminus-Rplus*minus[0])/cross)
    wplus = add(Cplus,t5)
    wminus = add(Cminus,t5)
    top5 = add(rotate((c72,r),phi), t5)
    bottom5 = add(rotate((c72,-r),phi), t5)
    top_n = tangent_normal(top3, -1)
    error = dot(top5,top_n)-r
    return error, (right, top3, bottom3, top5, bottom5, wplus, wminus, plus, minus, top_n, t5)


def other_root(w,n):
    return sub(scale(-1/dot(w,w),w),n)


def argument_right(n):
    return (n[1]/n[0]).atan()


def area_on_curve(a,y,phi):
    _, data = contacts(a,y,phi)
    right, top3, bottom3, top5, bottom5, wplus, wminus, plus, minus, top_n, _ = data
    zero = arb_series([0]); origin=(zero,zero)
    ap,p_m = argument_right(plus), argument_right(minus)
    d_p,d_m = argument_right(other_root(wplus,plus)), argument_right(other_root(wminus,minus))
    top = pi/2-(top_n[0]/top_n[1]).atan()
    # Each entry is ending normal angle and the arc (0-radius means point).
    arc_data = [
        (ap, 0, right), (d_p,r,origin), (3*pi/10+phi,1,wplus),
        (top,0,top5), (5*pi/6,0,top3), (pi+p_m,1,right),
        (pi+ap,r,origin), (7*pi/6,1,right), (3*pi/2,0,bottom3),
        (pi+top,1,top3), (3*pi/2+phi,1,top5), (17*pi/10+phi,0,bottom5),
        (2*pi+d_m,1,wminus), (2*pi+p_m,r,origin), (2*pi,0,right),
    ]
    first = zero
    twice_area = zero
    for j,(last,radius,center) in enumerate(arc_data):
        _, next_radius, next_center = arc_data[(j+1)%len(arc_data)]
        n_first,n_last = normal(first),normal(last)
        twice_area += radius*det(center,sub(n_last,n_first)) + radius**2*(last-first)
        end = add(center,scale(radius,n_last))
        start_next = add(next_center,scale(next_radius,n_last))
        twice_area += det(end,start_next)
        first=last
    return twice_area/2, top


def main():
    a=interval('-.012154558715405409','-.012154558715405408')
    zero,eps=arb_series([0]),arb_series([0,1])
    E_y=contacts(a,eps,zero)[0][1]
    E_phi=contacts(a,zero,eps)[0][1]
    assert E_y > q('.9')
    y_linear=-E_phi/E_y
    err2=contacts(a,eps*y_linear,eps)[0][2]
    y_quadratic=-err2/E_y
    y=arb_series([0,y_linear,y_quadratic])
    H,top=area_on_curve(a,y,eps)
    assert q('.005') < top[1] < q('.006')
    assert q('.15') < H[2] < q('.16')
    assert q('.834780945775202') < H[0] < q('.834780945775205')
    assert 0 in H[1]
    source=Path(__file__)
    print(json.dumps({
        'python':platform.python_version(),'flint':flint.__version__,'precision_bits':ctx.prec,
        'script_sha256':sha256(source.read_bytes()).hexdigest(),
        'a_bracket':['-.012154558715405409','-.012154558715405408'],
        'top_equation_y_derivative':str(E_y),'top_equation_phi_derivative':str(E_phi),
        'y_linear':str(y_linear),'y_quadratic':str(y_quadratic),
        'top_angle_linear':str(top[1]),
        'area_coefficients':[str(H[j]) for j in range(3)],
        'scope':'Analytic contact-manifold Taylor coefficients from the positive twist side. First-order vanishing is a separate analytic identity.'
    },indent=2))


if __name__=='__main__':
    main()
