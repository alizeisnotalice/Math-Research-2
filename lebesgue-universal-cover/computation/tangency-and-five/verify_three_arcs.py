"""Directed-Arb replay of three analytic angle intervals. No floating search."""
from fractions import Fraction as F
from hashlib import sha256
from pathlib import Path
import json
import sys
from flint import arb, ctx, __version__

ctx.prec = 192
def q(x):
    x = F(x)
    return arb(x.numerator) / x.denominator
def rad(x):
    return q(x)*arb.pi()/180
def u(x):
    return x.cos(), x.sin()
def n(x):
    return u(rad(x))
def dot(x, y):
    return x[0]*y[0] + x[1]*y[1]
def cross(x, y):
    return x[0]*y[1] - x[1]*y[0]
def scale(a, x):
    return tuple(a*z for z in x)
def add(x, y):
    return tuple(a+b for a,b in zip(x,y))
def bounds(x):
    def fraction(z):
        a,b = z.man_exp()
        return str(F(int(a))*F(2)**int(b))
    return [fraction(x.lower()),fraction(x.upper())]
def replay():
    half = q('1/2')
    root = arb(3).sqrt()
    r = 1/(2*rad(18).cos())
    b = r*(root*rad(6).cos()/2 + rad(72).cos())
    d = (r*rad(6).sin())**2
    s = d/(b+(b*b-d/4).sqrt())
    delta = r*rad(6).cos()+root*s/2-half
    assert 0<s<delta
    p = scale(s,n(330))
    ys = [add(scale(r,n(6+72*k)),p) for k in range(5)]
    assert 0<ys[0][1]<1/(2*root) and 0<r-s<1/root
    contacts = [scale(1/(half+delta),ys[k]) for k in (1,3)]
    for v,l,h in zip(contacts,(60,180),(120,240)):
        assert cross(n(l),v)>0 and cross(v,n(h))>0
    triples = [contacts+[n(300)], [n(a) for a in (84,228,300)], [n(a) for a in (90,234,306)]]
    intervals = [(6,9),(9,15),(15,18)]
    receipts = []
    for i,((lo,hi),v) in enumerate(zip(intervals,triples)):
        w = [cross(v[1],v[2]),cross(v[2],v[0]),cross(v[0],v[1])]
        assert all(a>0 for a in w)
        w = [a/sum(w) for a in w]
        def weighted(t):
            return sum(a*r*dot(z,u(t+rad(72*k))) for a,z,k in zip(w,v,(1,3,4)))
        sector = rad(lo).union(rad(hi))
        assert weighted(sector)>0
        hi_gap = weighted(rad(hi))-half-delta
        assert hi_gap>0
        if i:
            lo_gap = weighted(rad(lo))-half-delta
            assert lo_gap>0
        receipts.append(dict(angles=[lo,hi],lower_endpoint_surplus=bounds(lo_gap) if i else 'zero by contact',
                             upper_endpoint_surplus=bounds(hi_gap),concavity_margin=bounds(weighted(sector))))
    c = q('999/48000')
    T = 2*q('0.0001')/root
    budget = 2*c*rad('4.62').cos()*rad(2).sin() - T - q('0.001')
    assert budget>0
    def gibbs(delta_angle, zeta, subtraction):
        return half + 2*c*rad(F('4.62')+3*delta_angle).cos()*rad(30+2*delta_angle+zeta).cos()-subtraction
    return dict(exact_host_miss=bounds(delta),translation=bounds(s),pentagon_angles=receipts,
                gibbs_profile_exact=[bounds(gibbs(v,0,0)) for v in (0,10,15)],
                robust_H='1/10000',robust_slants='1/1000',robust_zeta_deg=2,
                robust_budget_margin=bounds(budget),
                gibbs_profile_robust=[bounds(gibbs(v,2,T)) for v in (0,10,15)],
                flint=__version__,python=sys.version,precision=ctx.prec,source_sha256=sha256(Path(__file__).read_bytes()).hexdigest())
if __name__=='__main__':
    j=replay()
    Path(sys.argv[1]).write_text(json.dumps(j,indent=2)+'\n')
    print({k:[float(F(z[0])) for z in j[k]] for k in ('gibbs_profile_exact','gibbs_profile_robust')})
    print([(z['angles'],float(F(z['upper_endpoint_surplus'][0]))) for z in j['pentagon_angles']])
