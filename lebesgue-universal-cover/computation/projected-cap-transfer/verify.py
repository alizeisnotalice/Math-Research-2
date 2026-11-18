"""Certify the R221 same-placement cap transfer with a projection-dependent height."""
import argparse
from fractions import Fraction as F
from hashlib import sha256
import json
from math import isqrt
from pathlib import Path

from flint import arb, ctx


def sqrt_upper(q):
    scale = 10**24
    num = q.numerator * scale**2
    lo = isqrt(num // q.denominator)
    return F(lo + int(lo**2 * q.denominator < num), scale)


def certify(comp, rho):
    ctx.prec = 192
    paths = [comp/'six-input-083745/points/normal-gap-083745-audit.json',
             comp/'five-global-minimum-single-worker-replay.json',
             comp/'finite-point-leaf-input.json']
    base, five, cloud = [json.loads(path.read_bytes()) for path in paths]
    s, a, G = [F(base[k]) for k in ('scale', 'wall', 'underlying_goal')]
    assert (s, a, G) == (F('.999999006'), F('.4742688878'), F(120689264254125501, 2**57))
    r, R = F(1, 2), F('.703')
    radius, center_x = F('.5286170'), F('-.01007206')
    points = [(F(x, cloud['denominator']), F(y, cloud['denominator'])) for x,y in cloud['points']]
    assert len(points) == 30
    assert max((x-center_x)**2+y*y for x,y in points) < radius**2
    diameter_squared = max((x-u)**2+(y-v)**2 for x,y in points for u,v in points)
    M = F('.8372950062210252434563')
    assert five['complete'] and arb(five['minimum_interval']) > arb(M.numerator)/M.denominator
    sr, eps = s*r, 1-s
    m = radius*rho/eps
    assert 0 < m < a < r
    deficit = G-s*s*M
    v = (1+s)*r*r+eps*m*m
    assert deficit > 0 and 1+s-(2*m*r+v)/(4*r*r) > 0
    budget = 4*(a-m)*s*s
    checked = []
    def cover(lo, hi, depth):
        assert depth <= 60
        # e <= x-r. Both branches of this elementary majorant increase in x.
        bound = (1+s)*(hi-r) + max(0, (v-2*m*hi)/(hi+r))
        if lo > r:
            # cap area at the same projected support x controls its own height e.
            height = min(hi-r, sqrt_upper(deficit*sqrt_upper(hi*hi-sr*sr)/sr))
            correction = max(0, (v-2*m*lo)/(lo*lo-r*r))
            bound = min(bound,height*(1+s+correction))
        if bound < budget:
            checked.append((lo,hi,budget-bound,depth))
        else:
            mid = (lo+hi)/2
            cover(lo,mid,depth+1)
            cover(mid,hi,depth+1)
    cover(r,R+m*eps,0)
    assert checked[0][0] == r and checked[-1][1] == R+m*eps
    assert all(left[1] == right[0] for left,right in zip(checked,checked[1:]))
    factor2 = 1/(diameter_squared*(1-rho)**2)
    assert s*s*factor2 < 1
    lower = factor2*G
    old = F(base['area_lower'])
    assert lower > old
    places = 24
    n = lower.numerator*10**places//lower.denominator
    return dict(complete=True,scale=str(s),wall=str(a),threshold=str(G),
        point_shrink=str(rho),center_x=str(center_x),radius_upper=str(radius),
        projection_rate=str(m),five_minimum_lower=str(M),point_diameter_squared=str(diameter_squared),
        coverage=[str(r),str(R+m*eps)],projection_slabs=len(checked),
        maximum_depth=max(z[3] for z in checked),minimum_comparison_slack=str(min(z[2] for z in checked)),
        factor_squared=str(factor2),area_lower=str(lower),gain_over_221=str(lower-old),
        strict_decimal_lower=f'{n//10**places}.{n%10**places:0{places}d}',
        sources={str(path.relative_to(comp)):sha256(path.read_bytes()).hexdigest() for path in paths})


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('computation',type=Path)
    parser.add_argument('output',type=Path)
    parser.add_argument('--rho',type=F,default=F('.000000872664'))
    args = parser.parse_args()
    assert __debug__
    result = certify(args.computation,args.rho)
    args.output.write_text(json.dumps(result,indent=2)+'\n')
    print(result['strict_decimal_lower'],result['projection_slabs'],result['maximum_depth'])
