"""Certify the 18,049-point lower bound using the published .83744 geometric tree."""
import argparse
from fractions import Fraction as F
from functools import reduce
from hashlib import sha256
import json
from math import isqrt
from pathlib import Path
import platform
import sys

import flint
import numpy as np


def serialized(value):
    return (json.dumps(value, separators=(',', ':')) + '\n').encode()


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('computation', type=Path)
    p.add_argument('output', type=Path)
    args = p.parse_args()
    assert __debug__, 'Run with Python assertions enabled.'
    comp = args.computation.resolve()
    sys.path.insert(0, str(comp))
    from certify_polygon_083731 import ar, arb, ctx, lower_radius
    from certify_shallow_cap_transfer import diameter_squared
    from tangent_polygons_finite_cover import generate

    def read(name):
        return json.loads((comp/name).read_bytes())

    def digest(name):
        return sha256((comp/name).read_bytes()).hexdigest()

    ctx.prec = 160
    goal = F(120687823102244742, 2**57)
    source = read('native-wall-083744-audit.json')
    assert source['complete'] and F(source['goal']) == goal
    native = read('finite-point-lift-083744.json')
    dense_name = native['dense']['filename']
    assert digest('finite-point-lift-083744.json') == source['ordinary_input_sha256']
    assert digest(dense_name) == source['dense_sha256']
    rows = native['rows']
    dense = read(dense_name)
    pure = set(map(tuple, native['five_only']))
    ordinary = set(tuple(row[:2]) for row in rows)
    mesh = set(map(tuple, dense['original_leaves']))
    assert (len(ordinary), len(pure), len(mesh)) == (2833,124,6)
    assert len(ordinary | pure | mesh) == 2963
    assert len({tuple(row[:2]) for row in rows if row[5] > 0} | mesh) == source['native_point_leaves'] == 1636
    assert digest('finite-point-leaf-input.json') == native['point_sha256']
    assert digest('finite-point-global-input.json') == native['bank_sha256']
    cloud = read('finite-point-leaf-input.json')
    bank = read('finite-point-global-input.json')
    directions = cloud['directions']
    assert len(directions) == 105 and len(bank['cells']) == 384
    active = sorted({j for cell in bank['cells'] for item in cell for j in item[:3]})
    assert len(active) == 54 and all(0 <= j < len(directions) for j in active)

    r, s, R = F(1,2), F('.999999006'), F('.703')
    eps, wall, norm = 1-s, F('.4742688878'), F('.02816')
    assert F(source['scale']) == s and F(source['max_b3_norm_squared']) < norm**2
    assert all(lower_radius(n,arb.pi()) > ar(wall) for n in (5,7,9)) and r > wall
    circle = 1/(2*(arb.pi()/6).cos())
    corners = []
    for j in range(3):
        sy,cx = (2*arb.pi()*j/3).sin_cos()
        corners.append((circle*cx,circle*sy))
    with np.load(comp/'five-body-radial-material.npz', allow_pickle=False) as z:
        assert list(z['counts'][:2]) == [96,96]
        pts = z['points'][96:192]
    inner = [(ar(F(int(x),2**28)), ar(F(int(y),2**28))) for x,y in pts]
    assert all((x-vx)**2 + (y-vy)**2 < 1 for x,y in inner for vx,vy in corners)
    direct, clipped, direct_lowers, clipped_uppers = [], [], [], []
    for j,(x,y) in enumerate(directions):
        length = arb(x*x+y*y).sqrt()
        nx,ny = arb(x)/length,arb(y)/length
        lower = reduce(arb.max, (nx*px+ny*py for px,py in inner))
        upper = 1+reduce(arb.min, (nx*vx+ny*vy for vx,vy in corners))
        if lower >= ar(wall):
            direct.append(j)
            direct_lowers.append(lower)
        else:
            assert ar(norm)+ar(s)*upper < ar(r-eps*wall)
            clipped.append(j)
            clipped_uppers.append(upper)
    assert len(direct)==64 and len(clipped)==41

    h, rho, q, D = F('.0121'), F('.000000867709'), F('.528619'), F('.495')
    m = rho*q/eps
    den = cloud['denominator']
    assert all((F(x,den)+F('.01007'))**2+F(y,den)**2 < q*q for x,y in cloud['points'])
    sr = s*r
    assert (arb.pi()-ar(sr/R).acos())*ar(sr)**2+ar(sr)*(ar(R**2-sr**2)).sqrt() > ar(goal)
    assert D*D > R*R-sr*sr and R*R < 2*sr*sr
    five = read('five-global-minimum-single-worker-replay.json')
    assert five['complete'] and arb(five['minimum_interval']) > ar(F('.837295'))
    height_margin = s*s*F('.837295') + sr*h*h/D - goal
    v = (1+s)*r*r+eps*m*m
    slope = 1+s-(2*m*r+v)/(4*r*r)
    margin = 4*(wall-m)*s*s-(1+s)*h-max(0,(v-2*m*(r+h))/(2*r+h))
    assert min(height_margin,slope,margin)>0 and m<wall
    dsq = diameter_squared(cloud)
    t = F(isqrt(10**36*dsq.denominator//dsq.numerator),10**18)
    f = t/(1-rho)
    groups,budget = generate((s*f)**2,F(1))
    assert s*f*F(budget['common_dilation'])<1
    groups[-1] = dict(body='P30',denominator=t.denominator*den,
                      points=[[t.numerator*x,t.numerator*y] for x,y in cloud['points']])
    final_diam = diameter_squared(groups[-1])
    assert final_diam == t*t*dsq <= 1
    checks = []
    for i,group in enumerate(groups):
        count = len(group['points'])
        denom = group['denominator']
        numerator = (budget['squared_diameter_numerators'][i] if i<5 else int(final_diam*denom**2))
        assert F(numerator,denom**2) <= 1
        checks.append(dict(body=group['body'],count=count,coordinate_denominator=denom,
                           pairs=count*(count-1)//2,maximum_squared_distance_numerator=numerator,
                           sha256=sha256(serialized(group)).hexdigest()))
    assert sum(z['count'] for z in checks)==18049
    coord = serialized(groups)
    lower = f*f*goal
    assert lower > F('.83744145457600263406')
    result = dict(area_lower=str(lower),display_lower='.83744145457600263406',underlying_goal=str(goal),
                  scale=str(s),height=str(h),point_shrink=str(rho),wall=str(wall),factor=str(f),
                  height_margin=str(height_margin),cap_margin=str(margin),cap_slope_lower=str(slope),
                  direct_normals=direct,disk_dominated_normals=clipped,active_normals=active,
                  direct_min_support=str(reduce(arb.min,direct_lowers)),
                  dominated_max_support=str(reduce(arb.max,clipped_uppers)),sampling=budget['sampling'],
                  common_dilation=str(budget['common_dilation']),groups=checks,points_sha256=sha256(coord).hexdigest(),
                  input_sha256={name:digest(name) for name in ('native-wall-083744-audit.json','finite-point-lift-083744.json',dense_name,'finite-point-leaf-input.json','finite-point-global-input.json','five-body-radial-material.npz','five-global-minimum-single-worker-replay.json')},
                  python=platform.python_version(),flint=flint.__version__,numpy=np.__version__,precision_bits=ctx.prec)
    assert result['points_sha256']=='8aa93e43d4a0b8d10aed3a59aea51b663d10252c8a9030923e07a859aecc9bc9'
    args.output.mkdir(parents=True,exist_ok=False)
    (args.output/'normal-gap-083744-points.json').write_bytes(coord)
    (args.output/'normal-gap-083744-audit.json').write_text(json.dumps(result,indent=2)+'\n')
    print(result['area_lower'])

if __name__=='__main__':
    main()
