"""Verify a displacement bound on every positive cap height for the R148 certificate."""
from fractions import Fraction as F
from hashlib import sha256
from math import isqrt
from pathlib import Path
import json
import sys

from certify_polygon_083731 import HERE, ORDERS, LAMBDA, arb, ar, ctx, diameter_squared, sampled_group


def checked(name, checksum):
    raw = (HERE/name).read_bytes()
    assert sha256(raw).hexdigest() == checksum
    return json.loads(raw)


def main(output_path):
    ctx.prec = 160
    replay = checked('tomar-six-083735-macbook-replay.json',
                     '100e0d535f6704925292eeaf3d6604514f5a18f66d80d4a4563abab3d6f7cdec')
    assert replay['complete'] and replay['original_leaves'] == 1046
    outer = checked('tomar-083735-stage-ledger.json', replay['outer_audit_sha256'])
    assert F(outer['target']) == F(replay['target']) == F('.83735')
    cloud = checked('finite-point-leaf-input.json', replay['geometry_sha256'])
    scale = LAMBDA
    rho = F('.000000555')
    epsilon = 1-scale
    base_radius = F(1,2)
    host_radius = F('.703')
    centered_support = F('.42264973')
    point_radius = F('.528619')
    point_loss = rho*point_radius
    ratio = point_loss/epsilon
    max_height = F('.204')
    original = {'denominator':cloud['denominator'],'points':cloud['points']}
    for x,y in original['points']:
        v = F(x,original['denominator']) + F('.01007')
        assert v*v + F(y,original['denominator'])**2 < point_radius**2
    for n in (3,5,7,9):
        assert 1-1/(2*(arb.pi()/(2*n)).cos()) > ar(centered_support)
    r, R = ar(scale*base_radius), ar(host_radius)
    spike = (arb.pi()-(r/R).acos())*r*r + r*(R*R-r*r).sqrt()
    goal = F(int(F('.83735')*2**57)+1,2**57)
    assert spike > ar(goal)
    assert host_radius + point_loss - base_radius < max_height
    assert host_radius**2 < 2*(scale*base_radius)**2
    v = (1+scale)*base_radius**2 + ratio**2*epsilon
    derivative = (1+scale)-(2*ratio*base_radius+v)/(4*base_radius**2)
    assert derivative > 0 and ratio < centered_support
    right = (1+scale)*max_height + max(F(0),(v-2*ratio*(base_radius+max_height))/(2*base_radius+max_height))
    margin = 4*(centered_support-ratio)*scale**2-right
    assert margin > F('.0067')
    # Normalize P30 rationally; contraction about its center differs only by translation.
    unit = 10**18
    d2 = diameter_squared(original)
    t = F(isqrt(unit**2*d2.denominator//d2.numerator),unit)
    expansion = t/(1-rho)
    groups = [sampled_group(n,arb.pi(),expansion*scale) for n in ORDERS]
    groups.append({'points':[[x*t.numerator,y*t.numerator] for x,y in original['points']],
                   'denominator':original['denominator']*t.denominator})
    diameters = [diameter_squared(g) for g in groups]
    assert all(v<=1 for v in diameters)
    floor = goal*expansion**2
    assert floor > F('.83735093')
    raw = (json.dumps(groups,separators=(',',':'))+'\n').encode()
    output_path.write_bytes(raw)
    return dict(threshold=str(goal), cap_height_bound=str(max_height), comparison_margin=str(margin),
                radius_area_lower=str(spike), rho=str(rho), sizes=[len(g['points']) for g in groups],
                squared_diameters=list(map(str,diameters)), area_floor=str(floor),
                point_sha256=sha256(raw).hexdigest())


if __name__ == '__main__':
    print(json.dumps(main(Path(sys.argv[1])),indent=2))
