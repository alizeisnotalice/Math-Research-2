"""Rational outer tangent vertices with one common similarity loss for five inputs."""
import argparse
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations
import json
from pathlib import Path

from flint import arb, ctx
from certify_nine import S, PI_UPPER, append_point, ar

HERE = Path(__file__).resolve().parent


def diameter_excess(body, cells):
    if body:
        x = (PI_UPPER/(2*body*cells))**2
        return x/(2-x)
    x = (PI_UPPER/cells)**2
    return 3*x/(8-4*x)


def radius_budget(body, target, lower):
    def enough(m):
        e = diameter_excess(body, m)
        denominator = S+4 + (S*e).__ceil__()
        return lower*(S-5)**2 > target*denominator**2
    lo, hi = 2, 8
    while not enough(hi):
        hi *= 2
    while hi-lo > 1:
        mid = (lo+hi)//2
        if enough(mid):
            hi = mid
        else:
            lo = mid
    return hi if body or hi%2 else hi+1


def tangent_points(n, m):
    pi = arb.pi()
    points = []
    if not n:
        assert m % 2
        radius = 1/(2*(pi/m).cos())
        for j in range(m):
            sn, cs = (2*pi*j/m).sin_cos()
            append_point(points, radius*cs, radius*sn)
        return sorted(set(points))
    radius = 1/(2*(pi/(2*n)).cos())
    assert 1-radius > ar(F(2,5))
    tip = 1/(pi/(2*n*m)).cos()
    for k in range(n):
        angle = 2*pi*k/n
        sn, cs = angle.sin_cos()
        x, y = radius*cs, radius*sn
        append_point(points,x,y)
        for j in range(m):
            normal = angle + pi*(ar(1-F(1,2*n)+F(2*j+1,2*n*m)))
            sn, cs = normal.sin_cos()
            append_point(points,x+tip*cs,y+tip*sn)
    return sorted(set(points))


def generate(target, lower):
    sizes = [radius_budget(n,target,lower) for n in (0,3,5,7,9)]
    excess = max(diameter_excess(n,m) for n,m in zip((0,3,5,7,9),sizes,strict=True))
    denominator = S+4+(S*excess).__ceil__()
    bound = lower*F(S-5,denominator)**2
    assert bound > target and arb.pi() < ar(PI_UPPER)
    groups, diameters = [], []
    for n,m in zip((0,3,5,7,9),sizes,strict=True):
        pts = tangent_points(n,m)
        assert len(pts)==(n*(m+1) if n else m)
        diam = max((x-a)**2+(y-b)**2 for (x,y),(a,b) in combinations(pts,2))
        assert diam <= denominator**2
        groups.append(dict(body=n,denominator=denominator,points=pts))
        diameters.append(diam)
    fixed = json.loads((HERE/'finite-point-leaf-input.json').read_bytes())
    groups.append(dict(body='P30',denominator=fixed['denominator'],points=fixed['points']))
    counts=[len(g['points']) for g in groups]
    return groups,dict(sampling=sizes,counts=counts,total=sum(counts),denominator=denominator,
                       squared_diameter_numerators=diameters,curved_lower=str(lower),
                       common_dilation=str(F(denominator,S-5)),area_lower=str(bound),target=str(target))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('coordinates',type=Path)
    parser.add_argument('--target',type=F,default=F('.837295'))
    parser.add_argument('--lower',type=F,default=F('.83730'))
    args=parser.parse_args()
    assert 0 < args.target < args.lower*F(S-5,S+5)**2
    ctx.prec=160
    proofs = {
        F('.83730'): ('finite-point-lift-083730-replay.json',
                      'd4a93ccaa6e8458c2427f26a533f3447f3961cef2e4b04721fdc396907034dd0', 541),
        F('.83731'): ('finite-point-lift-083731-replay.json',
                      'a395ae3c05008ce57a0adbd27c1dc47c71346b1eb119347f2a1d66050015bc96', 614),
        F('.83732'): ('finite-point-lift-083732-replay.json',
                      'cd21296af9d2d2e963588cd531aaee321364570ee182e0bf06a979b11d7f4846', 705),
    }
    name, expected, leaves = proofs[args.lower]
    proof = (HERE/name).read_bytes()
    summary = json.loads(proof)
    assert sha256(proof).hexdigest() == expected and summary['complete']
    if 'input_sha256' in summary:
        lift = (HERE/'finite-point-lift-083732.json').read_bytes()
        assert summary['input_sha256'] == sha256(lift).hexdigest()
        assert summary['covered'] == leaves and F(summary['threshold']) == args.lower
        geometry_hash = json.loads(lift)['point_input_sha256']
    else:
        assert summary['original_leaves'] == leaves and F(summary['target']) == args.lower
        geometry_hash = summary['geometry_sha256']
    assert geometry_hash == sha256((HERE/'finite-point-leaf-input.json').read_bytes()).hexdigest()
    groups,receipt=generate(args.target,args.lower)
    receipt.update(curved_proof=name, curved_proof_sha256=expected)
    raw=(json.dumps(groups,separators=(',',':'))+'\n').encode()
    args.coordinates.write_bytes(raw)
    receipt['sha256']=sha256(raw).hexdigest()
    print(json.dumps(receipt,indent=2))
