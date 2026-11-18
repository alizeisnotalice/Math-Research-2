"""Exclude every rigid placement of a width-one ten-piece curvature witness in the R203 host."""
from argparse import ArgumentParser
from fractions import Fraction
from hashlib import sha256
import json
from pathlib import Path
import sys

from flint import arb, arb_mat, ctx

GRID = 300
BREAKS = (0, 6, 7, 64, 65, 125, 126, 188, 248, 249, 300)
ROWS = ('K-C3', 'deep', 'I-E3', 'I-E3', 'I-E3', 'I-E3', 'I-E3', 'I-E3',
        'taper_30', 'K-C3', 'K-C3', 'K-C3')
FIXED = Fraction(-223359536107, 10**12)
COMPACT_SHA = 'ca131c3c44c07754f9abf7be22938e4d82d8c9942fc7b4c616f57fcd18b68733'


def run(comp):
    sys.path.insert(0, str(comp.resolve()))
    from certify_ah_tapered_minkowski_closure import geometry
    assert sha256((comp / 'certify_ah_tapered_release.py').read_bytes()).hexdigest() == COMPACT_SHA
    g, _ = geometry(comp)
    ctx.prec = 192
    pi = arb.pi()
    knots = [j*pi/GRID for j in BREAKS]
    half = arb(1)/2
    fixed = arb(FIXED.numerator)/FIXED.denominator

    def pieces(x):
        return (half,x[0],-half,fixed,half,x[1],-half,half,x[2],-half)

    def jump(y, v, f, dx):
        sn, cs = dx.sin_cos()
        return f + (y-f)*cs + v*sn, -(y-f)*sn + v*cs

    def states(x):
        z = [(arb(0), x[3])]
        for j, f in enumerate(pieces(x)):
            z.append(jump(*z[-1],f,knots[j+1]-knots[j]))
        return z

    def half_value(index, angle, z, fs):
        return jump(*z[index],fs[index],angle-knots[index])

    def at_integer(k,den,z,fs):
        k = k % (2*GRID*den)
        sign = 1 if k < GRID*den else -1
        k %= GRID*den
        ix = next(j for j in range(len(fs)) if k < den*BREAKS[j+1])
        y,v = half_value(ix,pi*k/(GRID*den),z,fs)
        return sign*y,sign*v

    def equations(x):
        z = states(x)
        fs = pieces(x)
        return (at_integer(100,1,z,fs)[0],at_integer(200,1,z,fs)[0],z[-1][0],z[-1][1]+x[3])

    # These four defining equations impose the other two hexagon contacts,
    # antiperiodicity of y at pi, and antiperiodicity of its derivative.
    c = equations([arb(0)]*4)
    columns = [equations([arb(int(i==j)) for i in range(4)]) for j in range(4)]
    mat = arb_mat([[columns[j][i]-c[i] for j in range(4)] for i in range(4)])
    det = mat.det()
    assert abs(det) > arb(1)/10**9
    solved = mat.solve(arb_mat([[-v] for v in c]))
    x = [solved[i,0] for i in range(4)]
    assert all(abs(v)<half for v in x[:3]) and abs(fixed)<half and abs(x[3])<arb(1)/5
    assert all(v.contains(0) for v in equations(x))
    fs,z = pieces(x),states(x)

    def input_value(direction):
        theta = arb.atan2(direction[1],direction[0])
        if theta < 0:
            theta += 2*pi
        sign = 1
        if theta > pi:
            theta -= pi
            sign = -1
        assert 0<theta<pi
        ix = next(j for j in range(len(fs)) if knots[j]<theta<knots[j+1])
        return half + sign*half_value(ix,theta,z,fs)[0]

    def compat(k,der=False):
        i = int(der)
        return sum((-1)**j*at_integer(k+1000*j,10,z,fs)[i] for j in range(3))

    mesh = 1000
    step = pi/(3*mesh)
    # On [0,pi], |y| <= |s| + (1/2) integral_0^pi sin(u) du < 1.2.
    # Consequently |F''| <= 3*(1/2+1.2) < 6 almost everywhere.
    eps = arb(6)*step**2/8
    top = max(compat(k).upper() for k in range(1,mesh)).upper()
    angle_drop = -(top+eps)
    slope_left = -(compat(0,True)+6*step)
    slope_right = compat(mesh,True)-6*step
    assert slope_left>arb(1)/100 and slope_right>arb(1)/100
    assert angle_drop>arb(1)/50000

    actual = {name:(n,h) for name,n,h in g['active']}
    poses = [(k,s) for k in range(0,360,60) for s in (-1,1)]
    margins = []
    for (k,s),label in zip(poses,ROWS,strict=True):
        normal,h = actual[label]
        pulled = g['preimage'](normal,k,s)
        violation = input_value(pulled)-h
        assert violation>arb(32)/10**6
        margins.append(dict(rotation_degrees=k,reflection_sign=s,wall=label,violation=str(violation)))
    return dict(complete=True,body='h=1/2+y, y(0)=0, y antiperiodic with derivative',
                breakpoints_in_pi_over_300=BREAKS,fixed_f_64=str(FIXED),defining_matrix_det=str(det),
                parameters=list(map(str,x)),grid_intervals=mesh,interior_negative_clearance=str(angle_drop),
                left_slope_clearance=str(slope_left),right_slope_clearance=str(slope_right),
                min_pose_violation=str(min(arb(row['violation']).lower() for row in margins)),poses=margins,
                host_source_sha256=COMPACT_SHA,
                host_receipt_sha256=sha256((comp/'ah-tapered-release.json').read_bytes()).hexdigest())


if __name__ == '__main__':
    parser=ArgumentParser(description=__doc__)
    parser.add_argument('computation',type=Path)
    parser.add_argument('output',type=Path)
    args=parser.parse_args()
    if not __debug__:raise SystemExit('Assertions must be enabled.')
    result=run(args.computation)
    result['verifier_sha256']=sha256(Path(__file__).read_bytes()).hexdigest()
    args.output.write_text(json.dumps(result,indent=2)+'\n')
    print(result['min_pose_violation'],result['interior_negative_clearance'])
