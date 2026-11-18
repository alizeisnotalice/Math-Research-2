#!/usr/bin/env python3
"""Replay positive similarity projections in the widened rational host."""
import argparse
from fractions import Fraction
import hashlib
import json
from math import lcm
from pathlib import Path

import audit_optimized_eh_chord as g
from certify_eh_general_relocation import exact_fraction
from certify_eh_no_interior_roots import verify
from certify_eh_polyline_obstruction import witness_points


def rational_pose(spec):
    z = Fraction(spec['half_tangent'])
    c = (1-z*z)/(1+z*z)
    s = 2*z/(1+z*z)
    tx,ty = map(Fraction,spec['translation'])
    assert c*c+s*s == 1 and spec['sign'] in (-1,1)
    return c,s,tx,ty,spec['sign']


def rational_point_margin(planes, points, pden, pose):
    c,s,tx,ty,sign = pose
    d = lcm(c.denominator,s.denominator,tx.denominator,ty.denominator)
    C,S,X,Y = map(int,(c*d,s*d,tx*d,ty*d))
    finite = [(C*x-S*sign*y+X*pden, S*x+C*sign*y+Y*pden) for x,y in points]
    D = d*pden
    minimum = min(Fraction(b*D-max(a*x+v*y for x,y in finite),D*(abs(a)+abs(v))) for a,v,b in planes)
    assert minimum > 0
    return minimum


def curvature(terms):
    assert terms and all(k >= 3 and k%2 for k,_ in terms)
    series = [(k,(1-k*k)*a) for k,a in terms]
    crude = Fraction(1,2)-sum(abs(a) for _,a in series)
    if crude > 0:
        return
    count = 1024
    tolerance = sum(k*k*abs(a) for k,a in series)*Fraction(44,7*count)**2/8
    for j in range(count):
        angle = 2*g.PI*j/count
        assert g.HALF+sum(exact_fraction(a)*(k*angle).cos() for k,a in series) > exact_fraction(tolerance)


def curved_margin(planes, spec, pose, source):
    c,s,tx,ty,sign = map(exact_fraction,pose)
    if spec['kind'] == 'fourier':
        terms = [(k,Fraction(a)) for k,a in (source['terms'] if 'terms' in source else spec['terms'])]
        curvature(terms)
    else:
        verify(source)
        points,_ = witness_points(source)
        image = [(c*x-s*sign*y+tx, s*x+c*sign*y+ty) for x,y in points.values()]
    gaps = []
    for a,b,h in planes:
        a,b,h = map(g.arb,(a,b,h))
        length = g.norm((a,b))
        if spec['kind'] == 'fourier':
            cosine = (a*c+b*s)/length
            polynomials = [g.ONE,cosine]
            for _ in range(2, max(k for k,_ in terms)+1):
                polynomials.append(2*cosine*polynomials[-1]-polynomials[-2])
            support = length*(g.HALF+sum(exact_fraction(coef)*polynomials[k] for k,coef in terms))+a*tx+b*ty
        else:
            support = max(g.dot((a,b),p).upper() for p in image)
        gap = (h-support)/length
        assert gap > 0,(spec['name'],str(gap))
        gaps.append(gap.lower())
    value = min(gaps)
    mantissa,exponent = map(int,value.man_exp())
    return Fraction(mantissa)*Fraction(2)**exponent


def verify_fits(host, tests, archive):
    hden = host['planes_denominator']
    planes = [tuple(row[1:]) for row in host['rational_planes']]
    planes += [(a*hden,b*hden,hden) for a,b in ((1,0),(-1,0),(0,1),(0,-1))]
    assert all(0 < h*h <= a*a+b*b for a,b,h in planes)
    results = []
    for test in tests:
        pose = rational_pose(test)
        source = json.loads((archive/test['source']).read_bytes()) if 'source' in test else test
        if 'record_index' in test:
            source = source[test['record_index']]
        if test['kind'] == 'rational_points':
            pden = source.get('denominator',source.get('point_denominator'))
            minimum = rational_point_margin(planes,source['points'],pden,pose)
            num_points = len(source['points'])
        else:
            assert test['kind'] in ('fourier','witness_points')
            minimum = curved_margin(planes,test,pose,source)
            num_points = 18 if test['kind']=='witness_points' else None
        # Scaling both the rotation coefficients and translation by this
        # factor yields a feasible point outside the unit circle of projection.
        rho = 1+Fraction((minimum*10**15/2)//1,10**15)
        assert rho > 1 and rho-1 < minimum
        results.append({'name':test['name'],'kind':test['kind'],'number_of_points':num_points,
                        'strict_clearance_lower_exact':str(minimum),'feasible_similarity_radius_exact':str(rho),
                        'rational_rigid_pose':[str(x) for x in pose]})
    return results


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('host',type=Path)
    parser.add_argument('test_poses',type=Path)
    parser.add_argument('archive',type=Path)
    parser.add_argument('output',type=Path)
    args=parser.parse_args()
    raw_host,raw_tests=(path.read_bytes() for path in (args.host,args.test_poses))
    tests = json.loads(raw_tests)
    fits=verify_fits(json.loads(raw_host),tests,args.archive)
    inputs = sorted({row['source'] for row in tests if 'source' in row})
    args.output.write_text(json.dumps({'host_sha256':hashlib.sha256(raw_host).hexdigest(),
                                      'test_poses_sha256':hashlib.sha256(raw_tests).hexdigest(),
                                      'source_sha256':{name:hashlib.sha256((args.archive/name).read_bytes()).hexdigest()
                                                       for name in inputs},
                                      'fits':fits},indent=2)+'\n')
