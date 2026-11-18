#!/usr/bin/env python3
"""Verify fixed rational poses of the R153 obstruction in low-area hulls."""
import argparse
from fractions import Fraction
from hashlib import sha256
import json
from math import isqrt, lcm
from pathlib import Path

from flint import arb, ctx

from certify_five_body_contact_stationary import BOX, interval, model
from certify_nine import S, append_point, ar, arc_body, hull


def inner_hull(pose, segments):
    pi = arb.pi()
    samples = []
    for j in range(8*segments):
        angle = pi*j/(4*segments)
        append_point(samples, angle.cos()/2, angle.sin()/2)
    arc_body(samples,segments,pi,3,arb(0),*pose[:2])
    for n,j in ((5,2),(7,5),(9,8)):
        arc_body(samples,segments,pi,n,*pose[j:j+3])
    # Every rounded grid point is within 2/S of its true body. The convex
    # hull contains the centered disk (1/2)B, so a common rescaling by
    # S/(S+4) puts this entire rational grid hull inside the true hull.
    return hull(samples), S+4


def positivity(poly, host_denominator, pose, cloud):
    sign = pose['sign']
    assert sign in (1,-1)
    q = Fraction(pose['half_tangent'])
    c,s = (1-q*q)/(1+q*q),2*q/(1+q*q)
    assert c*c+s*s == 1
    tx,ty = map(Fraction,pose['translation'])
    denominator = cloud['denominator']
    scale = lcm(c.denominator,s.denominator,tx.denominator,ty.denominator)
    C,V,X,Y = map(int,(c*scale,s*scale,tx*scale,ty*scale))
    image = [(C*x-V*sign*y+X*denominator,V*x+C*sign*y+Y*denominator)
             for x,y in cloud['points']]
    point_denom = denominator*scale
    gaps = []
    for (x,y),(u,v) in zip(poly,poly[1:]+poly[:1]):
        a,b = v-y,x-u
        assert a or b
        upper = max(a*p+b*q for p,q in image)
        norm_upper = isqrt(a*a+b*b)+1
        bound = Fraction((a*x+b*y)*point_denom-upper*host_denominator,
                         host_denominator*point_denom*norm_upper)
        assert bound > 0
        gaps.append(bound)
    return min(gaps)


def test(host, poses, cloud, segments):
    poly, d = inner_hull(host,segments)
    assert len(poly) >= 3
    return {'rational_inner_hull_vertices':len(poly),
            'strict_support_clearance_exact':{str(pose['sign']):str(positivity(poly,d,pose,cloud)) for pose in poses}}


def verify(archive, cloud, poses, segments):
    ctx.prec = 192
    n = cloud['denominator']
    points = cloud['points']
    assert len(points)==18 and n>0 and all(type(x) is int and type(y) is int for x,y in points)
    maximum = max((x-u)**2+(y-v)**2 for i,(x,y) in enumerate(points) for u,v in points[:i])
    assert maximum == cloud['diameter_squared_numerator'] and maximum < n*n
    assert all(sorted(item['sign'] for item in rows)==[-1,1] for rows in poses.values())
    files = {name:(archive/name).read_bytes() for name in ('free-p30-reoptimized-pose.json',
            'five-body-contact-basin.json','five-body-contact-stationary.json','competing-fits-input.json')}
    contact = json.loads(files['five-body-contact-basin.json'])
    stationary = json.loads(files['five-body-contact-stationary.json'])
    assert stationary['stationary_angle_box'] == [list(row) for row in BOX]
    _, fit, _, _ = model(*(interval(*row) for row in BOX))
    hosts = {'R155_escape':(list(map(ar,json.loads(files['free-p30-reoptimized-pose.json'])[:11])),'R155'),
             'five_body_contact_stationary':(fit,'contact'),
             'five_body_contact_box_center':(list(map(ar,contact['raw_center'])),'contact'),
             'R108_escape':(list(map(ar,json.loads(files['competing-fits-input.json'])['zero_gain_five']['params'])),'R108')}
    results = {name:test(body,poses[key],cloud,segments) for name,(body,key) in hosts.items()}

    # Transfer the same reflected placement from the box center to every
    # member of the larger raw rectangular contact basin.
    half = list(map(Fraction,contact['raw_halfwidth']))
    sqrt2_upper = Fraction('1.415')
    assert sqrt2_upper**2 > 2
    error = sqrt2_upper*max(half[:2])
    for size,start,lip in ((5,2,Fraction('.163')),(7,5,Fraction('.115')),(9,8,Fraction('.089'))):
        assert (arb.pi()/(2*size)).tan()/2 < ar(lip)
        error = max(error,lip*half[start]+sqrt2_upper*max(half[start+1:start+3]))
    center_margin = Fraction(results['five_body_contact_box_center']['strict_support_clearance_exact']['-1'])
    assert center_margin > error
    return {'segments_per_arc':segments,'strict_squared_diameter_exact':str(Fraction(maximum,n*n)),
            'checked_hosts':results,'uniform_contact_box':{'reflection_sign':-1,
                'Hausdorff_upper_exact':str(error),'remaining_clearance_exact':str(center_margin-error)},
            'archived_source_sha256':{name:sha256(raw).hexdigest() for name,raw in files.items()}}


if __name__=='__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('archive',type=Path)
    parser.add_argument('cloud',type=Path)
    parser.add_argument('placements',type=Path)
    parser.add_argument('output',type=Path)
    parser.add_argument('--segments',type=int,default=128)
    opts = parser.parse_args()
    cloud, poses = (path.read_bytes() for path in (opts.cloud,opts.placements))
    result = verify(opts.archive,json.loads(cloud),json.loads(poses),opts.segments)
    result['input_sha256'] = {path.name:sha256(raw).hexdigest() for path,raw in zip((opts.cloud,opts.placements),(cloud,poses))}
    opts.output.write_text(json.dumps(result,indent=2)+'\n')
