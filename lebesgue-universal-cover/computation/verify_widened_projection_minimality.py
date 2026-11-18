#!/usr/bin/env python3
"""Check that each single-point deletion of the 18-point obstruction fits."""
import argparse
from fractions import Fraction
import hashlib
import json
from pathlib import Path

from verify_widened_projected_counterexample import exact_projection


def verify_deletions(host, data, duals, poses):
    obstruction = exact_projection(host,data,duals)
    points = data['points']
    pden = data['denominator']
    assert len(points) == len(set(map(tuple,points))) == 18
    assert sorted(p['deleted_index'] for p in poses) == list(range(len(points)))
    hden = host['planes_denominator']
    planes = [tuple(row[1:]) for row in host['rational_planes']]
    planes += [(a*hden,b*hden,hden) for a,b in ((1,0),(-1,0),(0,1),(0,-1))]
    results = []
    for pose in poses:
        deleted = pose['deleted_index']
        sign = pose['sign']
        assert type(deleted) is int and sign in (1,-1)
        t = Fraction(pose['half_tangent'])
        c,s = (1-t*t)/(1+t*t),2*t/(1+t*t)
        tx,ty = map(Fraction,pose['translation'])
        assert c*c+s*s == 1
        images = [(Fraction(c*x-s*sign*y,pden)+tx,Fraction(s*x+c*sign*y,pden)+ty)
                  for j,(x,y) in enumerate(points) if j != deleted]
        clearance = min((h-max(a*x+b*y for x,y in images))/(abs(a)+abs(b)) for a,b,h in planes)
        assert clearance > 0, deleted
        results.append({'deleted_index':deleted,'remaining_points':len(images),'reflection_sign':sign,
                        'strict_support_clearance_lower_exact':str(clearance),
                        'rational_rigid_pose':[str(v) for v in (c,s,tx,ty)]})
    return {'full_set_projection_obstruction':obstruction,'every_proper_subset_fits':True,
            'lowest_deletion_clearance_exact':str(min(Fraction(r['strict_support_clearance_lower_exact']) for r in results)),
            'single_deletions':results}


if __name__=='__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('host',type=Path)
    parser.add_argument('points',type=Path)
    parser.add_argument('duals',type=Path)
    parser.add_argument('deletion_poses',type=Path)
    parser.add_argument('output',type=Path)
    args = parser.parse_args()
    paths = (args.host,args.points,args.duals,args.deletion_poses)
    raw = [path.read_bytes() for path in paths]
    result = verify_deletions(*map(json.loads,raw))
    result['input_sha256'] = {path.name:hashlib.sha256(data).hexdigest() for path,data in zip(paths,raw)}
    args.output.write_text(json.dumps(result,indent=2)+'\n')
