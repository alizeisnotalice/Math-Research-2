#!/usr/bin/env python3
"""Exact four-dimensional similarity constraints and two-dimensional exclusion."""
import argparse
from fractions import Fraction
import hashlib
import json
from pathlib import Path


BOUND = 4
DIRECTION_DEN = 10**9


def clip(vertices, normal, limit):
    a,b = normal
    result = []
    for p,q in zip(vertices,vertices[1:]+vertices[:1]):
        x = a*p[0]+b*p[1]-limit
        y = a*q[0]+b*q[1]-limit
        if x <= 0:
            result.append(p)
        if x*y < 0:
            lam = x/(x-y)
            result.append(tuple(u+(v-u)*lam for u,v in zip(p,q)))
    return result


def squared_diameter(points):
    return max((ax-bx)**2+(ay-by)**2 for i,(ax,ay) in enumerate(points) for bx,by in points[:i])


def exact_projection(host, input_data, cuts):
    points = input_data['points']
    pden = input_data['denominator']
    hden = host['planes_denominator']
    assert type(pden) is int and pden > 0 and type(hden) is int and hden > 0
    assert all(type(x) is int and type(y) is int for x,y in points)
    assert all(25*(x*x+y*y) < 9*pden*pden for x,y in points)
    largest = squared_diameter(points)
    assert largest == input_data['diameter_squared_numerator']
    assert 81*pden*pden < 100*largest < 100*pden*pden
    planes = [tuple(row[1:]) for row in host['rational_planes']]
    assert all(type(a) is int and type(b) is int and type(c) is int and c>0 and (a or b) for a,b,c in planes)
    assert len(cuts) == 2
    branches = []
    for sign, rows in zip((1,-1),cuts):
        poly = [(Fraction(x),Fraction(y)) for x,y in ((-BOUND,-BOUND),(BOUND,-BOUND),(BOUND,BOUND),(-BOUND,BOUND))]
        for normal, weights in rows:
            assert len(normal) == 2 and all(type(v) is int for v in normal) and any(normal)
            u,v = [Fraction(k,DIRECTION_DEN) for k in normal]
            remainder = [u,v,Fraction(0),Fraction(0)]
            limit = Fraction(0)
            for (i,j), value in weights:
                assert type(i) is int and type(j) is int and 0 <= i < len(planes) and 0 <= j < len(points)
                w = Fraction(value)
                assert w >= 0
                a,b,c = planes[i]
                x,y = points[j]
                y *= sign
                coefficients = (a*x+b*y,-a*y+b*x,a*pden,b*pden)
                remainder = [d-w*Fraction(e,hden*pden) for d,e in zip(remainder,coefficients)]
                limit += w*Fraction(c,hden)
            # Uniform variable bound makes the exact residual of any nonnegative
            # rational dual harmless. No exact LP optimum is presumed.
            limit += BOUND*sum(map(abs,remainder))
            poly = clip(poly,(u,v),limit)
            assert poly
        maximum = max(x*x+y*y for x,y in poly)
        assert maximum < 1
        numerator = -((-maximum*10**15)//1)
        branches.append({'reflection_sign':sign,'projected_outer_cut_count':len(rows),
                         'projected_outer_vertex_count':len(poly),
                         'max_squared_radius_upper':f'{numerator//10**15}.{numerator%10**15:015d}',
                         'max_squared_radius_exact':str(maximum)})
    return {'input_count':len(points),'strict_diameter_below_one':True,
            'squared_diameter_exact':str(Fraction(largest,pden*pden)),
            'no_euclidean_placement':True,'branches':branches}


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('host',type=Path)
    parser.add_argument('points',type=Path)
    parser.add_argument('duals',type=Path)
    parser.add_argument('output',type=Path)
    opts = parser.parse_args()
    paths = (opts.host,opts.points,opts.duals)
    data = [path.read_bytes() for path in paths]
    result = exact_projection(*map(json.loads,data))
    result['source_sha256'] = {path.name:hashlib.sha256(raw).hexdigest() for path,raw in zip(paths,data)}
    opts.output.write_text(json.dumps(result,indent=2)+'\n')
