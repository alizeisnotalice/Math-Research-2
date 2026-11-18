"""Certify all-orientation noncontainment of a smooth unit-width body."""

import argparse
import json
import math
import platform
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations
from pathlib import Path
from typing import NamedTuple

class Model(NamedTuple):
    source: str
    a3: F
    a5: F
    a15: F
    cells: int
    padding: F
    spare: F

MODELS = {
    "five": Model(
        "002-five-body-certificate.json", F(1,80), -F(19,1200), F(0),
        512, F(1,500), F(1,20000)
    ),
    "eight": Model(
        "013-eight-body-certificate.json", F(0), F(1,200), F(3,1792),
        1024, F(17,20000), F(1,40000)
    ),
}
GRID = 10**6
HOST_RADIUS_UPPER = F(581, 1000)
CHORD_X = F(2,5)
CHORD_HALF_LENGTH = F(3,10)
assert CHORD_X**2 + CHORD_HALF_LENGTH**2 == F(1,4)


def unit(w, sign):
    assert -1 <= w <= 1 and sign in (-1, 1)
    denom = 1 + w*w
    return sign*(1-w*w)/denom, sign*2*w/denom


def dot(x, y):
    return x[0]*y[0] + x[1]*y[1]


def cross(x,y):
    return x[0]*y[1]-x[1]*y[0]


def cosine3(x):
    return 4*x**3-3*x


def cosine5(x):
    return 16*x**5-20*x**3+5*x


def exact_support(x,model):
    # x is the cosine of the difference between the world normal and body angle.
    return F(1,2) + model.a3*cosine3(x) + model.a5*cosine5(x) + model.a15*cosine3(cosine5(x))


def balance(normals):
    u,v,w = normals
    weights = (cross(v,w),cross(w,u),cross(u,v))
    if all(x < 0 for x in weights):
        weights = tuple(-x for x in weights)
    assert min(weights) > 0
    total = sum(weights)
    weights = tuple(q/total for q in weights)
    assert all(sum(weight*normal[k] for weight, normal in zip(weights,normals)) == 0
               for k in range(2))
    return weights


def parse_host(raw):
    host = json.loads(raw)
    polygon = [tuple(F(p,host["coordinate_denominator"]) for p in row)
               for row in host["hull_grid_points"]]
    return polygon, F(host["epsilon"])


def host_support(polygon, eps, direction):
    return max(dot(direction,v) for v in polygon) + eps


def cells(n):
    # w in [-1,1] covers the right half. Multiplying by -1 covers the left half.
    # Angular radius <= 2/n, since angle derivative is 2/(1+w**2).
    for sign in (1,-1):
        for j in range(n):
            yield unit(F(2*j+1-n,n),sign)


def numerical_rows(normals, host_values, model):
    # Suggestions only. The following exact verifier recomputes all quantities.
    import numpy as np

    qs = np.asarray(normals,dtype=float)
    b = np.asarray(host_values,dtype=float)
    triples = np.array(list(combinations(range(len(qs)),3)))
    i,j,k = triples.T
    det = qs[:,None,0]*qs[None,:,1]-qs[:,None,1]*qs[None,:,0]
    w0,w1,w2 = det[j,k],det[k,i],det[i,j]
    positive = (w0>0)&(w1>0)&(w2>0)
    negative = (w0<0)&(w1<0)&(w2<0)
    nondegenerate = positive | negative
    triples = triples[nondegenerate]
    w = np.stack((w0[nondegenerate],w1[nondegenerate],w2[nondegenerate]))
    w /= w.sum(axis=0)
    i,j,k = triples.T
    rows = []
    for u in cells(model.cells):
        x = qs @ np.asarray(u,dtype=float)
        term = (0.5 + float(model.a3)*cosine3(x)+float(model.a5)*cosine5(x)
                + float(model.a15)*cosine3(cosine5(x)) - b)
        scores = w[0]*term[i]+w[1]*term[j]+w[2]*term[k]
        rows.append(triples[int(scores.argmax())].tolist())
    return rows


def direction_seeds(polygon):
    seeds = []
    for a,b in zip(polygon,polygon[1:]+polygon[:1],strict=True):
        x,y = float(b[1]-a[1]),float(a[0]-b[0])
        norm = math.hypot(x,y)
        x,y = x/norm,y/norm
        sign = 1 if x >= 0 else -1
        seed = round(GRID*sign*y/(1+sign*x))
        seeds.append([sign,seed])
    return seeds


def exact_verify(raw, certificate, model):
    assert certificate["host_sha256"] == sha256(raw).hexdigest()
    assert certificate["target_extra_minkowski_padding"] == str(model.padding)
    n = certificate["cells_per_half_circle"]
    assert n == model.cells
    lipschitz = 3*abs(model.a3) + 5*abs(model.a5) + 15*abs(model.a15)
    curvature = F(1,2) - (8*abs(model.a3) + 24*abs(model.a5) + 224*abs(model.a15))
    assert curvature > 0
    rows = certificate["balanced_normal_indices"]
    assert len(rows) == 2*n
    polygon, eps = parse_host(raw)
    assert all(dot(q,q)<(HOST_RADIUS_UPPER-eps)**2 for q in polygon)
    assert 0 < eps < HOST_RADIUS_UPPER
    seeds = certificate["normal_seeds"]
    assert all(type(sign) is int and type(k) is int for sign,k in seeds)
    normals = [unit(F(k,GRID),sign) for sign,k in seeds]
    assert all(dot(u,u)==1 for u in normals)
    b = [host_support(polygon,eps,u) for u in normals]
    worst = None
    for j,(rot,indices) in enumerate(zip(cells(n),rows,strict=True)):
        assert len(set(indices)) == 3 and all(type(i) is int and 0<=i<len(normals) for i in indices)
        qs = [normals[i] for i in indices]
        ws = balance(qs)
        exact_centre = sum(w*(exact_support(dot(q,rot),model)-b[i])
                           for q,i,w in zip(qs,indices,ws,strict=True))
        lower = exact_centre - lipschitz*2/n
        assert lower > model.padding + model.spare, (j,float(lower))
        worst = (lower,j) if worst is None else min(worst,(lower,j))
    return worst[0], worst[1], curvature


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("cover", choices=MODELS)
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--generate",type=Path)
    group.add_argument("--verify",type=Path)
    args=parser.parse_args()
    model=MODELS[args.cover]
    raw=Path(__file__).with_name(model.source).read_bytes()
    if args.generate:
        polygon,eps=parse_host(raw)
        seeds=direction_seeds(polygon)
        normals=[unit(F(k,GRID),sign) for sign,k in seeds]
        host_values=[host_support(polygon,eps,u) for u in normals]
        certificate = {
            "host_sha256":sha256(raw).hexdigest(),
            "target_extra_minkowski_padding":str(model.padding),
            "cells_per_half_circle":model.cells,
            "normal_seeds":seeds,
            "balanced_normal_indices":numerical_rows(normals,host_values,model),
        }
    else:
        certificate = json.loads(args.verify.read_bytes())
    lower, index, curvature = exact_verify(raw,certificate,model)
    summary = {
        "python":platform.python_version(),
        "worst_cell":index,
        "worst_certified_margin":str(lower),
        "worst_certified_margin_display":float(lower),
        "target_extra_minkowski_padding":str(model.padding),
        "certified_cells":2*model.cells,
        "curvature_radius_lower_bound":str(curvature),
        "coordinate_halfwidth_for_robust_exclusion":str(model.padding*F(2,5)),
        "extra_minkowski_padding_still_excluded_on_this_neighborhood":str(model.spare),
        "host_radial_upper_bound":str(HOST_RADIUS_UPPER),
        "area_increase_lower_bound_for_adding_any_copy_to_fixed_true_hull":str(
            CHORD_HALF_LENGTH*model.padding**2/(HOST_RADIUS_UPPER+model.padding-CHORD_X)
        ),
    }
    if args.cover == "eight":
        coefficient_radius = F(1,50000)
        assert curvature > 224*coefficient_radius
        assert lower - coefficient_radius > model.padding
        summary["coefficient_l1_radius_for_closed_family"] = str(coefficient_radius)
        summary["extra_padding_still_excluded_over_coefficient_family"] = str(model.padding)
    if args.generate:
        args.generate.write_text(json.dumps(certificate,indent=2)+"\n")
    print(json.dumps(summary,indent=2))


if __name__ == "__main__":
    main()
