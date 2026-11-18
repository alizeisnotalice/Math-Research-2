"""Certify all translations and all angles for a rational point set against ten outer half-planes."""
import argparse
from fractions import Fraction
from hashlib import sha256
from itertools import combinations
import json
from math import gcd
from pathlib import Path


def balanced_triples(normals):
    """Integer nonnegative vertices of the normalized Farkas cone, duplicates removed."""
    def det(a, b):
        return a[0]*b[1] - a[1]*b[0]

    candidates = set()
    for i, j, k in combinations(range(len(normals)), 3):
        a, b, c = normals[i], normals[j], normals[k]
        weights = [det(b, c), det(c, a), det(a, b)]
        if max(weights) <= 0:
            weights = [-v for v in weights]
        if min(weights) < 0 or sum(weights) == 0:
            continue
        divisor = gcd(*weights)
        row = tuple((r, w//divisor) for r, w in zip((i, j, k), weights) if w)
        assert all(sum(w*normals[r][axis] for r, w in row) == 0 for axis in (0, 1))
        candidates.add(row)
    return sorted(candidates)


def one_quadrant(points, walls, triples, point_den, wall_den, quadrant, parity, target):
    # For t in [0,1], the rotation is pi*quadrant/2 + 2*atan(t).
    vectors = []
    for x, y in points:
        y *= parity
        for _ in range(quadrant):
            x, y = -y, x
        vectors.append((x, y))
    ab = [[(a*x+b*y, b*x-a*y) for x, y in vectors] for _, a, b, _ in walls]
    offsets = [point_den*c for _, _, _, c in walls]

    def endgap(t, a, b, c, norm):
        v, d = t.numerator, t.denominator
        return Fraction((d*d-v*v)*a + 2*d*v*b - (d*d+v*v)*c, (d*d+v*v)*norm)

    ledger = []
    def visit(left, right, depth):
        mid = (left+right)/2
        v, d = mid.numerator, mid.denominator
        cosnum, sinnum, den = d*d-v*v, 2*d*v, d*d+v*v
        def score(pair):
            return cosnum*pair[0] + sinnum*pair[1]
        pairs = [max(choices, key=score) for choices in ab]
        excess = [score(pair) - den*c for pair, c in zip(pairs, offsets)]
        best = max(triples, key=lambda row: Fraction(sum(w*excess[r] for r,w in row),sum(w for _,w in row)))
        a = sum(w*pairs[r][0] for r,w in best)
        b = sum(w*pairs[r][1] for r,w in best)
        c = sum(w*offsets[r] for r,w in best)
        scale = sum(w for _, w in best)*point_den*wall_den
        lower = min(endgap(t, a,b,c,scale) for t in (left,right))
        if lower > target:
            # Both endpoint projections are positive. Every intermediate unit direction
            # is a positive combination of the endpoints with coefficient sum >= 1.
            assert c >= 0 and target >= 0
            ledger.append((depth, lower))
            return
        if depth >= 30 or endgap(mid,a,b,c,scale) <= target:
            raise ValueError(f"unclosed angle interval {quadrant},{parity},{left},{right}")
        visit(left, mid, depth+1)
        visit(mid, right, depth+1)
    visit(Fraction(0),Fraction(1),0)
    return len(ledger), max(d for d,_ in ledger), min(g for _,g in ledger)


def check_outer(comp, walls, wall_den):
    """Every point of the actual host satisfies the frozen rational half-planes."""
    import sys
    sys.path.insert(0, str(comp.resolve()))
    from certify_ah_tapered_minkowski_closure import geometry
    from flint import arb, ctx
    g, _ = geometry(comp)
    ctx.prec = 192
    actual = dict((name,(u,b)) for name,u,b in g['active'])
    actual.update((f'hex_{i}',(u,arb(1)/2)) for i,u in enumerate(g['hexn']))
    values = []
    for name,a,b,c in walls:
        n,h = actual[name]
        err = abs(arb(a)/wall_den-n[0]) + abs(arb(b)/wall_den-n[1])
        # The host is inside the normalized hexagon, hence both coordinates have magnitude <1.
        room = arb(c)/wall_den - (h+err)
        assert room > 0
        values.append(room.lower())
    return str(min(values))


def verify(path, comp=None):
    raw = path.read_bytes()
    doc = json.loads(raw)
    p, d = doc['points'], doc['coordinate_denominator']
    walls, n = doc['walls'], doc['normal_denominator']
    diam2 = max((x-X)**2+(y-Y)**2 for (x,y),(X,Y) in combinations(p,2))
    assert diam2 < d*d
    assert all(c > 0 and a*a+b*b <= (n*1000001//1000000)**2 for _,a,b,c in walls)
    triples = balanced_triples([(a,b) for _,a,b,_ in walls])
    goal = Fraction(doc.get('target','0.000002'))
    grid = [one_quadrant(p,walls,triples,d,n,q,s,goal) for s in (-1,1) for q in range(4)]
    least = min(m for _,_,m in grid)
    out = dict(input_sha256=sha256(raw).hexdigest(),points=len(p),diameter_squared=str(Fraction(diam2,d*d)),
               planar_rows=len(walls),balanced_rows=len(triples),target=str(goal),
               closed_angle_leaves=sum(k for k,_,_ in grid),max_dyadic_depth=max(h for _,h,_ in grid),
               leaves_by_hand_quadrant=[k for k,_,_ in grid],uniform_signed_overflow_lower=str(least),
               euclidean_offset_excluded=str(least*Fraction(1000000,1000001)))
    if comp is not None:
        out['minimum_outward_rounding_slack'] = check_outer(comp, walls, n)
    return out


if __name__ == '__main__':
    if not __debug__:
        raise SystemExit('Assertions must be enabled.')
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('input',type=Path)
    ap.add_argument('output',type=Path)
    ap.add_argument('--host',type=Path)
    args = ap.parse_args()
    data = verify(args.input,args.host)
    args.output.write_text(json.dumps(data,indent=2)+'\n')
    print(data['points'],data['closed_angle_leaves'],float(Fraction(data['uniform_signed_overflow_lower'])))
