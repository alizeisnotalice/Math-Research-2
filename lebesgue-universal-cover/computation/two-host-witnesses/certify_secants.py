"""Propagate balanced certificates from both endpoints of each orientation cell."""
import argparse
import json
import sys
from fractions import Fraction
from hashlib import sha256
from pathlib import Path

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))
from cross_common import A, C, Host, arb, balance, ctx, dot, hf, hp, pi, q, rot


def fraction_floor(value, scale=10**12):
    integer = (value.lower() * scale).floor().unique_fmpz()
    return Fraction(int(integer), scale)


def preimage(direction, angle, sign):
    x, y = rot(-angle, direction)
    return x, sign * y


def smooth_cell(rows, weights, endpoints, sign, support, cosine):
    host = sum(w * row[2] for w, row in zip(weights, rows))
    values = [
        sum(w * support(preimage(row[1], angle, sign)) for w, row in zip(weights, rows))
        for angle in endpoints
    ]
    lower = min(value.lower() for value in values)
    assert lower < 1
    return (1 - (1 - lower) / cosine - host).lower()


def material_cell(rows, weights, endpoints, sign, points):
    # Select one material point per normal; validity holds for any selection.
    center = sum(endpoints) / 2
    host = sum(w * row[2] for w, row in zip(weights, rows))
    selected = []
    for row in rows:
        direction = preimage(row[1], center, sign)
        selected.append(max(points, key=lambda point: dot(point, direction).lower()))
    corners = [
        sum(w * dot(p, preimage(row[1], angle, sign))
            for w, p, row in zip(weights, selected, rows)) for angle in endpoints
    ]
    least = min(c.lower() for c in corners)
    assert least > 0
    return (least - host).lower()


def certify(output):
    raw = json.loads((ROOT / 'certificate-labels.json').read_bytes())
    finite_source = ROOT / 'sparse' / 'points.json'
    points = [tuple(q(Fraction(z)) for z in row) for row in json.loads(finite_source.read_bytes())['points']]
    five_source = ROOT / 'smallest' / 'points5.json'
    five_rows = json.loads(five_source.read_bytes())
    selected = five_rows['indices_from_original_0based']
    parent = json.loads(finite_source.read_bytes())['points']
    assert [parent[i] for i in selected] == five_rows['points'] and len(set(selected)) == 5
    five_points = [tuple(q(Fraction(z)) for z in row) for row in five_rows['points']]
    delta_source = ROOT / 'smallest' / 'five-labels.json'
    changes = json.loads(delta_source.read_bytes())
    diameter = {}
    for name, coords in [('sparse',parent),('five-points',five_rows['points'])]:
        xy = [tuple(map(Fraction,row)) for row in coords]
        largest = max((x-u)**2 + (y-v)**2 for x,y in xy for u,v in xy)
        assert largest < 1
        diameter[name] = str(largest)
    K = raw['K']
    assert K == 2048 and len(points) == 11
    assert len(raw['seeds']) == 512
    assert ((arb.pi()/2) > 2*pi/K) and all(len(v) == K for v in raw['rows'].values())
    get = lambda name: json.loads((C/name).read_bytes())
    hosts = {
        'r155':Host(get('free-p30-reoptimized-pose.json'),False),
        'r175':Host(get('seven-point-families-pose.json')['pose'],True)
    }
    # Moving the whole grid by this angle avoids the atan2 cut at material knots.
    offset = A('0.000000000001')
    rows_by_host = {name:[h.row(tuple(seed)) for seed in raw['seeds']] for name,h in hosts.items()}
    bounds = {}
    for name, lookup in rows_by_host.items():
        for sign in (1,-1):
            for label, support, subs in [('r151',hp,1),('r216',hf,1),('sparse',None,4),('five-points',None,4)]:
                source = 'r151' if support is None else label
                key = f'{name}-{source}-{sign}'
                labels = raw['rows'][key]
                overrides = changes[key] if label == 'five-points' else {}
                material = five_points if label == 'five-points' else points
                assert all(str(int(k)) == k and 0 <= int(k) < K for k in overrides)
                cosine = (pi / (K*subs)).cos()
                least = None
                for i, indices in enumerate(labels):
                    indices = overrides.get(str(i),indices)
                    assert all(type(j) is int and 0 <= j < len(lookup) for j in indices)
                    assert len(indices) == len(set(indices)) == 3
                    rows = [lookup[j] for j in indices]
                    weights = tuple(q(w) for w in balance([row[0] for row in rows]))
                    for k in range(subs):
                        start = i * subs + k
                        endpoints = [2*pi*(start+e)/(K*subs)+offset for e in (0,1)]
                        if support is None:
                            gain = material_cell(rows, weights, endpoints, sign, material)
                        else:
                            gain = smooth_cell(rows, weights, endpoints, sign, support, cosine)
                        least = gain if least is None or gain < least else least
                rounded = fraction_floor(least)
                assert rounded > 0
                bounds[f'{name}-{label}-{sign}'] = {
                    'lower':str(rounded), 'cells':len(labels)*subs, 'method':
                        'fixed material sine interpolation' if support is None else 'curvature upper bound one'
                }
                print(name,label,sign,rounded,flush=True)
    sources = [ROOT/'cross_common.py',ROOT/'source-manifest.json',ROOT/'certificate-labels.json',finite_source,five_source,delta_source,Path(__file__)]
    data = {'diameter_squared':diameter,'precision_bits':ctx.prec,'grid_shift_radians':str(Fraction(1,10**12)),'results':bounds,
            'sources':{str(p.relative_to(ROOT)):sha256(p.read_bytes()).hexdigest() for p in sources}}
    output.write_text(json.dumps(data,indent=2)+'\n')


if __name__=='__main__':
    arguments = argparse.ArgumentParser()
    arguments.add_argument('computation',type=Path)
    arguments.add_argument('output',type=Path)
    args=arguments.parse_args()
    if not __debug__:
        raise RuntimeError('The certificate requires enabled assertions')
    certify(args.output)
