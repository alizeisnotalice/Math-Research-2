"""Enclose the analytic positive-part curvature dual over two full parameter bands."""
from argparse import ArgumentParser
from fractions import Fraction as F
from hashlib import sha256
import json
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from certify_b1_classification import Geometry, PINS, rat
from dual_analytic import calculate
from dual_derivative import derivative_bound


SEGMENTS = {
    49: ('.3828125', '.3867'), 50: ('.3866', '.3945'),
    51: ('.3944', '.4023'), 52: ('.4022', '.4092'),
    410: ('.4091', '.4109'), 411: ('.4108', '.4119'), 53: ('.4118', '.4145'),
    57: ('.4442', '.4478'), 449: ('.4477', '.4501'), 58: ('.4500', '.4570'),
    59: ('.4569', '.4648'), 60: ('.4647', '.4727'),
    61: ('.4725', '.4805'), 62: ('.4803', '.484375'),
}


def certify(geometry, weights, key):
    lo,hi = map(F,SEGMENTS[key])
    stack = [(lo+(hi-lo)*k/256,lo+(hi-lo)*(k+1)/256,0) for k in range(256)]
    nodes = leaves = depth = 0
    upper = None
    while stack:
        a,b,d = stack.pop()
        nodes += 1
        depth = max(depth,d)
        mid = (a+b)/2
        try:
            center = calculate(geometry,weights,mid,mid)
            slope = derivative_bound(geometry,weights,a,b)
            bound = center.upper()+abs(slope).upper()*rat((b-a)/2)
        except (AssertionError, ZeroDivisionError):
            bound = None
        if bound is not None and bound < 0:
            leaves += 1
            value = bound.upper()
            upper = value if upper is None else max(upper,value)
            continue
        assert d < 22, (key,a,b)
        stack.extend(((mid,b,d+1),(a,mid,d+1)))
    return {'j':key,'lo':str(lo),'hi':str(hi),'nodes':nodes,'leaves':leaves,
            'depth':depth,'strict_upper':str(upper)}


if __name__ == '__main__':
    parser = ArgumentParser(description=__doc__)
    parser.add_argument('computation', type=Path)
    parser.add_argument('output', type=Path)
    args = parser.parse_args()
    geometry = Geometry(args.computation.resolve())
    files = ('dual-parameters.json','additional-dual-parameters.json')
    weights = {}
    for name in files:
        weights.update(json.loads((HERE/name).read_bytes()))
    assert set(map(int,weights)) == set(SEGMENTS)
    rows = [certify(geometry,weights[str(k)],k) for k in SEGMENTS]
    union = []
    for row in sorted(rows,key=lambda x:F(x['lo'])):
        a,b = F(row['lo']),F(row['hi'])
        if union and a <= union[-1][1]:
            union[-1][1] = max(b,union[-1][1])
        else:
            union.append([a,b])
    assert union == [[F(49,128),F('.4145')],[F('.4442'),F(62,128)]]
    sources = (*files,'dual_analytic.py','dual_derivative.py','certify_b1_classification.py',Path(__file__).name)
    result = {'segments':rows,'bands':[[str(a),str(b)] for a,b in union],
              'geometry_sha256':PINS,'source_sha256':{n:sha256((HERE/n).read_bytes()).hexdigest() for n in sources}}
    args.output.write_text(json.dumps(result,indent=2)+'\n')
    print(result['bands'],sum(row['leaves'] for row in rows),flush=True)
