"""Transfer balanced three-normal padding witnesses to gains outside each old hull."""
import json
from fractions import Fraction as F
from itertools import combinations
from pathlib import Path

from certify_multi_host_points import arb, balance, dot, lower, points, q, two_hulls

HERE = Path(__file__).resolve().parent
R = q(F(707, 1000))
r = q(F(1, 2))
assert R**2 < 2*r**2
disk_point_area = arb.pi()*r**2 + r*(R**2-r**2).sqrt() - r**2*(r/R).acos()
assert disk_point_area > q(F('.838'))

raw = (HERE/'competing-fits-input.json').read_bytes()
certificate = json.loads((HERE/'multi-host-points.json').read_bytes())
cloud, _ = points(certificate)
result = {}
for name, host in two_hulls(raw).items():
    bounds = []
    conflicts = 0
    for sign, cells in certificate['proofs'][name].items():
        intervals = []
        for depth, idx, rows in cells:
            unit = [host.direction(tuple(row[:2]))[0] for row in rows]
            weights = balance(unit)
            assert sum(weights) == 1 and min(weights) > 0
            assert all(dot(u,u) == 1 for u in unit)
            acute = [(i,j) for i,j in combinations(range(3),2) if dot(unit[i],unit[j]) > 0]
            assert len(acute) <= 1
            cost = sum(w*w for w in weights)
            if acute:
                i,j = acute[0]
                cost += 2*weights[i]*weights[j]
                conflicts += 1
            d = lower(host, cloud, int(sign), depth, idx, rows)
            assert d > 0
            bounds.append((d*d/q(cost)).lower())
            intervals.append((F(idx,2**depth), F(idx+1,2**depth)))
        intervals.sort()
        assert intervals[0][0] == 0 and intervals[-1][1] == 1
        assert all(a[1]==b[0] for a,b in zip(intervals,intervals[1:]))
    minimum = min(bounds)
    # High-area branch uses the published actual host area upper bounds.
    host_upper = F('.836491') if name == 'Kstar' else F('.837657369')
    assert q(F('.838')-host_upper) > minimum
    result[name] = dict(cells=len(bounds), cells_with_acute_pair=conflicts,
                        area_increment_lower=str(minimum), old_area_upper=str(host_upper))
print(json.dumps(result, indent=2))
