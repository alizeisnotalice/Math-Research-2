#!/usr/bin/env python3
"""Cross-check family evaluations against independently rebuilt primary inputs."""
import copy
import json
import sys
from fractions import Fraction as F
from pathlib import Path

import periodic_halo
from vertical_scaling import VerticalFamily

data = json.loads(Path(sys.argv[1]).read_text())
family = VerticalFamily(data)
points = (F(999,1000), F(1999,2000), F(1), F(2001,2000), F(1001,1000))
whole_upper = family.interval_upper(points[0], points[-1])
results = {}
for beta in points:
    scaled = copy.deepcopy(data)
    scaled['Py'] = str(F(data['Py'])*beta)
    scaled['pts'] = [[x, str(F(y)*beta), w] for x,y,w in data['pts']]
    primary = periodic_halo.calculate(scaled)
    inner, outer = family.evaluate(beta)
    # The primary routine has an additional safe pruning test, so enclosure
    # containment is sufficient even if the exact boundary rectangles differ.
    pi, po = F(primary['inner_area']), F(primary['outer_area'])
    assert inner == pi and outer == po
    assert outer <= whole_upper
    results[str(beta)] = {'inner_area':str(inner), 'outer_area':str(outer), 'primary_equal':True}
for left, right in zip(points, points[1:]):
    upper = family.interval_upper(left, right)
    assert family.evaluate(left)[1] <= upper and family.evaluate(right)[1] <= upper
    assert upper <= whole_upper
    results[str(left)+':'+str(right)] = {'interval_upper':str(upper)}
print(json.dumps({'points':results, 'whole_upper':str(whole_upper)},indent=2))
