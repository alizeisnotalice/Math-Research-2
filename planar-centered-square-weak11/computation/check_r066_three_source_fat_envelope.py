"""Exact arithmetic control for the three-source fat envelope obstruction.

This script checks source geometry and rational analytic lower certificates.
It does not compute a halo union; the production and independent geometry
programs have separate outputs and scope.
"""
from fractions import Fraction as F
from pathlib import Path
import hashlib
import json
import os
import platform
import time


def main():
    assert os.path.ismount('/Volumes/data')
    started = time.perf_counter()
    alpha, epsilon = F(1023, 1024), F(1, 1024)
    points = [(F(0), F(0)), (alpha, alpha), (alpha+F(1,8), alpha)]
    masses = [F(1,5), F(4,5), epsilon]
    rows = []
    for mask in range(1, 8):
        ids = [i for i in range(3) if mask >> i & 1]
        lo = [min(points[i][j] for i in ids) for j in range(2)]
        hi = [max(points[i][j] for i in ids) for j in range(2)]
        closure = sum(1 << i for i in range(3)
                      if all(lo[j] <= points[i][j] <= hi[j] for j in range(2)))
        mass = sum((masses[i] for i in ids), F(0))
        span2 = max(hi[j]-lo[j] for j in range(2))**2
        rows.append({'mask':mask, 'mass':str(mass), 'span_squared':str(span2),
                     'closure_mask':closure, 'eligible':span2 < mass})
    eligible = [r for r in rows if r['eligible']]
    assert [r['mask'] for r in eligible] == [1, 2, 3, 4, 6]
    maximum = max(F(r['mass']) for r in eligible)
    assert maximum == 1
    assert [r['mask'] for r in eligible if F(r['mass']) == maximum] == [3]
    high = [r['mask'] for r in eligible if r['mask'] & 3
            and r['closure_mask'] == r['mask'] and F(r['mass']) > F(3,4)]
    assert high == [2,3,6]
    reached = 1
    while True:
        newer = reached
        for row in eligible:
            if row['mask'] & reached:
                newer |= row['mask']
        if newer == reached:
            break
        reached = newer
    assert reached == 7
    domain_margin = (alpha+F(1,8))**2-F(35,33)*(1+epsilon)
    assert domain_margin == F(6982433,34603008) > 0
    assert F(3,4) < alpha < 1
    assert F(433,250)**2 < 3
    assert F(97,250)**2 > F(3,20)
    budget_margin = alpha**2+F(93,80)-epsilon+alpha*F(433,250)-F(97,250)-F(7,2)
    strong_margin = F(33,80)-epsilon-F(97,250)-F(7,4)+F(433,250)
    assert budget_margin == F(246397,131072000) > 0
    assert strong_margin == F(707,128000) > 0
    assert F(9,64)**2 < F(3,16)  # residual halo lies inside C
    assert F(9,64)**2 < F(1,5)   # residual halo lies inside heavy Q
    deletion_upper = 1+(1-alpha)**2
    assert deletion_upper == 1+F(1,1048576) < F(7,2)
    root = Path(__file__).resolve().parents[1]
    input_path = root/'computation/inputs/r066-three-source-fat-envelope.json'
    input_obj = {'name':'r066-three-source-fat-envelope',
                 'points':[[str(t) for t in p] for p in points],
                 'masses':[str(m) for m in masses],
                 'initial_masses':[str(m) for m in masses],
                 'S_indices':[0,1], 'S_mask':3, 'theta':'3/4',
                 'scope':'Unique fat maximum block, connected and nonterminal; fixed theta envelope fails, true deletion remains small.'}
    input_path.write_text(json.dumps(input_obj, indent=2)+'\n')
    output = {'status':'PASS_EXACT_SEVEN_SUBSET_AND_RATIONAL_CERTIFICATES',
              'python':platform.python_version(),
              'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
              'input_sha256':hashlib.sha256(input_path.read_bytes()).hexdigest(),
              'source_subsets':rows, 'high_closed_masks':high,
              'unique_maximum_mask':3, 'connected':True, 'both_axes_fat':True,
              'nonterminal_margin':str(domain_margin),
              'strict_budget_excess_rational_lower':str(budget_margin),
              'strict_strong_criterion_excess_rational_lower':str(strong_margin),
              'true_deletion_analytic_upper':str(deletion_upper),
              'seconds':time.perf_counter()-started,
              'scope':'Exact source and rational analytic checks only; no halo area union is computed by this program.'}
    output_path = root/'computation/results/round066/three-source-fat-arithmetic.json'
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(json.dumps(output, indent=2)+'\n')
    print(json.dumps({k:v for k,v in output.items() if k!='source_subsets'}, indent=2))


if __name__ == '__main__':
    main()
