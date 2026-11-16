"""Rational sanity checks; the all-height arguments are proved separately."""
from fractions import Fraction as F
from pathlib import Path
import json

cases = 0
for degree in range(1, 17):
    for n in range(4 * degree + 2, 4 * degree + 102):
        m = (n - 1) // (2 * degree)
        lower = sum((F(n - 1 - degree*q, n - q - 1)
                     for q in range(1, m+1)), F())
        assert lower >= F(m, 4) >= F(n-1, 8*degree)-F(1, 4) >= F(n, 16*degree)
        cases += 1
small = 1 + 2 + F(4, 3)/(1-F(2, 7))
pruned = 1 + 3 + 3 + 3*F(2, 5)/(1-F(2, 5))
assert small < 5 and pruned == 9
product = F(1)
for j in range(1, 65):
    product *= 1-F(1, 2**j)
    assert product > F(1, 4)
out = dict(scope='Finite rational checks, not a substitute for the all-height proof',
           promotion_cases=cases, small_binary_majorant=str(small),
           pruned_majorant=str(pruned), finite_G_product_64_gt_quarter=True,
           B_upper_constant=720, degree3_A2_lower_denominator=48,
           ratio_lower_denominator=720*48, original_grid_problem='unresolved')
Path(__file__).with_name('constants-audit.json').write_text(json.dumps(out, indent=2)+'\n')
print(json.dumps(out))
