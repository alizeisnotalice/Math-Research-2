#!/usr/bin/env python3
"""Pure-integer witness for the n=100 endpoint in the analytic threshold proof."""
from fractions import Fraction
import json

n = 100
N = n * n
m = N - 4
a = 2 * n + 6
# Equivalent to 2*N^(-a) > (4/5)^m.
lhs = 2 * pow(5, m)
rhs = pow(4, m) * pow(N, a)
assert lhs > rhs
# The lower analytic bound after log estimates is positive at n=100.
poly = Fraction(n * n, 45) - Fraction(3 * n, 5) - Fraction(8, 9)
assert poly > 0
print(json.dumps({"status": "passed", "n": n,
                  "integer_inequality": "2*5^(N-4) > 4^(N-4)*N^(2n+6)",
                  "positive_rational_log_lower_without_log2": str(poly)}))
