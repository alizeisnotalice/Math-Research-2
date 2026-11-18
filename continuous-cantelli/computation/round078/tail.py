"""Exact Taylor-coefficient check for the two formal odd-flux series."""
import json
import math

def odd_product(n):
    return math.prod(range(1, n + 1, 2))

orders = list(range(1, 81))
for n in orders:
    a = odd_product(2*n - 1)
    s = odd_product(2*n + 1)
    assert -2*s + 4*(n - 1)*a == -6*a
print(json.dumps({"orders_checked": len(orders), "maximum_order": max(orders),
                  "first_two_normalized_left_coefficients": [-6, -18]}))
