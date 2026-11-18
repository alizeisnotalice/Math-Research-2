"""Deterministic checks of the resolvent and Gaussian thin-valley exponents."""
import json
import platform
from fractions import Fraction as F
from math import exp, log, sqrt

def logsum(xs):
    m = max(xs)
    return m + log(sum(exp(v - m) for v in xs))

def surrogate_rate(z, high, low, alpha):
    # Gaussian background exactly integrated; thin plateau density sampled at
    # its center, times exp(-alpha*n^2). This is a surrogate, not an exact graph.
    terms = [(1 - high) * z*z/2]
    terms += [-low*z*z/2 + z*n - (.5+alpha)*n*n for n in range(1, 2*z+1)]
    return -2*logsum(terms)/z**2

high, alpha = F(8,5), F(1)
rows = []
for low in [F(19,20), F(17,20)]:
    predicted = min(high-1, low-1/(1+2*alpha))
    observed = [(z, surrogate_rate(z,float(high),float(low),float(alpha)))
                for z in [12,48,192]]
    assert abs(observed[-1][1] - float(predicted)) < .0001
    rows.append(dict(low=str(low), predicted=str(predicted), rates=observed))

# Verify difference of the half-line kernels in the reflected-resolvent identity.
from math import sinh
errors = []
for x in [1.3,1.7,2.2]:
    y1,y2,k = 1.2,1.7,.8
    left = exp(-k*(y2-y1))*2*sinh(k*(x-y1))
    if x>y2:
        left -= 2*sinh(k*(x-y2))
    right = exp(-k*abs(x-y2))-exp(-k*(x+y2-2*y1))
    errors.append(abs(left-right))
assert max(errors)<1e-14

print(json.dumps(dict(python=platform.python_version(), sign_examples=rows,
                      max_kernel_error=max(errors)),indent=2))
