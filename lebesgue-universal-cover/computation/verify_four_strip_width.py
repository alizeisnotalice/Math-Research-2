"""Exact squared radii of the four-direction signed Fourier sums in Q(sqrt(2))."""
from collections import Counter
from itertools import product
import json
import platform

# Directions exp(-3 i j*pi/4). Write 2*Re and 2*Im as a+b sqrt(2).
histogram = Counter()
for v0,v1,v2,v3 in product((-1,1),repeat=4):
    x,y = (2*v0,-v1+v3),(2*v2,-v1-v3)
    rational = x[0]**2+y[0]**2+2*(x[1]**2+y[1]**2)
    radical = 2*(x[0]*x[1]+y[0]*y[1])
    assert rational % 4 == 0 and radical % 4 == 0
    histogram[rational//4,radical//4] += 1
assert histogram == {(4,-2):8,(4,2):8}
print(json.dumps({"python":platform.python_version(),"squared_moduli":[
    {"rational":a,"sqrt2_coefficient":b,"multiplicity":n}
    for (a,b),n in sorted(histogram.items())],"passed":True},indent=2))
