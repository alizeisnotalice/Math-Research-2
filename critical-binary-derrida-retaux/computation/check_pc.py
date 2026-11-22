"""Exact checks for the zero-mass convention; run with Python 3."""
from fractions import Fraction as F

def threshold(law):
    assert sum(law.values()) == 1
    a = sum((k-1)*2**k*w for k, w in law.items())
    return a/(1+a)

for law, expected in [
    ({2: F(1)}, F(4, 5)),
    ({3: F(1)}, F(16, 17)),
    ({1: F(1, 2), 2: F(1, 2)}, F(2, 3)),
]:
    p = threshold(law)
    assert p == expected
    for probe, sign in [(p/2, 1), (p, 0), ((1+p)/2, -1)]:
        delta = -probe + (1-probe)*sum((k-1)*2**k*v for k,v in law.items())
        assert (delta > 0) - (delta < 0) == sign
    print(law, p)

for r in (F(3,5), F(2,3), F(3,4), F(9,10)):
    x = 2*(1-r)
    a = 2*r*x/(1-x)**2
    p = a/(1+a)
    assert p == 4*r*(1-r)
    print("geometric", r, p)

# Deterministic stress test for the H3 truncation certificate.
cases = 0
for top in range(2, 65):
    for numerator in (1, 3, 7):
        mass = F(numerator, 8)
        law = {1: 1-mass, top: mass}
        exact = threshold(law)
        moment = sum(k**3*2**k*v for k,v in law.items())
        for cut in range(1, top+1):
            lower = sum(((k-1)*2**k*v for k,v in law.items() if k <= cut), F(0))
            upper = lower + moment/(cut+1)**2
            assert lower/(1+lower) <= exact <= upper/(1+upper)
            cases += 1
print("exact interval checks", cases)
