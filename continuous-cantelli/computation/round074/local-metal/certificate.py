"""Exact real-arithmetic error bounds for the Cantor level calculation."""
from fractions import Fraction as Q
from math import factorial


def exp_lower(x, terms):
    return sum(x**j/factorial(j) for j in range(terms+1))


# Uniform geometric upper bound at depth 96. At fixed lambda it decreases
# with p, so use p=1/2. At that exponent the bound increases up to lambda=.461.
lam=Q(461,1000)
assert 144-291*lam>0
assert Q(68,100)**2>lam
assert Q(198,1000)**2>(1-2*lam)/2
assert Q(22,7)**2<10
ratio=2*lam*Q(68,100)
geometric=20*1024**2*32*(1-2*lam)*Q(198,1000)*ratio**96/(Q(3,2)*(1-ratio))
assert geometric<Q(1,10**12)
fine=96*Q(2,10**14)
heat=1/exp_lower(Q(42),100)
expm1=Q(6,100)**7/factorial(7)
sinc=Q(1,100)**6/factorial(7)
assert geometric+fine+heat+expm1+sinc<Q(4,10**12)
for label,value in (("geometric",geometric),("fine",fine),("heat",heat),
                    ("expm1",expm1),("sinc",sinc),
                    ("total_before_quadrature",geometric+fine+heat+expm1+sinc)):
    print(label,float(value))
