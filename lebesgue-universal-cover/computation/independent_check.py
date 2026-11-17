"""Small independent standard-library check, suitable for python3 - over SSH."""
from fractions import Fraction as F
from itertools import combinations
import json, sys, math
V=[(F(1,2),F(1,4)),(F(1,6),F(1,2)),(-F(1,6),F(1,2)),(-F(1,2),F(1,4)),(-F(1,2),-F(1,4)),(-F(1,6),-F(1,2)),(F(1,6),-F(1,2)),(F(1,2),-F(1,4))]
area=F(0);perimeter=F(0)
for a,b in zip(V,V[1:]+V[:1]):
    area+=(a[0]*b[1]-a[1]*b[0])/2
    q=sum((x-y)**2 for x,y in zip(a,b))
    z=F(math.isqrt(q.numerator),math.isqrt(q.denominator));assert z*z==q
    perimeter+=z
assert area==F(5,6) and perimeter==F(10,3) and 4*area==perimeter
for n in range(1,101):
    h=F(1,n);assert 2*h*(1+h)/(1+2*h)<=2*h
print(json.dumps({'python':sys.version,'area':str(area),'perimeter':str(perimeter),'rational_scaling_checks':100,'status':'passed','scope':'one tangential octagon and finite scalar checks only'},indent=2))
