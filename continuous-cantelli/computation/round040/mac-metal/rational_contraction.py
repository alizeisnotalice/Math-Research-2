"""Universal dyadic center-collapse/width-contraction obstruction.

For n >= 5 and a=2**n, set beta=(7/5)n+2. On x=a-z/a, 0<=z<=1,
require x/f(x) >= a+(beta-kappa*z)/a for any fixed 0<=kappa<=1,
or equivalently f(x) <= x/(a+(beta-kappa*z)/a).
Choose t in a+[-1/16,1/16]. The positive action is bounded below by
  beta/2 - beta/(16*a) - beta**2/(2*a**2) - (1/a+1/16)**2/2.
All three losses decrease with n, since beta(n+1)<=2 beta(n).
At n=5 they sum to 63/1024<1/16.
"""
from fractions import Fraction as F
import json
from math import factorial

a=F(32)
beta=F(7,5)*5+2
d=F(1,16)
loss=beta*d/a+beta**2/(2*a**2)+(1/a+d)**2/2
assert loss == F(63,1024) < F(1,16)
assert beta+F(7,5)<=2*beta

def taylor(q,degree):
    return sum(q**j/factorial(j) for j in range(degree+1))

assert taylor(F(7,10),3)>2
normal_lower=F(199,500)
assert normal_lower**2*F(44,7)<1  # pi < 22/7
ratio_lower=normal_lower*taylor(F(15,16),4)
assert ratio_lower==F(265681517,262144000)>1
print(json.dumps({'dyadicCenterExponentAtLeast':5,'betaSlope':'7/5','betaIntercept':2,
   'compressionRange':[0,1],'scaledSourceWidth':1,'observationHalfwidth':str(d),
   'lossUpperAtN5':str(loss),'strictDualRatioLower':str(ratio_lower)},indent=2))
