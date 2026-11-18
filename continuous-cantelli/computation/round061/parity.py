"""Verify parity and noninteger exponents in the signed monomial integral."""
from fractions import Fraction as F
from itertools import product
import json
n=0
for den,num in product(range(1,13),range(1,109)):
    alpha=F(num,den)
    if alpha<=2: continue
    for right,left in product(range(-2,3),repeat=2):
        if not (right or left): continue
        R,L=right/(1+alpha),-left/(1+alpha)
        matches=alpha.denominator==1 and L==(-1)**(alpha.numerator+1)*R
        assert matches==(alpha.denominator==1 and left==(-1)**alpha.numerator*right)
        n+=1
print(json.dumps(dict(denominators=12,largest_numerator=108,checked_pairs=n)))
