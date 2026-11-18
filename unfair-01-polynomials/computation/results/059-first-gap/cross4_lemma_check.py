"""Symbolic identities and exact rational regression for the general cross4 lemma."""
from fractions import Fraction as Fr
from itertools import product
from pathlib import Path
import json
import sys
import sympy as sp

x = sp.symbols('t')
m = 1-x+x*x
b = x*(2*x-1)+m
R = x*(1-x*m)+b*b/4
D = (1-x)*(5*x**3+3*x*(1-x)+3)
assert sp.expand(4*(1-R)-D)==0
assert sp.expand(m-b/2-(1-x*x)/2)==0
count = 0
# These exact grid checks exercise both endpoints of the A2 and A3 budget
# segments, and do not replace the proof over the real domain.
for k in range(1,41):
    t=Fr(k,41)
    A1,B1=t,1-t
    m2=1-t*(1-t)
    for u,v in product(range(17),repeat=2):
        A2=m2*u/16
        B2=m2-A2
        m3=1-A1*B2-B1*A2
        A3=m3*v/16
        B3=m3-A3
        interior=A1*B3+B1*A3+A2*B2
        assert 0<interior<1
        count+=1
out={'status':'verified','factorized_4_times_gap':str(sp.factor(4*(1-R))),
     'exact_grid_rows':count,'python':sys.version,'sympy':sp.__version__}
with Path(sys.argv[1]).open('x') as f:f.write(json.dumps(out,sort_keys=True,indent=2)+'\n')
print(json.dumps(out,sort_keys=True))
