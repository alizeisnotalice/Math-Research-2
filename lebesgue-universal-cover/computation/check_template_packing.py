"""Integer consequences of the width-one packing construction."""
from fractions import Fraction as F
import json

def lower(m):
    # ceil(16*pi*m) <= ceil((352/7)*m); valid since pi<22/7.
    cap = (352*m+6)//7
    orbit_bound = 2*49**2*cap
    count = -(-(2**m)//orbit_bound)
    return dict(m=m, excess=str(F(1,32*m*m)), minimum_templates=count,
                orbit_bound=orbit_bound)

rows = [lower(m) for m in (64,128,176)]
assert rows[2]["minimum_templates"]>10**44
print(json.dumps(rows,indent=2))
