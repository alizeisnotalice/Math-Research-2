"""Join the disjoint new counts with the refined modular denominator."""
from decimal import Decimal, ROUND_FLOOR, localcontext
from fractions import Fraction
from hashlib import sha256
from math import factorial
from pathlib import Path
import json

D=Path(__file__).resolve().parent
R51=D.parents[1]/'r051/n12'
den=json.loads((D/'r052-denominator-refined.json').read_text())
num=json.loads((D.parent/'construction/result.json').read_text())
old_path=R51/'n12-conditional-final.json'
signature=sha256(old_path.read_bytes()).hexdigest()
assert den['r051_numerator_sha256']==num['source_R51_sha256']==signature
old=json.loads(old_path.read_text())
N=int(den['denominator'])
rows=[]
for a,b in zip(num['thresholds'],old['thresholds'],strict=True):
    cutoff=a['strict_threshold']
    assert cutoff==b['strict_threshold']
    old_count=int(b['sum_parts']['total'])
    add_count=sum(int(z['additive_count']) for z in num['new_strata'] if z['distance']>cutoff)
    new_count=old_count+add_count
    assert old_count==int(a['old_R51']) and new_count==int(a['combined'])
    value=Fraction(new_count,N)
    assert 0<value<=1
    with localcontext() as ctx:
        ctx.prec=25
        ctx.rounding=ROUND_FLOOR
        display=str(Decimal(value.numerator)/Decimal(value.denominator))
    rows.append({'strict_threshold':cutoff,'combined_numerator':str(new_count),
                 'lower_exact':str(value),'lower_decimal_rounded_down':display})
upper=Fraction(N,factorial(144))
print(json.dumps({'chosen_layout':den['selected_layout'],'A2_upper':str(N),
                  'A2_over_144factorial_upper_exact':str(upper),'tail_lower_bounds':rows},indent=2))
