#!/usr/bin/env python3
"""Independent Fraction side-length/L1-margin certificate for reflection-column proof."""
from fractions import Fraction as F
from pathlib import Path
import hashlib,json,sys
lower2=F(353,250);upper2=F(3,2);upper6=F(49,20)
assert lower2**2<2<upper2**2 and upper6**2>6
assert 2*upper2-1==2
assert 3-lower2<2
assert 2*upper6-lower2-2==F(186,125)<2
assert (5*lower2-7)/3==F(1,50)
assert F(1,50)-F(12,1000)==F(1,125)
rows=[]
for N in range(2,7):
    K=2*N+2*(N-1)+(N-1)+(N-2)
    assert K==6*N-5
    loss=2*F(1,1000)*K/N
    assert loss<F(12,1000)
    rows.append(dict(N=N,K=K,ratio_loss=str(loss)))
p=Path(__file__).resolve()
out=dict(status='PASS',python=sys.version,arithmetic='Fraction, exact rational comparisons',random_seed=None,floating_tolerance=None,sqrt2_lower_squared=str(lower2**2),sqrt2_upper_squared=str(upper2**2),sqrt6_upper_squared=str(upper6**2),q2_side_sum_strict_upper=str(2*upper6-lower2-2),eta='1/1000',remaining_N3_gap_strict_lower='1/125',cases=rows,scope='Independent rational side-length and smoothing-margin checks, general-N statements proved analytically in review',script_sha256=hashlib.sha256(p.read_bytes()).hexdigest())
result=p.parent/'results/reflection-column-independent/margins.json'
result.parent.mkdir(parents=True,exist_ok=True)
result.write_text(json.dumps(out,ensure_ascii=False,indent=2)+'\n')
print(json.dumps(out,ensure_ascii=False,indent=2))
