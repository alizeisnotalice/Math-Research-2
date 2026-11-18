"""Exact combined lower bound from a critical-origin band and a distant patch."""
from fractions import Fraction as Q
from math import factorial
from pathlib import Path
import json
import platform


def exp_partial(x, n):
    return sum(x**i/factorial(i) for i in range(n+1))


def exp_upper(x, n):
    ratio=x/(n+2)
    assert ratio<1
    return exp_partial(x,n)+x**(n+1)/factorial(n+1)/(1-ratio)


t=Q(4096)
ylo,yhi=Q(10),Q(16)
slo,shi=Q(7,2),Q(9,2)
d0=Q(13,2)
# For ell=log(1/x) and x≤20/t, these inequalities certify
# ell + log(2 ell)/2 > 13/2 in the origin model.
assert exp_upper(Q(133,25),40)<t/20
assert exp_upper(Q(59,25),40)<2*Q(133,25)
assert Q(133,25)+Q(59,25)/2==d0
# On the origin band and observation window, the squared cost is <2/5.
distance=max(abs((1+slo/t)*ylo-2*d0),abs((1+shi/t)*yhi-2*d0))
assert distance**2/(4*d0)<Q(2,5)
assert exp_partial(d0,22)>600
# sqrt(2pi)<51/20 by pi<22/7 and exp(-x²/2)>1-x²/2.
assert Q(44,7)<Q(51,20)**2
assert (1-(yhi/t)**2/2)*Q(20,51)>Q(39,100)
origin=(yhi-ylo)/t*600*Q(3,5)*Q(39,100)
# For every observation s∈[t+7/2,t+9/2], the tail patch [t+2,t+6]
# contains its shifted interval s+[-3/2,3/2]. Upper Mills on both
# complementary Gaussian tails uses exp(9/8)>3 and sqrt(2pi)>12/5.
assert exp_partial(Q(9,8),4)>3
assert 4*sum(Q((-1)**i,2*i+1) for i in range(8))>3
assert Q(12,5)**2<6
tail=1-2/(3*Q(3,2)*Q(12,5))
assert origin+tail>1
# If these two source bands are saturated by equality, their
# contribution at the nominal t is below one. Other sources are not
# included in this upper bound.
assert exp_partial(Q(121,20),25)>t/ylo
assert exp_partial(Q(5,2),12)>Q(121,10)
assert Q(121,20)+Q(5,2)/2==Q(73,10)
assert exp_upper(Q(73,10),45)<1500
assert exp_partial(Q(2),5)>7
point_upper=(yhi-ylo)/t*1500*Q(5,12)+1/(7*2*Q(12,5))
assert point_upper<1
out=dict(reference_t=str(t),source_origin_band=[str(ylo/t),str(yhi/t)],local_D_lower=str(d0),
         origin_model='u^2/[2 (log(1/u)+log(2 log(1/u))/2)]',
         observation_offsets=[str(slo),str(shi)],source_tail_offsets=['2','6'],
         upper_tail_f='1',strict_origin_lower=str(origin),strict_tail_lower=str(tail),
         strict_combined_lower=str(origin+tail),rigorous_gap=str(origin+tail-1),
         strict_nominal_saturated_two_source_upper=str(point_upper),
         arithmetic='Fraction and exact exponential series bounds',python=platform.python_version())
Path('certificate.json').write_text(json.dumps(out,indent=2)+"\n")
print(json.dumps(out,indent=2))
