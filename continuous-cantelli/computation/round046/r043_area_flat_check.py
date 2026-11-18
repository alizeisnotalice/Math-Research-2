"""Exact checks of interpolated exponents and narrow-spike log-area scales."""
from fractions import Fraction as F
import json, platform
import mpmath as mp

rows=[]
for s in (F(1,2),F(3,4),F(9,10),F(99,100)):
    n=3
    while F(1)-F(2,n)<=s:
        n+=1
    t=1-F(2,n)
    # (u^2/m)^t (u^-(N-2))^(1-t) loses all u powers.
    assert 2*t-(n-2)*(1-t)==0
    theta=s/t
    exponent=(n-2)*(1-theta)
    assert theta*t==s and exponent>0
    rows.append({'s':str(s),'N':n,'t':str(t),'remaining_distance_power':str(exponent)})
mp.mp.dps=60
gamma=mp.mpf('1.5')
spikes=[]
for j in (10,20,40,80):
    m=mp.power(2,-j)
    # Dominating triangular bump has exactly half its area inside its first
    # radius; all further areas are bounded by m exp(-4^{j+1}) .
    log_main=mp.log(m/2)-4**j
    ratio=mp.power(m,-2)/mp.power(m,-mp.mpf('.9'))
    spikes.append({'j':j,'log_half_current_spike_area':mp.nstr(log_main,17),
       'spike_action_over_height_to_minus_0.9':mp.nstr(ratio,15),
       'log_modulus_action_over_height_to_minus_0.9':mp.nstr((2/m)**(1/gamma)/(m**(-mp.mpf('.9'))),15)})
print(json.dumps({'python':platform.python_version(),'mpmath':mp.__version__,'digits':mp.mp.dps,
    'exact_exponents':rows,'Gaussian_half_window_probability':mp.nstr(mp.erf(1/mp.sqrt(2))/2,24),
    'spikes':spikes},indent=2))
