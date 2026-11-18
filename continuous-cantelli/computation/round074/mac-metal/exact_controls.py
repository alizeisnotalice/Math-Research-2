"""Exact guards for a Cantor-distance exclusion and a digital Haar-test control."""
from fractions import Fraction as F
import json
from pathlib import Path
import sys
import mpmath as mp

assert F(1,6)<F(1,4)  # d^alpha <= sqrt(1/6) < 1/2 for alpha >= 3/4.
assert 3**7>F(46,7)**4
r_upper=F(7,23)
assert r_upper/(1-r_upper)==F(7,16)
# Main central first harmonic <= -(7/16)(h/2) I, with all other gaps strictly smaller.

p,q=7,11
B=2*p*q
even=[2*q*j for j in range(p)]
odd=[1+2*p*j for j in range(q)]
D=sorted(even+odd)
assert set(even).isdisjoint(odd)
assert all(n%p and n%q for n in range(1,7))
assert p%q and q%p and max(p,q)<B/2
# Each orbit sum is exactly zero unless its orbit length divides the frequency.
# At n=p,q the corresponding first factor survives, and later factors
# are confined to an open semicircle. These nonzero modes have gcd 1.
with mp.workdps(85):
    checks=[]
    for depth in (1,3,7):
        for n in (*range(1,7),p,q):
            P=mp.mpc(1)
            for j in range(1,depth+1):
                P*=sum(mp.exp(2j*mp.pi*n*d/B**j) for d in D)
            P*=(mp.exp(2j*mp.pi*n/B**depth)-1)/(2j*mp.pi*n)
            mag=abs(P)*(-mp.expm1(-mp.mpf(n*n)/2))
            if n<=6: assert mag<mp.mpf('1e-79')
            else: assert mag>mp.mpf('1e-12')
            checks.append({'depth':depth,'frequency':n,'magnitude':mp.nstr(mag,40)})
report={'middleThirdsDistance':{'minPower':'3/4','heightRange':'(0,1]','tailRadialMomentRatioUpper':str(r_upper/(1-r_upper))},
 'digitalStep':{'base':B,'digits':D,'zeroModes':list(range(1,7)),'detectingModes':[p,q],
 'mpmathValidation':checks},'mpmath':mp.__version__}
Path(sys.argv[1]).write_text(json.dumps(report,indent=2)+'\n')
print('exact guards passed; example signs/zeros numerically checked at three depths')
