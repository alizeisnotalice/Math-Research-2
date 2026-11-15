#!/usr/bin/env python3
"""Independent exact audit of the R56 shared-background transport obstruction.
No numerical result is used as a substitute for the accompanying proof.
"""
from fractions import Fraction as F
import hashlib, json, math, platform, sys
from pathlib import Path

def overlap(a,b,c,d):
    return max(F(0), min(b,d)-max(a,c))

def mass(L,s,t,r):
    bg=overlap(s-r,s+r,-L,L)*overlap(t-r,t+r,0,L)
    excess=overlap(s-r,s+r,-1,1)*overlap(t-r,t+r,0,1)
    return bg+excess

# Exact polynomial identities. q=1/t^2 in (0,1/16].
# (4(t+1/t)^2+2)/t^2 = 4+10q+4q^2.
q=F(1,16)
assert 4+10*q+4*q*q==F(297,64)<5
assert 10+4*q==F(41,4)
assert F(3,2)/(4*F(297,64))==F(8,99)
assert F(3,2)/20==F(3,40)
assert (64+F(41,4))/(64+2)==F(9,8)
checks=0
band_data=[]
for n in range(1,13):
    L=F(16*2**n)
    bands=[]
    for k in range(n):
        T=F(4*2**k)
        assert 2*T<=L/4
        bands.append({'T':str(T),'area':str(F(3,2)*T*T),
                      'denominator_bound':str(F(297,16)*T*T),
                      'load_lower_bound':str(F(8,99))})
        for t in [T,F(5,4)*T,F(3,2)*T,F(7,4)*T,2*T]:
            for s in [-t/2,-t/4,F(0),t/4,t/2]:
                assert mass(L,s,t,t)==4*t*t+2
                # Every radius t-1<r<t+1/t contains the whole S.
                for r in [t-F(1,2),t,t+F(1,2)/t]:
                    assert r>t-1
                    assert s-r<=-1 and s+r>=1
                    assert t-r<1 and t+r>2
                    assert mass(L,s,t,r)<=4*r*r+2
                    if mass(L,s,t,r)>4*r*r:
                        assert r<t+1/t
                        assert mass(L,s,t,r)<=(4*t*t+F(41,4))
                        assert mass(L,s,t,r)<=F(297,64)*t*t
                    checks+=1
    band_data.append({'n':n,'L':str(L),'rational_total_lower':str(F(8*n,99)),
                     'simple_total_lower':str(F(3*n,40)),
                     'J0':math.log((float(L*L/4)+2)/66)/8,
                     'J1':math.log((float(L*L/4)+10.25)/74.25)/8,
                     'Jstar':(lambda b,a: math.log(b/a)/4+(b**-4-a**-4)/16)
                         ((float(L)/4+math.sqrt(float(L*L)/16+4))/2,2+math.sqrt(5))})
# Rational parametrization audits the exact restricted-Omega minimax formula.
# Put t=R-1/R; then R(t)=R without any floating square root.
minimax_checks=0
for R in [F(5),F(9,2),F(6),F(17,3),F(20),F(100)]:
    t=R-1/R
    assert t>=4 and R*(R-t)==1
    L=8*t
    s=t/2
    assert 2*R*(R+t)+2==4*R*R
    # Exact chain-rule identity for the proposed antiderivative.
    assert (t/(4*R*R))*(1+1/(R*R))==F(1,4)*(1/R-1/R**5)
    old=None
    for m in [1,2,3,10,100,10000]:
        r=t+(1-F(1,m))*(R-t)
        M=mass(L,s,t,r)
        assert t<=r<R and M>4*r*r
        assert M==2*r*(r+t)+2 and M<4*R*R
        if old is not None: assert M>old
        old=M
        minimax_checks+=1
# Exact interval length shows witness strictness is essential but grazing excess
# does not avoid S; test r=t-1+1/j for j arbitrarily large along finite audit.
for j in [2,3,10,100,10**6]:
    t=F(4); r=t-1+F(1,j)
    assert -r < -1 and r>1 and t-r<1 and t+r>2
    checks+=1
result={'status':'PASS','python':sys.version,'platform':platform.platform(),
        'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        'exact_sample_checks':checks,'exact_minimax_checks':minimax_checks,'exact_symbolic_coefficients':
        {'quadratic_q_coefficients':[4,10,4],'q_upper_bound':'1/16',
         'normalized_mass_upper':'297/64','band_load_lower':'8/99',
         'canonical_minus_optimal_gap_upper_ratio':'9/8'},
        'scope':'Exact rational finite audits plus separately written analytic all-parameter proof; no global C2 bound.',
        'bands':band_data,'floating_gap_limit':math.log(9/8)/8,
        'randomness':'none','tolerance':'none for rational assertions; J0/J1 float diagnostics only'}
out=Path(__file__).with_name('shared-background-proportional-independent-output.json')
out.write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,indent=2))
