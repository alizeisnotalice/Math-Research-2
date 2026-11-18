"""Independent checks of the Stieltjes clock and null dimension-one heat model."""

import itertools
import json
import sys
from pathlib import Path

import mpmath as mp
import numpy as np
import sympy as sp


def exact():
    k,u,v = sp.symbols('k u v',positive=True)
    tail = sp.erfc(k/sp.sqrt(2*(1+u)))/2
    derivative = k*sp.exp(-k*k/(2*(1+u)))/(2*sp.sqrt(2*sp.pi)*(1+u)**sp.Rational(3,2))
    assert sp.simplify(sp.diff(tail,u)-derivative)==0
    th,x = sp.symbols('theta x',positive=True)
    assert sp.simplify(sp.diff(th*x**(1+th)/(1+th),x)-x*sp.diff(x**th,x))==0
    return str(derivative)


def kant(n,beta,age):
    # Restrict to the leftmost level-(n-m) ancestor; the omitted
    # heat tail is negligible because each scaled unit cell has mass <= 3.
    m,k=7,10
    offsets=np.array([0.0])
    ls=lambda j: 2.**(n-j)*((n+1)/(j+1))**beta
    for level in range(n-m,n+k):
        offsets=np.concatenate((offsets,offsets+ls(level)-ls(level+1)))
    length=ls(n+k)
    z=offsets+length/2
    heat=np.exp(-(z+1)**2/(2*age))/np.sqrt(2*np.pi*age)
    value=(2**m)*float(np.mean(heat))
    midpoint_abs_error=(2**m)*length**2/(8*np.sqrt(2*np.pi)*age**1.5)
    return value,midpoint_abs_error


def fractional(kappa):
    K=mp.mpf(1)
    half=mp.mpf('.5')
    power=10

    def left(y):
        v=half*y**power
        return power*half**(1-kappa)*y**(power*(1-kappa)-1)*(1-v)**(kappa-half)*mp.exp(-K*K*(1-v)/2)

    def right(y):
        w=half*y**power
        return power*half**(kappa+half)*y**(power*(kappa+half)-1)*(1-w)**(-kappa)*mp.exp(-K*K*w/2)

    via_quadrature=-K/(2*kappa*mp.sqrt(2*mp.pi))*(mp.quad(left,(0,1))+mp.quad(right,(0,1)))
    via_beta=-K/(2*kappa*mp.sqrt(2*mp.pi))*mp.exp(-K*K/2)*mp.beta(1-kappa,kappa+mp.mpf('.5'))*mp.hyp1f1(1-kappa,mp.mpf('1.5'),K*K/2)
    return via_quadrature,via_beta


def staircase_clock(L,theta):
    def g(x):
        return mp.sqrt(x)+mp.mpf('.8')*mp.sin(mp.sqrt(x))
    def gp(x):
        return (1+mp.mpf('.8')*mp.cos(mp.sqrt(x)))/(2*mp.sqrt(x))
    scaled=mp.quad(lambda w: theta*gp(L+w)*mp.exp(-2*w+theta*(g(L+w)-g(L))),(0,1,4,16,48,mp.inf))
    D=mp.exp(mp.mpf('.9'))
    bound=mp.expm1(theta*mp.log(D))/(1-mp.exp(-2)*D**theta)
    assert scaled<bound
    return scaled,bound


def main():
    derivative=exact()
    with mp.workdps(64):
        rows=[]
        for n,beta,age in itertools.product((16,64,256,1024),(1,2),(1,2,4,16)):
            v,e=kant(n,beta,age)
            target=mp.erfc(1/mp.sqrt(2*age))/2
            rows.append({'level':n,'beta':beta,'age_ratio':age,'normalized_heat_midpoint':v,
                         'midpoint_quadrature_error_upper_bound':e,'halfline_limit':mp.nstr(target,26),
                         'error_to_limit':mp.nstr(v-target,15)})
        frac=[]
        for kappa in (mp.mpf('.2'),mp.mpf('.5'),mp.mpf('.8')):
            numerical,closed=fractional(kappa)
            error=abs(numerical-closed)
            assert error<mp.mpf('1e-60'),error
            frac.append({'kappa':str(kappa),'limiting_fractional_difference':mp.nstr(closed,40),
                         'cross_method_error':mp.nstr(error,15)})
        clocks=[]
        for L,t in itertools.product((1,100,10000),('.05','.2')):
            v,b=staircase_clock(L,mp.mpf(t))
            clocks.append({'L':L,'theta':t,'short_clock_normalized':mp.nstr(v,36),
                           'unit_ratio_bound':mp.nstr(b,36)})
        out={'scope':'Cantor endpoint heat-kernel calibration, exact derivative and Stieltjes algebra; no exact variance graph claim.',
             'python':sys.version,'sympy':sp.__version__,'mpmath':mp.__version__,'numpy':np.__version__,
             'positive_heat_tail_derivative':derivative,
             'cantor_level_rows':rows,'negative_fractional_limits':frac,'short_clock_rows':clocks}
    p=Path(sys.argv[1]) if len(sys.argv)>1 else Path(__file__).with_name('gauge-results.json')
    p.write_text(json.dumps(out,indent=2)+'\n')
    print(len(rows),*[r['limiting_fractional_difference'] for r in frac])


if __name__=='__main__':
    main()
