"""Reference call residuals for exact fourth-moment matched cosine profiles."""
import math
from pathlib import Path
import json
import sys

import mpmath as mp
import numpy as np
import scipy
from scipy.integrate import quad
from scipy.optimize import brentq, minimize_scalar
from scipy.special import ndtr, erfcx, roots_legendre

norm=math.sqrt(2*math.pi)
u,w=roots_legendre(16)
u=(u+1)/2
w/=2


def setup(c):
    k,l,z,m=[c[x] for x in ('k','l','coefficient','mean')]
    rk,rl=math.exp(-k*k/2),math.exp(-l*l/2)
    cov=(math.exp(-(k-l)**2/2)+math.exp(-(k+l)**2/2))/2-rk*rl
    V=(1-rk*rk)**2/2+z*z*(1-rl*rl)**2/2+2*z*cov
    L=2*(k*k*rk+z*l*l*rl)/V
    center=rk+z*rl
    return k,l,z,m,L,center,V


def funs(c):
    k,l,z,m,L,center,V=setup(c)
    target=math.sqrt(1+m)
    def derivative(y):
        a=y*target
        def inside(x):
            dv=L*(math.cos(k*x)+z*math.cos(l*x)-center)
            noise=m+u*dv
            local=(np.exp(-0.5*(a-x)**2/noise)+np.exp(-0.5*(a+x)**2/noise))/np.sqrt(noise)
            return math.exp(-x*x/2)*dv*float(np.dot(w,local))
        return quad(inside,0,12,epsabs=2e-14,epsrel=5e-12,limit=300)[0]/(4*math.pi*target)
    def call(v,x,a):
        s=math.sqrt(v)
        b=(x-a)/s
        g=math.exp(-b*b/2)
        value=g/norm+b*ndtr(b) if b>=0 else g*(1/norm+b*erfcx(-b/math.sqrt(2))/2)
        return s*value
    def direct(y):
        a=y*target
        def inside(x):
            v=m+L*(math.cos(k*x)+z*math.cos(l*x)-center)
            d=call(v,x,a)-call(m,x,a)+call(v,-x,a)-call(m,-x,a)
            return math.exp(-x*x/2)*d
        return quad(inside,0,12,epsabs=2e-14,epsrel=5e-12,limit=300)[0]/(norm*target)
    return derivative,direct


def moments(c):
    k,l,z,m,L,center,V=setup(c)
    one=[(0.,-L*center),(k,L/2),(-k,L/2),(l,z*L/2),(-l,z*L/2)]
    two=[(a+b,u*v) for a,u in one for b,v in one]
    three=[(a+b,u*v) for a,u in two for b,v in one]
    expect=lambda poly,j:sum(w*(-t*t)**j*math.exp(-t*t/2) for t,w in poly)
    return {'varianceNoise':L*L*V,'meanVariancePerturbation':expect(one,0),
            'standardizedHermite4':(6*expect(one,1)+3*expect(two,0))/(1+m)**2,
            'standardizedHermite6':15*(expect(one,2)+3*expect(two,1)+expect(three,0))/(1+m)**3}


def high(c,y):
    with mp.workdps(60):
        k,l,z,m=map(lambda key:mp.mpf(str(c[key])),('k','l','coefficient','mean'))
        rk,rl=mp.exp(-k*k/2),mp.exp(-l*l/2)
        cov=(mp.exp(-(k-l)**2/2)+mp.exp(-(k+l)**2/2))/2-rk*rl
        V=(1-rk*rk)**2/2+z*z*(1-rl*rl)**2/2+2*z*cov
        L=2*(k*k*rk+z*l*l*rl)/V
        center=rk+z*rl
        target=mp.sqrt(1+m)
        strike=mp.mpf(str(y))*target
        nf=mp.sqrt(2*mp.pi)
        r2=mp.sqrt(2)
        def call(v,x):
            s=mp.sqrt(v)
            b=(x-strike)/s
            return s*(mp.exp(-b*b/2)/nf+b*mp.erfc(-b/r2)/2)
        def inside(x):
            v=m+L*(mp.cos(k*x)+z*mp.cos(l*x)-center)
            d=call(v,x)-call(m,x)+call(v,-x)-call(m,-x)
            return mp.exp(-x*x/2)*d/nf
        return mp.quad(inside,list(range(15)))/target


def main():
    path=Path(sys.argv[1]);dest=Path(sys.argv[2])
    payload=json.loads(path.read_text())
    verified={}
    dim=payload['dimensions'][4]
    strikes=np.linspace(0,7,dim)
    for label,item in payload['representatives'].items():
        core=item['core']
        f,g=funs(core)
        diag={}
        for tag in ('positive','negative'):
            strike=item[tag+'Strike']
            delta=f(strike)
            check=g(strike)
            assert abs(delta-check)<1e-11
            optimized=minimize_scalar(lambda x:(-1 if tag=='positive' else 1)*f(x),
                   bounds=(max(0,strike-7/(dim-1)),min(7,strike+7/(dim-1))),method='bounded',
                   options={'xatol':1e-10})
            optimized_at=max((float(optimized.x),strike),key=lambda t:(1 if tag=='positive' else -1)*f(t))
            diag[tag]={'atMetalStrike':delta,'directConditionalCallCheck':check,
                       'refinedStrike':optimized_at,'refinedValue':f(optimized_at)}
        vals=[f(float(t)) for t in strikes]
        threshold=max(1e-6*max(abs(x) for x in vals),1e-12)
        visible=[(float(t),v) for t,v in zip(strikes,vals) if abs(v)>=threshold]
        roots=[brentq(f,left,right,xtol=1e-10) for (left,v),(right,w) in zip(visible,visible[1:]) if v*w<0]
        diag['visibleSignChangeStrikes']=roots
        diag['moments']=moments(core)
        if label in ('leastScaledResidual','mostSignChanges'):
            diag['mp60']=mp.nstr(high(core,item['positiveStrike']),45)
            assert abs(float(diag['mp60'])-diag['positive']['atMetalStrike'])<1e-12
        verified[label]={'metal':item,**diag}
        print(label,'roots',roots,'positive',diag['positive']['atMetalStrike'],
            'negative',diag['negative']['atMetalStrike'],'H6',diag['moments']['standardizedHermite6'],flush=True)
    dest.write_text(json.dumps({'python':sys.version,'scipy':scipy.__version__,'mpmath':mp.__version__,'verified':verified},indent=2)+'\n')


if __name__=='__main__': main()
