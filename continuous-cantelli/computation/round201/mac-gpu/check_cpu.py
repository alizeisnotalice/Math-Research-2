"""Multiprecision Kummer roots and independent Euler/time Marchaud integrals."""
import json
import math
from pathlib import Path
import platform
import struct
import sys

import mpmath as mp
from scipy.optimize import brentq
from scipy.special import hyp1f1
import scipy

mp.mp.dps=108

def f32(x):
    return struct.unpack('f',struct.pack('f',x))[0]

def root(k,x0):
    # Factoring k outside the hypergeometric function preserves it even
    # when the generic implementation rounds the small parameter -k to 0.
    fun=lambda x:1-2*k*x*mp.hyper([1,1-k],[2,mp.mpf('1.5')],x)
    deriv=lambda x:-2*k*mp.hyp1f1(1-k,mp.mpf('1.5'),x)
    return mp.findroot(fun,x0,df=deriv,solver='newton',tol=mp.mpf('1e-103'))

def regularized_integral(k,fn,constant):
    # Integrate [fn(u)-constant] on the left. The removed constant is for
    # the kernel including (1-u)^(k-.5), so its left integral is explicit.
    e=1-k
    power=k+mp.mpf('.5')
    h=mp.mpf('.5')
    def left(t):
        if not t:return 0
        return t**(-k)*((1-t)**(k-mp.mpf('.5'))*fn(t)-constant)
    def right(v):
        t=1-h*v**(1/power)
        return h**power/power*t**(-k)*fn(t)
    return constant*h**e/e+mp.quad(left,[0,mp.mpf('.125'),mp.mpf('.375'),h])+mp.quad(right,[0,mp.mpf('.25'),mp.mpf('.75'),1])

def beta_positive_residual(x,k):
    e=1-k
    def j_minus(t):
        if t==0:return 0
        return mp.expm1(x*t)/t-x
    val=regularized_integral(k,j_minus,0)
    residual=2*k*(x+val/mp.beta(e,k+mp.mpf('.5')))-1
    return residual/e if k>mp.mpf('.7') else residual

def time_residual(x,k):
    e=1-k
    def ff(u):
        return mp.exp(x*u)*(1-2*x*(1-u))
    return regularized_integral(k,ff,1-2*x)/mp.beta(e,k+mp.mpf('.5'))

def near_one():
    x=mp.mpf('.5')
    term=2*x*x/3
    total=term
    for n in range(2,100):
        term*=((n-1)*x)/((n+mp.mpf('.5'))*(n+1))
        total+=term
    return (1-total)/2

def main():
    path,out=map(Path,sys.argv[1:3]); raw=json.loads(path.read_text())
    reports=[]
    for bucket in raw['rows']:
        for hit in bucket['selected']:
            if hit['nearOne']:
                e=mp.exp(mp.mpf(f32(hit['logE'])))
                k=1-e
            else:
                k=mp.exp(mp.mpf(f32(hit['logK'])))
                e=1-k
            x=root(k,mp.mpf(f32(hit['rootX'])))
            y=mp.sqrt(2*x)
            stable=2*(x-mp.mpf('.5'))/e/(y+1)
            by={'group':bucket['group'],'index':hit['index'],'input':hit,
                'kappa':mp.nstr(k,54),'oneMinusKappa':mp.nstr(e,54),
                'rootX':mp.nstr(x,54),'rootY':mp.nstr(y,54),'scaledYMinusOne':mp.nstr(stable,54),
                'absoluteMetalRootYError':float(abs(y-f32(hit['rootY']))),
                'absoluteMetalScaledGapError':float(abs(stable-f32(hit['scaledOneGap'])))}
            if k>mp.mpf('1e-70'):
                kummer=mp.hyp1f1(-k,mp.mpf('.5'),x)
                assert abs(kummer)<mp.mpf('1e-68')
                by['directKummerResidual']=mp.nstr(kummer,24)
            if hit['group']==0 or hit['index'] in (0,raw['perStratum']-1):
                numerical=beta_positive_residual(x,k)
                assert abs(numerical)<mp.mpf('1e-68'),(hit,numerical)
                by['normalizedBetaIntegralResidual']=mp.nstr(numerical,24)
                if k>mp.mpf('.009') and e>mp.mpf('1e-8'):
                    v=time_residual(x,k)
                    assert abs(v)<mp.mpf('1e-68'),(hit,v)
                    by['normalizedTimeIntegralResidual']=mp.nstr(v,24)
            if k>mp.mpf('1e-6') and e>mp.mpf('1e-6'):
                double=brentq(lambda v:hyp1f1(-float(k),.5,v),.5,max(2.,f32(hit['rootX'])+2.),xtol=2e-14)
                by['scipyAbsoluteRootYError']=float(abs(mp.sqrt(2*double)-y))
            reports.append(by)
    exact=[]
    for text in ('.01','.1','.25','.5','.75','.9','.99'):
        k=mp.mpf(text); ell=-mp.log(k); x=root(k,max(mp.mpf('.5'),ell+mp.log(ell/mp.pi)/2))
        exact.append({'kappa':text,'rootY':mp.nstr(mp.sqrt(2*x),61),
            'timeResidual':mp.nstr(time_residual(x,k),25)})
    far=[]
    for ell in (16,64,128,256,600,1200):
        l=mp.mpf(ell)
        k=mp.exp(-l); q=mp.log(l/mp.pi)
        x=root(k,l+q/2)
        b0=2*l+q; b1=b0+(q-2)/(2*l)
        far.append({'minusLogKappa':ell,'rootY':mp.nstr(mp.sqrt(2*x),42),
                    'leadingY2Error':mp.nstr(2*x-b0,30),'refinedY2Error':mp.nstr(2*x-b1,30)})
    out.write_text(json.dumps({'python':platform.python_version(),'mpmath':mp.__version__,'scipy':scipy.__version__,
        'nearOneScaledGapLimit':mp.nstr(near_one(),61),'exactAnchors':exact,'smallKappaAsymptotics':far,'checked':reports},indent=2)+'\n')
if __name__=='__main__':main()
