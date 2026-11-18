"""Independent Float64 log-CDF and high-precision verification of window scans."""
import json
import math
from pathlib import Path
import sys

import mpmath as mp
import numpy as np
import scipy
from scipy.integrate import quad
from scipy.optimize import brentq
from scipy.special import log_ndtr


def get(p):
    return [p[k] for k in ('a','w','r','length','b')]


def constants(a,w,r,L,b):
    rt=math.sqrt(2)
    phi=lambda y: math.exp(-y*y/2)/math.sqrt(2*math.pi)
    pv=math.erf(w/(a*rt))
    pb=math.erf((r+L)/rt)-math.erf(r/rt)
    A=b*pb/pv
    m2=1+2*(r*phi(r)-(r+L)*phi(r+L))/pb
    D=(1+1/A)*math.log1p(A)+(1-b)/b*math.log1p(-b)
    return A,m2,D


def logprob(lo,hi,shift):
    x,y=lo-shift,hi-shift
    if x>=0:
        p,q=log_ndtr(-x),log_ndtr(-y)
    else:
        p,q=log_ndtr(y),log_ndtr(x)
    return p+math.log(-math.expm1(q-p))


def double_fun(p, derivative):
    a,w,r,L,b=get(p)
    A,m2,D=constants(a,w,r,L,b)
    scale=D if derivative else math.log(a)
    v=w/a
    rt=math.sqrt(2)
    pv=math.erf(v/rt)
    pb=math.erf((r+L)/rt)-math.erf(r/rt)
    norm=math.sqrt(2/math.pi)
    def evaluate(nu):
        s=nu*scale/(a*m2)
        if abs(s)<=4:
            shift=math.exp(-s*s/2)
            delta=math.expm1(-s*s/2)
            def shape(y):
                return math.exp(-y*y/2)*(delta+2*shift*math.sinh(s*y/2)**2)
            dv=norm*quad(shape,0,v,epsabs=1e-24,epsrel=1e-12)[0]/pv
            db=norm*quad(shape,r,r+L,epsabs=1e-24,epsrel=1e-12)[0]/pb
            ds=(s/a)*(2+s/a)
            ga=(1+1/A)*math.expm1(ds*math.log1p(A))
            gb=-(1-b)/b*math.expm1(ds*math.log1p(-b))
            raw=(ga+dv+ga*dv)-(gb+db+gb*db)
            return raw/((1+ga)*(1+dv)+(1+gb)*(1+db))
        q=(1+s/a)**2
        gain=q*math.log1p(A)
        loss=q*math.log1p(-b)
        if q==0:
            return 0.0
        lp=gain+math.log(-math.expm1(-gain))+logprob(-v,v,s)
        ln=math.log(-math.expm1(loss))+float(np.logaddexp(logprob(r,r+L,s),logprob(-r-L,-r,s)))
        return math.tanh((lp-ln)/2)
    return evaluate,scale,m2,D


def mp_fun(p, derivative):
    a,w,r,L,b=map(lambda x:mp.mpf(str(x)),get(p))
    rt=mp.sqrt(2)
    norm=mp.sqrt(2*mp.pi)
    pv=mp.erf(w/(a*rt))
    pb=mp.erf((r+L)/rt)-mp.erf(r/rt)
    A=b*pb/pv
    phi=lambda y:mp.exp(-y*y/2)/norm
    m2=1+2*(r*phi(r)-(r+L)*phi(r+L))/pb
    D=(1+1/A)*mp.log1p(A)+(1-b)/b*mp.log1p(-b)
    scale=D if derivative else mp.log(a)
    v=w/a
    def prob(lo,hi,s):
        x,y=lo-s,hi-s
        if x>=0:
            return (mp.erfc(x/rt)-mp.erfc(y/rt))/2
        return (mp.erfc(-y/rt)-mp.erfc(-x/rt))/2
    def evaluate(nu):
        s=nu*scale/(a*m2)
        q=(1+s/a)**2
        if q==0:
            return mp.mpf(0)
        pos=mp.expm1(q*mp.log1p(A))*prob(-v,v,s)
        neg=-mp.expm1(q*mp.log1p(-b))*(prob(r,r+L,s)+prob(-r-L,-r,s))
        return (pos-neg)/(pos+neg)
    return evaluate


def check(path, derivative):
    data=json.loads(Path(path).read_text())
    results={}
    for key,hit in data['extremalCrossings'].items():
        p=hit['profile']
        lo,hi=hit['lowerNu'],hit['upperNu']
        fun,scale,m2,D=double_fun(p,derivative)
        original=[fun(lo),fun(hi)]
        step=(hi-lo)/3
        left=max(lo-step,lo/2,1e-12)
        right=hi+step
        if fun(left)<=0 or fun(right)>=0:
            with mp.workdps(80):
                mf=mp_fun(p,derivative)
                endpoint=[mp.nstr(mf(mp.mpf(str(x))),30) for x in (lo,hi)]
            results[key]={'gpuBracket':[lo,hi],'doubleSignRatiosAtGpuBracket':original,
                          'mpmathSignRatiosAtGpuBracket':endpoint,'literalGpuBracketConfirmed':False,
                          'profile':p,'rejectedBracket':True}
            print(f'{path}/{key}: reference rejected GPU bracket {lo,hi}',flush=True)
            continue
        root=brentq(fun,left,right,xtol=1e-12)
        with mp.workdps(80):
            mf=mp_fun(p,derivative)
            high=mp.findroot(mf,(mp.mpf(str(root)),mp.mpf(str(root+0.0001))),tol=mp.mpf('1e-66'))
            assert abs(float(high)-root)<1e-5
            rec={'gpuBracket':[lo,hi],'doubleSignRatiosAtGpuBracket':original,
                 'literalGpuBracketConfirmed':original[0]>=0 and original[1]<=0,
                 'doubleRoot':root,'mpmathRoot':mp.nstr(high,55),
                 'rootInDerivativeScale':float(high)*scale/D,
                 'shift':float(high)*scale/(p['a']*m2),'profile':p}
            results[key]=rec
            print(f'{path}/{key}: ν={float(high):.12g}, D-scaled={rec["rootInDerivativeScale"]:.12g}',flush=True)
    return results


def main():
    directory=Path(sys.argv[1])
    runs=[('extension.json',False),('derivative/derivative.json',True),('finite/finite.json',True),('far_repair/far_repair.json',True),('large_far_repair/large_far_repair.json',True)]
    checks={k:check(directory/k,d) for k,d in runs if (directory/k).exists()}
    stress={'a':8,'w':1/128,'r':0.25,'length':0.125,'b':0.98}
    fun,scale,m2,D=double_fun(stress,True)
    with mp.workdps(85):
        mf=mp_fun(stress,True)
        obj=[]
        for nu in (-2,1,4,8,12,16):
            v=mf(mp.mpf(nu))
            assert abs(float(v)-fun(nu))<1e-9
            obj.append({'nu':nu,'physicalShift':nu*scale/(stress['a']*m2),'signedRelative':mp.nstr(v,40)})
    out={'python':sys.version,'scipy':scipy.__version__,'mpmath':mp.__version__,'runs':checks,'largeShiftStress':{'profile':stress,'rows':obj}}
    Path(sys.argv[2]).write_text(json.dumps(out,indent=2)+'\n')


if __name__=='__main__': main()
