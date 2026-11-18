"""Far return crossing check using stable Gaussian endpoint probabilities."""
import json
import math
from pathlib import Path
import sys

import mpmath as mp
import scipy
from scipy.optimize import brentq
from scipy.special import erfcx


def cn(p):
    return [p[k] for k in ('a','w','r','length','b')]


def parameter_double(p):
    a,w,r,L,b=cn(p)
    v=w/a
    pb=math.erf((r+L)/math.sqrt(2))-math.erf(r/math.sqrt(2))
    ell=math.log1p(b*pb/math.erf(v/math.sqrt(2)))
    return a,v,r,L,b,ell


def double_fun(p):
    a,v,r,L,b,ell=parameter_double(p)
    R=r+L
    result_scale=R*a*a/ell
    def relative_mass(s,u,d):
        su=s-u
        lead=math.log(erfcx(su/math.sqrt(2)))
        minus=math.log(erfcx((su+d)/math.sqrt(2)))
        return lead+math.log(-math.expm1(-su*d-d*d/2+minus-lead))

    def fn(factor):
        s=factor*result_scale
        q=(1+s/a)**2
        g=q*ell
        pos=relative_mass(s,v,2*v)
        band=relative_mass(s,R,L)
        neg_band=relative_mass(s,-r,L)
        log_reverse=-s*(R+r)+(R*R-r*r)/2+neg_band-band
        return (g+math.log(-math.expm1(-g))-math.log(-math.expm1(q*math.log1p(-b)))
                -s*(R-v)+(R*R-v*v)/2+pos-band-math.log1p(math.exp(log_reverse)))/(s*R)
    return fn


def mp_fun(p):
    a,w,r,L,b=map(lambda x:mp.mpf(str(x)),cn(p))
    v=w/a
    R=r+L
    rt=mp.sqrt(2)
    pv=mp.erf(v/rt)
    pb=mp.erf(R/rt)-mp.erf(r/rt)
    ell=mp.log1p(b*pb/pv)
    result_scale=R*a*a/ell
    def relative_mass(s,u,d):
        z=(s-u)/rt
        zs=(s-u+d)/rt
        lead=mp.log(mp.erfc(z))+z*z
        minus=mp.log(mp.erfc(zs))+zs*zs
        return lead+mp.log(-mp.expm1(-(s-u)*d-d*d/2+minus-lead))

    def fn(factor):
        s=factor*result_scale
        q=(1+s/a)**2
        g=q*ell
        pos=relative_mass(s,v,2*v)
        band=relative_mass(s,R,L)
        neg_band=relative_mass(s,-r,L)
        log_reverse=-s*(R+r)+(R*R-r*r)/2+neg_band-band
        return (g+mp.log(-mp.expm1(-g))-mp.log(-mp.expm1(q*mp.log1p(-b)))
                -s*(R-v)+(R*R-v*v)/2+pos-band-mp.log1p(mp.exp(log_reverse)))/(s*R)
    return fn,a,v,R,ell,result_scale


def solve(p,lo,hi):
    d=double_fun(p)
    # Roots close to a grid edge require extra room for Float32 error.
    x=brentq(d,max(0.1,lo-0.02),hi+0.02,xtol=1e-13)
    with mp.workdps(95):
        f,a,v,R,ell,scale=mp_fun(p)
        root=mp.findroot(f,(mp.mpf(str(x)),mp.mpf(str(x+0.0001))),tol=mp.mpf('1e-79'))
        assert abs(float(root)-x)<1e-10
        k=R-v-2*ell/a
        C=ell+(R*R-v*v)/2
        discr=k*k-4*ell*C/a**2
        proxy=a*a*(k+mp.sqrt(discr))/(2*ell) if k>0 and discr>=0 else None
        physical=root*scale
        return {'doubleFactor':x,'factor':mp.nstr(root,57),'physicalShift':mp.nstr(physical,57),
                'gpuBracket':[lo,hi],'doubleScoreAtGpuBracket':[d(lo),d(hi)],
                'quadraticProxy':mp.nstr(proxy,57) if proxy is not None else None,
                'shiftMinusProxyTimesProxy':mp.nstr((physical-proxy)*proxy,40) if proxy is not None else None,
                'parameters':p}


def main():
    p=Path(sys.argv[1])
    cases={}
    if p.exists():
        data=json.loads(p.read_text())
        for key,hit in data['extremalCrossings'].items():
            cases[key]=solve(hit['profile'],hit['lowerNu'],hit['upperNu'])
            print(key,cases[key]['factor'],flush=True)
    fixed_profile={'w':0.5,'r':1,'length':1,'b':0.5}
    asymptotic=[]
    for a in (32,128,1024,65536,10**8):
        out=solve({'a':a,**fixed_profile},0.7,1.1)
        asymptotic.append(out)
        print('a=',a,'factor=',out['factor'],'proxy-scaled remainder=',out['shiftMinusProxyTimesProxy'],flush=True)
    Path(sys.argv[2]).write_text(json.dumps({'python':sys.version,'scipy':scipy.__version__,'mpmath':mp.__version__,
       'meshCases':cases,'asymptoticChecks':asymptotic},indent=2)+'\n')


if __name__=='__main__': main()
