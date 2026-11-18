"""Independent intersections, adaptive integrals, and scaled large-a checks."""
import json
import math
from pathlib import Path
import sys

import mpmath as mp
import scipy
from scipy.integrate import quad
from scipy.optimize import brentq


def moments(a,w,r,L,b):
    rt=math.sqrt(2)
    normal=math.sqrt(2*math.pi)
    phi=lambda y: math.exp(-y*y/2)/normal
    pv=math.erf(w/(a*rt))
    pb=math.erf((r+L)/rt)-math.erf(r/rt)
    A=b*pb/pv
    m2=1+2*(r*phi(r)-(r+L)*phi(r+L))/pb
    return pv,pb,A,m2


def direct(a,w,r,L,b,nu):
    pv,pb,A,m2=moments(a,w,r,L,b)
    loga=math.log(a)
    s=nu*loga/(a*m2)
    q=(1+s/a)**2
    nv=math.expm1(q*math.log1p(A))
    nb=-math.expm1(q*math.log1p(-b))
    norm=math.sqrt(2*math.pi)
    def signed(y):
        return math.exp(-(y-s)**2/2)+math.exp(-(y+s)**2/2)
    v=quad(signed,0,w/a,epsabs=1e-13,epsrel=1e-13)[0]/norm
    p=quad(signed,r,r+L,epsabs=1e-13,epsrel=1e-13)[0]/norm
    return (nv*v-nb*p)/(b*pb)*a*a*m2/(loga*loga)


def high_mp(profile,nu):
    with mp.workdps(max(65, int(3*math.log10(float(profile["a"])))+30)):
        a,w,r,L,b=map(lambda k: mp.mpf(str(profile[k])),('a','w','r','length','b'))
        nu=mp.mpf(str(nu))
        loga=mp.log(a)
        rt=mp.sqrt(2)
        norm=mp.sqrt(2*mp.pi)
        phi=lambda y: mp.exp(-y*y/2)/norm
        pv=mp.erf(w/(a*rt))
        pb=mp.erf((r+L)/rt)-mp.erf(r/rt)
        A=b*pb/pv
        m2=1+2*(r*phi(r)-(r+L)*phi(r+L))/pb
        s=nu*loga/(a*m2)
        q=(1+s/a)**2
        cdf=lambda x: mp.erfc(-x/rt)/2
        pvs=cdf(w/a-s)-cdf(-w/a-s)
        pbs=cdf(r+L-s)-cdf(r-s)+cdf(-r-s)-cdf(-r-L-s)
        return (mp.expm1(q*mp.log1p(A))*pvs+mp.expm1(q*mp.log1p(-b))*pbs)/(b*pb)*a*a*m2/(loga*loga)


def main():
    payload=json.loads(Path(sys.argv[1]).read_text())
    checked={}
    for key,hit in payload['extremalCrossings'].items():
        p=hit['profile']
        args=[p[x] for x in ('a','w','r','length','b')]
        lo,hi=hit['lowerNu'],hit['upperNu']
        f=lambda nu: direct(*args,nu)
        assert f(lo)>0 and f(hi)<0
        root=brentq(f,lo,hi,xtol=1e-12)
        value_mp=high_mp(p,root)
        assert abs(float(value_mp)) < 1e-6
        pv,pb,A,m2=moments(*args)
        D=(1+1/A)*math.log1p(A)+(1-p['b'])/p['b']*math.log1p(-p['b'])
        root_physical=root*math.log(p['a'])/(p['a']*m2)
        checked[key]={'parameters':hit,'scipyBracket':[f(lo),f(hi)],'rootNu':root,'physicalShift':root_physical,
                      'rootWithDerivativeScaling':root*math.log(p['a'])/D,'mpmathResidualAtRoot':mp.nstr(value_mp,45)}
        print(f'{key}: nu={root:.12g} shift={root_physical:.10g} derivative-scaled={root*math.log(p["a"])/D:.10g}')
    fixed={'w':0.5,'r':1,'length':1,'b':0.5}
    asympt=[]
    for a in (1e3,1e6,1e12,1e24,1e48):
        for nu in (-2,1,3,4,6):
            val=high_mp({'a':a,**fixed},nu)
            asympt.append({'a':a,'nu':nu,'value':mp.nstr(val,32),'limitingValue':2*nu-nu**2/2})
    Path(sys.argv[2]).write_text(json.dumps({'python':sys.version,'scipy':scipy.__version__,'mpmath':mp.__version__,
       'crossings':checked,'fixedAsymptoticProfile':fixed,'asymptoticChecks':asympt},indent=2)+'\n')


if __name__=='__main__': main()
