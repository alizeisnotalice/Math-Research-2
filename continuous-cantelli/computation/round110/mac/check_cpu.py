"""Independent real-transform, interval Jensen, and late-crossing checks."""
import argparse
from fractions import Fraction
import json
from pathlib import Path

import mpmath as mp
import sympy as sp
from scipy.integrate import quad
import numpy as np
import scipy

mp.mp.dps = 65
mp.iv.dps = 55
LAMS = (0.25, 0.5, 0.75, 1, 1.5, 2, 2.5, 3, 4, 6, 8, 12)

def profile(seed, ctx=mp):
    """Return the mean- and fourth-moment-matched profile in ctx arithmetic."""
    a, b, c, u = map(lambda x: x if hasattr(x, '_mpi_') else ctx.mpf(str(x)), seed)
    A, B, C = a+b+c, a*b+a*c+b*c, a*b*c
    g = [(1+j*u)**(-ctx.mpf('0.5')) for j in range(11)]
    d = (g[5]-A*g[4]+B*g[3]-C*g[2])/(g[4]-A*g[3]+B*g[2]-C*g[1])
    p = [ctx.mpf(0), C*d, -C-B*d, B+A*d, -A-d, ctx.mpf(1)]
    H2 = sum(-j*u*g[j]**3*p[j] for j in range(1,6))
    norm = sum(p[i]*p[j]*g[i+j] for i in range(1,6) for j in range(1,6))
    return (a,b,c,d), u, g, p, 2*H2/norm

def product(s, roots):
    result = s
    for r in roots:
        result *= s-r
    return result

def jensen(data, lam, ctx=mp):
    _,u,g,p,_ = data
    t = ctx.mpf(str(lam))**2/2
    kappa = u/(1+u)
    return sum(p[j]*g[j]*ctx.exp(-t*(j*u/(1+j*u)-kappa)) for j in range(1,6))

def minima(data):
    roots,u,_,p,epsilon = data
    roots = sorted(roots)
    def find(left,right):
        for _ in range(170):
            s=(left+right)/2
            if 1/s+sum(1/(s-r) for r in roots)>0:
                left=s
            else:
                right=s
        return (left+right)/2
    result=[]
    for i in (0,2):
        s = find(roots[i], roots[i+1])
        x = mp.sqrt(-2*mp.log(s)/u)
        value = 1+epsilon*product(s, roots)
        second = epsilon*sum(j*(j-1)*p[j]*s**(j-2) for j in range(2,6))*(u*x*s)**2
        result.append((value,x,second))
    return sorted(result)

def residual64(data, lam):
    roots,u,_,_,epsilon=data
    l, u, eps = map(float,(lam,u,epsilon))
    roots = list(map(float,roots))
    t = l*l/2
    k = u/(1+u)
    h = float(jensen(data,l))
    marks = sorted(set([-18.,0.,18.]+[j*np.sqrt(-2*np.log(s)/u)-l for s in roots for j in (-1,1) if abs(j*np.sqrt(-2*np.log(s)/u)-l)<18]))
    def integrand(y):
        v = -y*y/2+t*k
        w = -t*eps*product(np.exp(-u*(l+y)**2/2),roots)
        if abs(w)<.02:
            rem = (w*w/2)*(1+w*(1/3+w*(1/12+w*(1/60+w/360))))
            return np.exp(v)*rem
        return np.exp(v+w)-np.exp(v)*(1+w)
    n = sum(quad(integrand,a,b,epsabs=2e-12,epsrel=2e-11)[0] for a,b in zip(marks,marks[1:]))
    return -h+n/(np.sqrt(2*np.pi)*t*eps)

def residual_mp(data,lam):
    roots,u,_,_,eps=data
    l=mp.mpf(str(lam));t=l*l/2;k=u/(1+u)
    def fn(y):
        v=-y*y/2+t*k;s=mp.exp(-u*(l+y)**2/2);w=-t*eps*product(s,roots)
        return mp.exp(v)*(mp.expm1(w)-w)
    marks = sorted(set(map(mp.mpf,(-18,-9,0,9,18))) |
                   {j*mp.sqrt(-2*mp.log(s)/u)-l for s in roots for j in (-1,1) if abs(j*mp.sqrt(-2*mp.log(s)/u)-l)<18})
    return -jensen(data,l)+mp.quad(fn,marks)/(mp.sqrt(2*mp.pi)*t*eps)

def interval_profile(seed):
    roots,u,g,p,eps = profile(seed,mp.iv)
    ordered=sorted(roots,key=lambda r:float(r.mid))
    assert ordered[0]>0 and ordered[-1]<1 and eps>0
    assert all(a.b<b.a for a,b in zip(ordered,ordered[1:]))
    def frac_iv(f):
        return mp.iv.mpf(f.numerator)/f.denominator
    def slope(x):
        return 1/x+sum(1/(x-r) for r in ordered)
    certificates = []
    for j in (0,2):
        a,b=ordered[j:j+2]
        lo=Fraction(str(float(a.b)+1e-8)); hi=Fraction(str(float(b.a)-1e-8))
        assert frac_iv(lo)>a.b and frac_iv(hi)<b.a and slope(frac_iv(lo))>0 and slope(frac_iv(hi))<0
        l0,h0,l1,h1=lo,hi,lo,hi
        for _ in range(110):
            m=(l0+h0)/2
            if slope(frac_iv(m))>0: l0=m
            else: h0=m
            m=(l1+h1)/2
            if slope(frac_iv(m))<0: h1=m
            else: l1=m
        s=mp.iv.mpf([frac_iv(l0).a,frac_iv(h1).b])
        value=1+eps*product(s,ordered)
        assert value>0
        certificates.append((s,value))
    return (roots,u,g,p,eps),certificates

def interval_jensen(seed,lam):
    data,certificates=interval_profile(seed)
    ans=jensen(data,lam,mp.iv)
    u,epsilon=data[1],data[4]
    t=mp.iv.mpf(str(lam))**2/2
    log_lower=-t*epsilon*mp.iv.exp(-t*u/(1+u))*ans
    return {'lambda':lam,'lower':float(ans.a),'upper':float(ans.b),'interval':str(ans),
            'negative':bool(ans.b<0),'logJensenLowerInterval':str(log_lower),
            'positiveMinima':[str(v) for _,v in certificates]}

def crossing_certificate():
    seed=('0.19346052','0.6506234','0.71849674')
    bounds=('1.8654002890','1.8654002892')
    q_edges=[]
    for u in bounds:
        data, mins=interval_profile((*seed,u))
        assert mins[0][1]<mins[1][1]
        q_edges.append(mins[0][1]-data[1]/(1+data[1]))
    assert q_edges[0]>0 and q_edges[1]<0
    _,mins=interval_profile((*seed,mp.iv.mpf(list(bounds))))
    assert mins[0][1]<mins[1][1]
    return {'fixedRootTriple':seed,'criticalUBracket':bounds,
            'qEdgeIntervals':list(map(str,q_edges)),'minimaOverWholeBracket':[str(v) for _,v in mins]}

def verify_symbols():
    z,eps = sp.symbols('lambda epsilon')
    x, h = sp.symbols('x h')
    poly = sp.Poly(sp.series(sp.exp(z*x-z*z*(1+eps*h)/2),z,0,7).removeO(),z)
    H2,H4,H6=(sp.hermite_prob(j,x) for j in (2,4,6))
    assert sp.expand(poly.nth(4) - H4/24 + eps*h*H2/4 - eps**2*h*h/8)==0
    formula=H6/720-eps*h*H4/48+eps**2*h*h*H2/16-eps**3*h**3/48
    assert sp.expand(poly.nth(6)-formula)==0
    return {'fourth':str(sp.factor(poly.nth(4))),'sixth':str(sp.factor(poly.nth(6)))}

def run(root,output):
    scan=json.loads(Path(root).read_text())
    report={}
    for name,hit in scan['hits'].items():
        if not hit:
            continue
        seed=hit['roots']
        data=profile(seed)
        roots,u,g,p,eps=data
        assert abs(sum(p[j]*g[j] for j in range(6)))<mp.mpf('1e-60')
        ms=minima(data)
        fmin,xstar,curvature=ms[0]
        q=fmin-u/(1+u)
        replay=[residual64(data,l) for l in LAMS]
        jv=[float(jensen(data,l)) for l in LAMS]
        selected=sorted({int(np.argmax(np.abs(replay))),int(np.argmin(jv)),0})
        precision={str(LAMS[i]):str(residual_mp(data,LAMS[i])) for i in selected}
        delta=max(abs(float(precision[str(LAMS[i])])-replay[i]) for i in selected)
        out=dict(r4=str(roots[-1]),epsilon=str(eps),fmin=str(fmin),q=str(q),xstar=str(xstar),
            curvature=str(curvature), alternateMinimum=str(ms[1][0]),metal=hit,raw64=replay,
            jensen=jv,precision=precision,maxAbsoluteCpuCheckError=delta)
        if min(jv)<0:
            out['intervalJensen']=interval_jensen(seed,LAMS[int(np.argmin(jv))])
        report[name]=out
    Path(output).write_text(json.dumps({'lambdas':LAMS,'mpmath':mp.__version__,'scipy':scipy.__version__,
        'sympy':sp.__version__,'symbols':verify_symbols(),'crossingCertificate':crossing_certificate(),
        'checks':report},indent=2)+'\n')

if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('scan')
    parser.add_argument('out')
    args=parser.parse_args()
    run(args.scan,args.out)
