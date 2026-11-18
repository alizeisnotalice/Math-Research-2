"""Cross-check Metal crossing brackets with adaptive physical-coordinate integrals."""
import json
import math
from pathlib import Path
import sys

import mpmath as mp
import scipy
from scipy.integrate import quad
from scipy.optimize import brentq


def smooth(s):
    if s <= 0: return 0.0
    if s >= 1: return 1.0
    f = math.exp(-1/s)
    return f/(f+math.exp(-1/(1-s)))


def physical(hit, kappa):
    p = hit['profile']
    l,h,m,d,w = map(p.__getitem__, ('lambda','h','m','downWidth','plateauWidth'))
    c = hit['curvature']
    normal = math.sqrt(2*math.pi)
    rout = 1/math.sqrt(c)
    rin = rout/2
    s = math.sqrt(1-m)
    anchor = l*(1+s) + kappa*math.log(l)/(l*s)
    d, w = d/(l*s), w/(l*s)
    left, right = anchor-d, anchor+w+d
    assert left > l > rout

    def corrected(x,v):
        e = l*x - x*x/2 - l*l*v/2
        return math.exp(e)*(-math.expm1(-l*l*(h-v)/2))/normal

    def origin(x):
        cutoff = smooth(2*math.sqrt(c)*abs(x)-1)
        return corrected(x,h*((c*x*x)*(1-cutoff)+cutoff))

    def island(x):
        if x < anchor:
            taper = smooth((anchor-x)/d)
        elif x > anchor+w:
            taper = smooth((x-anchor-w)/d)
        else:
            taper = 0.0
        return corrected(x,m+(h-m)*taper)

    def integrate(f, edges):
        return sum(quad(f,a,b,epsabs=1e-12,epsrel=1e-12,limit=250)[0] for a,b in zip(edges[:-1],edges[1:]))

    compact = integrate(origin,[-rout,-rin,0,rin,rout])
    band = integrate(island,[left,anchor,anchor+w,right])
    background = math.exp(-(h-1)*l*l/2)
    return compact, band, background


def high_precision(hit, kappa, digits=55):
    with mp.workdps(digits):
        p=hit['profile']
        l,h,m,d,w=(mp.mpf(str(p[key])) for key in ('lambda','h','m','downWidth','plateauWidth'))
        kappa=mp.mpf(str(kappa))
        c=mp.mpf(str(hit['curvature']))
        ri=1/(2*mp.sqrt(c))
        ro=2*ri
        norm=mp.sqrt(2*mp.pi)
        s=mp.sqrt(1-m)
        cur=kappa*mp.log(l)
        penalty=l*l*(h-m)/2

        def step(u):
            if u <= 0: return mp.mpf(0)
            if u >= 1: return mp.mpf(1)
            return 1/(1+mp.exp(1/u-1/(1-u)))

        def origin(x):
            v=step(2*mp.sqrt(c)*abs(x)-1)
            f=h*((1-v)*c*x*x+v)
            return mp.exp(l*x-x*x/2-l*l*f/2)*(-mp.expm1(-l*l*(h-f)/2))/norm

        def island(t):
            bend=step(-t/d) if t<0 else step((t-w)/d) if t>w else mp.mpf(0)
            exponent= -cur-t-(cur+t)**2/(2*(l*s)**2)-penalty*bend
            return mp.exp(exponent)*(-mp.expm1(-penalty*(1-bend))) / (l*s*norm)

        return mp.quad(origin,[-ro,-ri,0,ri,ro]) + mp.quad(island,[-d,-d/2,0,w,w+d/2,w+d])+mp.exp(-(h-1)*l*l/2)


def main():
    grid=json.loads(Path(sys.argv[1]).read_text())
    values={}
    for regime,hit in grid['bestByRegime'].items():
        low,high=hit['lowerKappa'],hit['upperKappa']
        f=lambda k: sum(physical(hit,k))
        low_log,high_log=math.log(f(low)),math.log(f(high))
        root=brentq(lambda k: f(k)-1,low-0.02,high+0.02,xtol=1e-12)
        oracle=high_precision(hit,root)
        assert abs(float(oracle)-1)<1e-10
        p=hit['profile']
        l,m=p['lambda'],p['m']
        compact,band,background=physical(hit,root)
        z=1+math.sqrt(1-m)+root*math.log(l)/(l*l*math.sqrt(1-m))
        values[regime]={'parameters':hit,'adaptiveLowerLogF':low_log,'adaptiveUpperLogF':high_log,
                        'adaptiveCrossingKappa':root,'mpmathFAtCrossing':mp.nstr(oracle,40),
                        'crossingScaledPlateauStart':z,'compact':compact,'tail':band,'background':background}
        print(f'{regime}: kappa={root:.12g}; compact={compact:.9g}; tail={band:.9g}; z={z:.12g}')
    Path(sys.argv[2]).write_text(json.dumps({'scipy':scipy.__version__,'mpmath':mp.__version__,
       'python':sys.version,'crossings':values},indent=2)+'\n')


if __name__=='__main__': main()
