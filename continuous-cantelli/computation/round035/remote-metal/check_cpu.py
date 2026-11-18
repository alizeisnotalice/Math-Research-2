"""Check GPU-selected ratios by the Gaussian dual and original positive kernel."""
import json
import math
from pathlib import Path
import sys

import mpmath as mp
import scipy
from scipy.integrate import quad
from scipy.special import ndtr


def unpack(hit):
    c=hit['core']
    return (c['center'],c['width'],c['power'],c['depth'],hit['offset'],hit['halfWidth'])


def dual(hit):
    a,w,p,h,offset,d=unpack(hit)

    def integrand(y):
        gap=h*(max(0,1-(y/w)**2))**p/a**2
        f=1-gap
        mean_offset=(y+a*gap)/f
        lo=(offset-d-mean_offset)*math.sqrt(f)
        hi=(offset+d-mean_offset)*math.sqrt(f)
        probability=ndtr(-lo)-ndtr(-hi) if lo>=0 else ndtr(hi)-ndtr(lo)
        return math.exp((a+y)**2*gap/(2*f))/math.sqrt(f)*probability

    return quad(integrand,-w,w,epsabs=1e-12,epsrel=1e-12,points=[0])[0]/(2*d)


def direct(hit):
    a,w,p,h,offset,d=unpack(hit)
    norm=math.sqrt(2*math.pi)
    def outer(u):
        def inner(y):
            shape=max(0,1-(y/w)**2)**p
            action=h*shape*((a+u)/a)**2/2 - (y-u)**2/2
            return math.exp(action)/norm
        return quad(inner,-w,w,epsabs=1e-12,epsrel=1e-12,points=[0])[0]
    return quad(outer,offset-d,offset+d,epsabs=1e-12,epsrel=1e-12)[0]/(2*d)


def high_precision(hit):
    with mp.workdps(55):
        a,w,p,h,offset,d=map(lambda v:mp.mpf(str(v)),unpack(hit))
        norm=mp.sqrt(2)
        def integrand(y):
            gap=h*(max(0,1-(y/w)**2))**p/a**2
            f=1-gap
            mean_offset=(y+a*gap)/f
            lo=(offset-d-mean_offset)*mp.sqrt(f)
            hi=(offset+d-mean_offset)*mp.sqrt(f)
            mass=(mp.erfc(lo/norm)-mp.erfc(hi/norm))/2
            return mp.exp((a+y)**2*gap/(2*f))*mass/mp.sqrt(f)
        return mp.quad(integrand,[-w,-w/2,0,w/2,w])/(2*d)


def main():
    src=json.loads(Path(sys.argv[1]).read_text())
    stages={}
    for group in ('bestByRegime','minimumDepthByRegime'):
        results={}
        for name,hit in src[group].items():
            value=dual(hit)
            primal=direct(hit)
            assert abs(primal-value) < 1e-10*max(1,value)
            mpcheck=high_precision(hit)
            assert abs(float(mpcheck)-value) < 1e-10*max(1,value)
            results[name]={'parameters':hit,'scipyDualRatio':value,'scipyPrimalRatio':primal,
                           'mpmathDual55Digits':mp.nstr(mpcheck,48),'metalRelativeError':(hit['ratio']-value)/value}
            print(f'{group} / {name}: metal {hit["ratio"]:.9g}; double {value:.11g}; depth {hit["core"]["depth"]}',flush=True)
        stages[group]=results
    Path(sys.argv[2]).write_text(json.dumps({'python':sys.version,'scipy':scipy.__version__,'mpmath':mp.__version__,'checks':stages},indent=2)+'\n')


if __name__=='__main__':
    main()
