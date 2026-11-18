"""Verify contracted-center positive interval ratios in two coordinate systems."""
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
    return (c['center'],c['width'],c['compression'],c['drift'],hit['offset'],hit['halfWidth'])


def check(hit):
    a,w,k,eta,offset,d=unpack(hit)
    beta=2*math.log(a)+eta
    a2=a*a
    def dual(z):
        xs=1-z/a2
        f=xs/(1+(beta-k*z)/a2)
        sigma_inv=math.sqrt(f)
        mu=(beta-k*z)/a
        lo=(offset-d-mu)*sigma_inv
        hi=(offset+d-mu)*sigma_inv
        probability=ndtr(-lo)-ndtr(-hi) if lo>=0 else ndtr(hi)-ndtr(lo)
        exponent=xs*(beta+(1-k)*z)/2
        return math.exp(exponent)*probability/sigma_inv/(2*d*a)

    def primal(u):
        def integrate(z):
            gap=(beta+(1-k)*z)/(a2+beta-k*z)
            action=(-z/a-u)**2*(-0.5)+gap*(a+u)**2/2
            return math.exp(action)/math.sqrt(2*math.pi)
        return quad(integrate,0,w,epsabs=1e-12,epsrel=1e-12)[0]/(2*d*a)

    v_dual=quad(dual,0,w,epsabs=1e-12,epsrel=1e-12)[0]
    v_primal=quad(primal,offset-d,offset+d,epsabs=1e-12,epsrel=1e-12)[0]
    assert abs(v_dual-v_primal)<1e-10*max(1,v_dual)
    with mp.workdps(55):
        a,w,k,eta,offset,d=(mp.mpf(str(v)) for v in unpack(hit))
        beta=2*mp.log(a)+eta
        def mp_dual(z):
            xs=1-z/a**2
            f=xs/(1+(beta-k*z)/a**2)
            sigma_inv=mp.sqrt(f)
            mu=(beta-k*z)/a
            lo=(offset-d-mu)*sigma_inv/mp.sqrt(2)
            hi=(offset+d-mu)*sigma_inv/mp.sqrt(2)
            probability=(mp.erfc(lo)-mp.erfc(hi))/2
            return mp.exp(xs*(beta+(1-k)*z)/2)*probability/sigma_inv/(2*d*a)
        mp_value=mp.quad(mp_dual,[0,w/2,w])
        assert abs(float(mp_value)-v_dual)<1e-10*max(1,v_dual)
    return {'parameters':hit,'scipyDual':v_dual,'scipyPrimal':v_primal,'mpmathDual':mp.nstr(mp_value,48),
            'metalRelativeError':(hit['ratio']-v_dual)/v_dual}


def main():
    payload=json.loads(Path(sys.argv[1]).read_text())
    groups={}
    for group in ('bestByRegime','minimumDriftByRegime'):
        selected={}
        for name,hit in payload[group].items():
            result=check(hit)
            selected[name]=result
            print(f'{group}/{name}: ratio {result["scipyDual"]:.12g}, eta {hit["core"]["drift"]}',flush=True)
        groups[group]=selected
    # Limiting drift thresholds for the centered J with fixed observation width.
    prob=ndtr(1/16)-ndtr(-1/16)
    thresholds=[]
    for w in (0.25,0.5,1,2,4):
        for k in (0,0.5,1):
            b=(1-k)/2
            scale=w if b==0 else math.expm1(b*w)/b
            thresholds.append({'w':w,'kappa':k,'halfWidth':1/16,'etaThreshold':2*math.log((1/8)/(prob*scale))})
    Path(sys.argv[2]).write_text(json.dumps({'python':sys.version,'scipy':scipy.__version__,'mpmath':mp.__version__,
        'checks':groups,'limitingThresholds':thresholds},indent=2)+'\n')


if __name__=='__main__': main()
