"""Calibrate tail-well dual witness thresholds in the large-center limit."""
import json
import math
from pathlib import Path

import scipy
from scipy.integrate import quad
from scipy.optimize import brentq
from scipy.special import ndtr

from check_cpu import dual


def limit_ratio(w, p, h, d, offset=0):
    def integrand(y):
        if d == 0:
            average = math.exp(-(offset-y)**2/2)/math.sqrt(2*math.pi)
        else:
            average = (ndtr(offset+d-y)-ndtr(offset-d-y))/(2*d)
        return math.exp(h*(1-(y/w)**2)**p/2)*average
    return quad(integrand,-w,w,epsabs=1e-13,epsrel=1e-13)[0]


def main():
    results=[]
    for w in (0.25,0.5,1,2,6):
        for p in (1,2,4):
            zero=brentq(lambda h: limit_ratio(w,p,h,0)-1,0,30,xtol=1e-13)
            finite=brentq(lambda h: limit_ratio(w,p,h,1/32)-1,0,30,xtol=1e-13)
            results.append({'width':w,'power':p,'pointwiseDepth':zero,'intervalDepth':finite})
    scaled=[]
    for offset in (0,0.5):
        limiting=limit_ratio(1,2,4,1/8,offset=offset)
        for a in (16,32,64,128,256,512):
            hit={'core':{'center':a,'width':1,'power':2,'depth':4},'offset':offset,'halfWidth':1/8}
            value=dual(hit)
            scaled.append({'center':a,'offset':offset,'ratio':value,'limitingRatio':limiting,
                           'errorTimesA':(value-limiting)*a,'errorTimesASquared':(value-limiting)*a*a})
    Path('limit_check.json').write_text(json.dumps({'scipy':scipy.__version__,'thresholds':results,'finiteCenterScaling':scaled},indent=2)+'\n')
    for w in results:
        if w['power']==2:print(w)
    print(scaled[-7],scaled[-1])


if __name__=='__main__':main()
