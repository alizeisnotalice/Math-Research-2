"""Check global Hermite-to-call asymptotic under added Gaussian variance."""
import json
from pathlib import Path
import math
import sys
from check_cpu import funs,moments,high


def main():
    base={'k':3.,'l':4.,'coefficient':0.,'mean':1.}
    delta6=8*moments(base)['standardizedHermite6']
    assert 5<delta6<6
    rows=[]
    for m in (1.,4.,16.,64.,256.):
        f,d=funs({**base,'mean':m})
        M=1+m
        for y in (0.,1.4,2.6):
            a=f(y);b=d(y)
            assert abs(a-b)<1e-11
            scaled=a*720*M**3/delta6
            exact=(y**4-6*y*y+3)*math.exp(-y*y/2)/math.sqrt(2*math.pi)
            rows.append({'meanNoiseVariance':m,'standardizedStrike':y,'callResidual':a,
                         'rescaledResidual':scaled,'targetHermiteShape':exact,'rescaledError':scaled-exact})
    for y in (0.,1.4):
        mpcheck=float(high({**base,'mean':256.},y))
        row=next(x for x in rows if x['meanNoiseVariance']==256. and x['standardizedStrike']==y)
        row['mp60']=mpcheck
        assert abs(row['callResidual']-mpcheck)<1e-14
    Path(sys.argv[1]).write_text(json.dumps({'base':base,'rawSixthMomentExcess':delta6,'rows':rows},indent=2)+'\n')
    print('raw m6 excess',delta6,'last layer',rows[-3:],flush=True)


if __name__=='__main__':main()
