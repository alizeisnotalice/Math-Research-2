"""Independent Cantor gap quadrature and explicit finite-level checks."""
import json
from pathlib import Path
import sys
import mpmath as mp
import numpy as np
import scipy
from scipy.special import roots_legendre


def digits(B,mask):
    return [d for d in range(B) if (mask>>d)&1]


def components(row, hit):
    B=row['base']
    D=digits(B,hit['mask'])
    return B,D,hit['depth'],row['exponent'],hit['height']


def basic(D):
    return [(left+1,right) for left,right in zip(D,D[1:]) if right>left+1]


def nodes(count, midpoint=False):
    if midpoint:
        z=(np.arange(count)+0.5)/(2*count)
        return z,np.full(count,1/count)
    t,w=roots_legendre(count)
    t=(1+t)/2
    return t**4/2,2*t**3*w


def gap_eval(n,power,height,centers,widths,rule):
    z,weights=rule
    x=widths[:,None]*z[None,:]
    f=-np.expm1(-0.5*n*n*height*x**power)
    loss=widths*(f@weights)
    shift=widths[:,None]*(z[None,:]-0.5)
    gain=widths*((f*np.cos(2*np.pi*n*shift))@weights)
    return gain*np.exp(2j*np.pi*n*centers),loss


def factored(params,rule):
    B,D,N,power,height=params
    pairs=np.array(basic(D),dtype=float)/B
    output=[];mass=[]
    for n in range(1,7):
        amp=1+0j;total=0j;z=0.; scale=1.; live=1.
        for _ in range(N):
            width=scale*(pairs[:,1]-pairs[:,0]);center=scale*(pairs[:,0]+pairs[:,1])/2
            g,l=gap_eval(n,power,height,center,width,rule)
            total+=amp*live/scale*g.sum();z+=live/scale*l.sum()
            amp*=np.mean(np.exp(2j*np.pi*n*np.array(D)*scale/B))
            scale/=B;live*=len(D)/B
        output.append(-total);mass.append(z)
    return np.array(output),np.array(mass)


def enumerate_check(params,rule):
    B,D,N,power,height=params
    full=np.array([0],dtype=np.int64)
    centers=[];widths=[]
    for k in range(N):
        denom=B**(k+1)
        for left,right in basic(D):
            centers.append((B*full+0.5*(left+right))/denom)
            widths.append(np.full(full.shape,(right-left)/denom))
        if k+1<N:full=(B*full[:,None]+np.array(D)[None,:]).ravel()
    center=np.concatenate(centers);width=np.concatenate(widths)
    values=[];mass=[]
    for n in range(1,7):
        g,l=gap_eval(n,power,height,center,width,rule)
        values.append(-g.sum());mass.append(l.sum())
    return np.array(values),np.array(mass)


def score(out,mass):
    v=np.abs(out)
    f=v/mass
    return [float(f.max()),float(f[:3].max()),float(v.max()),float(v[0])]


def main():
    source=json.loads(Path(sys.argv[1]).read_text())
    rules=nodes(128),nodes(256),nodes(source['panels'],midpoint=True)
    verified=[];explicit=0
    for row in source['rows']:
        for family in ('best3','best6'):
            hit=row[family];p=components(row,hit);B,D,N,alpha,h=p
            a,b=factored(p,rules[0]);x,y=factored(p,rules[1]);u,v=factored(p,rules[2])
            independence=float(max(np.max(np.abs(a-x)),np.max(np.abs(b-y))))
            assert independence<1e-10
            current=score(x,y);mid=score(u,v)
            rec={'base':B,'digits':D,'depth':N,'exponent':alpha,'height':h,'densityIndex':row['densityIndex'],'selector':family,
                 'metal':hit['scores'],'cpuScores':current,'doubleAtMetalMidpoint':mid,'legendreMaxDifference':independence,
                 'metalScoreRoundoffVsFloat64Midpoint':float(np.max(np.abs(np.array(mid)-hit['scores']))),
                 'midpointVsLegendreScores':float(np.max(np.abs(np.array(current)-mid)))}
            if explicit<24 and (len(D)**(N-1))*len(basic(D))<=5000:
                j,k=enumerate_check(p,rules[0])
                err=float(max(np.max(np.abs(j-x)),np.max(np.abs(k-y))))
                assert err<1e-10
                rec['explicitIntervalEnumerationDifference']=err
                explicit+=1
            tail=np.array([(len(D)/B)**N*(n*n*h/2)*(B**(-N)/2)**alpha for n in range(1,7)])
            rec['tailAbsoluteErrorBound']=tail.tolist()
            rec['largestRawMagnitudeMinusTailBound']=float(np.max(np.abs(x)-tail))
            rec['cpuFourier']=[[float(z.real),float(z.imag)] for z in x]
            verified.append(rec)
    report={'python':sys.version,'scipy':scipy.__version__,'numpy':np.__version__,'explicitChecks':explicit,'verified':verified}
    Path(sys.argv[2]).write_text(json.dumps(report,indent=2)+'\n')
    print('checks',len(verified),'explicit',explicit,'largest midpoint discrepancy',max(r['midpointVsLegendreScores'] for r in verified))
    print('largest Float32 roundoff',max(r['metalScoreRoundoffVsFloat64Midpoint'] for r in verified),flush=True)


if __name__=='__main__':main()
