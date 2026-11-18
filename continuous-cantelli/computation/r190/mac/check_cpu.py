"""Independent real-space difference-law audit of weighted Cantor Fourier windows."""
import json
import math
from pathlib import Path
import platform
import sys

import mpmath as mp
import numpy as np
from scipy.fft import irfft, next_fast_len, rfft
import scipy

mp.mp.dps=64
D=math.log(2)/math.log(3)
CRITICAL=1-D

def intervals(n):
    j=np.array([0],dtype=np.int64)
    for k in range(n):
        j=np.r_[j,j+3**k]
    return np.sort(j)

def difference_rule(n):
    r=np.array([0.])
    w=np.array([1.])
    for j in range(1,n+1):
        delta=2*3.**(-j)
        r=np.r_[r-delta,r,r+delta]
        w=np.r_[.25*w,.50*w,.25*w]
    return r,w

def spatial_window(case, extra, n_tail=7):
    m=case['m']; W=case['W']; offset=case['offset']
    N=min(13,m+extra)
    n=3**N; inv=1/n
    i=intervals(N)
    centers=(2*i-(n-1)//2)*inv
    g=np.exp(-.5*(2+centers)**2)/math.sqrt(2*math.pi)
    j=np.zeros((n+1)//2)
    j[i]=g/len(i)
    nfft=next_fast_len(2*len(j)-1,real=True)
    ft=rfft(j,nfft)
    corr=irfft(ft.real**2+ft.imag**2,nfft)[:len(j)]
    r,w=difference_rule(n_tail)
    cutoff=min(len(j),math.ceil((9/W+inv)/(2*inv))+1)
    ctotal=total=0.
    for lo in range(0,cutoff,192):
        k=np.arange(lo,min(cutoff,lo+192),dtype=np.int64)
        lag=2*inv*k[:,None]
        u=lag+inv*r
        slope=1-inv*lag*r/2
        base=np.exp(-.5*(W*u)**2)*slope
        base=np.multiply(base,w)
        phase=np.zeros((len(k),1)) if m>=N else 2*math.pi*((2*k)%(3**(N-m)))[:,None]/(3**(N-m))
        phase=phase+offset*lag+(case['xi']+offset)*inv*r
        multiplier=np.where(k==0,1.,2.)*np.maximum(corr[k],0.)
        total+=float(np.dot(multiplier,np.sum(base*np.cos(phase),axis=1)))
        ctotal+=float(np.dot(multiplier,np.sum(base,axis=1)))
    factor=W/math.sqrt(2*math.pi)
    return {'N':N,'tailDepth':n_tail,'I':total*factor,'I0':ctotal*factor,
        'ratio':total/ctotal,'lagCount':cutoff}

def point_power(m,frequency_offset=0):
    N=17;n=3**N
    i=intervals(N)
    j=2*i-(n-1)//2
    x=j/float(n)
    g=np.exp(-.5*(2+x)**2)/math.sqrt(2*math.pi)
    res=3**(N-m)
    phase=(2*math.pi)*(j%res)/res+frequency_offset*x
    bare=np.mean(g*(np.cos(phase)+1j*np.sin(phase)))
    T=2*math.pi/res+frequency_offset/n
    tail=math.prod(math.cos(T*3.**(-k)) for k in range(1,25))
    center=abs(bare*tail)**2
    return center,float(np.mean(g))**2

def main():
    path,out=map(Path,sys.argv[1:3])
    j=json.loads(path.read_text())
    subset=[a for a in j['rows'] if a['m'] in (2,7,12) and a['step']==.5
                and abs(a['c']-CRITICAL)<1e-12 and a['A']==.5]
    out_rows=[]
    for a in subset:
        finite=spatial_window(a,4,7)
        refined=spatial_window(a,5,8)
        center,unshift=point_power(a['m'])
        row={'gpu':a,'realSpaceCoarse':finite,'realSpaceRefined':refined,'centerPowerByAddresses':center,
             'zeroPowerByAddresses':unshift,'absRelativeWindowError':abs(refined['I']-a['I'])/refined['I'],
             'absRelativeWindowDepthDelta':abs(refined['I']-finite['I'])/refined['I'],
             'absRelativeCenterWindowError':abs(refined['I0']-a['I0'])/refined['I0'],
             'absCenterPowerError':abs(a['centerPower']-center)}
        assert row['absRelativeWindowError']<3e-4,(a,row)
        assert row['absRelativeCenterWindowError']<3e-4,(a,row)
        out_rows.append(row)
    limits=mp.fprod(mp.cos(2*mp.pi/3**n) for n in range(1,96))
    by=[]
    for a in j['rows']:
        if a['step']>=.5:continue
        pair=next(b for b in j['rows'] if b['m']==a['m'] and b['beta']==a['beta'] and b['c']==a['c']
            and b['A']==a['A'] and b['j']==a['j'] and b['step']==.5)
        by.append({'beta':a['beta'],'m':a['m'],'IHalfStep':a['I'],'IFullStep':pair['I'],
                   'relativeDifference':abs(a['I']-pair['I'])/a['I']})
    out.write_text(json.dumps({'python':platform.python_version(),'mpmath':mp.__version__,
        'numpy':np.__version__,'scipy':scipy.__version__,'persistentProduct':mp.nstr(limits,52),
        'limitingPowerRatio':mp.nstr(limits**2,52),'checked':out_rows,'dyadicStepCheck':by},indent=2)+'\n')
if __name__=='__main__':main()
