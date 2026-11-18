from pathlib import Path
import json

import numpy as np
from scipy.integrate import simpson

import bound_primitives as b

ROOT=Path(__file__).parent
H=float(b.H)
WEIGHTS=[(kind,float(t),float(y),w/b.WEIGHT_DEN) for kind,t,y,w in b.TESTS]


def values(x,f):
    out=f.copy()
    for kind,t,y,w in WEIGHTS:
        if kind=='t':
            out+=w*np.exp(t*x-.5*t*t*f)
        elif kind=='h':
            var=H+t-f
            out+=w*np.exp(-.5*(x-y)**2/var)/np.sqrt(var)
    return out


def strip(x,a):
    amin=np.where(np.abs(x)<=.75,np.maximum(0,a-np.abs(x-.5)),0)
    phases=np.linspace(0,1,257)
    grid=amin[:,None]+(H-amin[:,None])*phases
    xcol=x[:,None]
    q=values(xcol,grid)
    index=np.argmin(q,axis=1)
    q_old=q[np.arange(len(x)),index]
    low=grid[np.arange(len(x)),np.maximum(index-1,0)]
    high=grid[np.arange(len(x)),np.minimum(index+1,len(phases)-1)]
    for _ in range(31):
        m=(2*low+high)/3
        n=(low+2*high)/3
        q0=values(x,m)
        q1=values(x,n)
        low=np.where(q0>q1,m,low)
        high=np.where(q0<=q1,n,high)
    return np.minimum(q_old,values(x,(low+high)/2))*np.exp(-.5*x*x)/np.sqrt(2*np.pi)


def margin(a,n):
    # Integrate the smooth pieces separately at the known modulus boundaries.
    z=0.
    for left,right in [(-8,-.75),(-.75,.5),(.5,.75),(.75,8)]:
        x=np.linspace(left,right,n+1)
        if left==.75:
            x[0]=np.nextafter(x[0],np.inf)
        if right==-.75:
            x[-1]=np.nextafter(x[-1],-np.inf)
        z+=simpson(np.concatenate([strip(y,a) for y in np.array_split(x,16)]),x=x)
    target=sum(w*(np.exp(-y*y/(2*(H+t)))/np.sqrt(H+t) if kind=='h' else 1)
               for kind,t,y,w in WEIGHTS)
    return z-target


def main():
    anchors=[437/256,6999/4096,1.71,1.711]
    rows=[dict(anchor=a,coarse=margin(a,1024),fine=margin(a,4096)) for a in anchors]
    (ROOT/'calibration.json').write_text(json.dumps(rows,indent=2)+'\n')
    print(json.dumps(rows))


if __name__=='__main__':
    main()
