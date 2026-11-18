#!/usr/bin/env python3
"""Independent confluent-hypergeometric evaluation of the kernel and shell sums."""
import json
import sys
from pathlib import Path

import numpy as np
import scipy
from scipy.optimize import brentq
from scipy.special import beta,hyp1f1

N=np.array([0,1,2,4,8])


def evaluate(z,k):
    a=k+.5
    w=-.5*z*z
    return beta(a,1-k)/(2*k)*(hyp1f1(a,1.5,w)-(a/1.5)*z*z*hyp1f1(a+1,2.5,w))


def endpoints(m):
    r=m**-2
    a=np.r_[0.,np.linspace(.9,1-r,m-1)]
    return r,a,a.mean()/(1-r)


def depths(m,depth):
    r,a,mean=endpoints(m)
    starts=np.zeros(1)
    for i in range(depth):
        starts=(starts[:,None]+r**i*a[None,:]).ravel()
    return starts+r**depth*mean


def scaled(m,D,k,depth,refine=False):
    r,a,mean=endpoints(m)
    base=float(evaluate(D*depths(m,depth),k).mean())
    suffix=depths(m,1) if refine else np.array([mean])
    y=a[1:,None]+r*suffix[None,:]
    sums=[base]
    for j in range(1,9):
        outer=(m**(j-1))*np.sum(evaluate(D*r**(-j)*y,k))/len(suffix)
        sums.append(sums[-1]+outer)
    return np.asarray(sums)[N]


def main():
    data=json.loads(Path(sys.argv[1]).read_text())
    rows=data["rows"]
    assert np.isfinite([r["scales"] for r in rows]).all()
    selected=[]
    values=np.zeros(5)
    spectral=[]
    for r in rows:
        m,D,k,d=(r[x] for x in ("m","D","kappa","depth"))
        if not (d==2 or D in (5.,10.,15.)):
            continue
        ref=scaled(m,D,k,d)
        assert np.isfinite(ref).all()
        diff=np.abs(ref-np.array(r["scales"]))
        values=np.maximum(values,diff)
        record=dict(m=m,D=D,kappa=k,depth=d,gpu=r["scales"],cpu=ref.tolist(),maximum_error=float(diff.max()))
        if d==2 and m<=24 and D in (5.,10.,15.):
            record["cpu_third_depth"]=scaled(m,D,k,3,True).tolist()
        if d==1 and m in (32,64,128,256,512):
            record["cpu_second_depth"]=scaled(m,D,k,2,True).tolist()
        selected.append(record)
    for m in sorted({r["m"] for r in rows}):
        for k in (.25,.5):
            region=[r for r in rows if r["depth"]==1 and r["m"]==m and r["kappa"]==k]
            array=np.array([r["scales"][-1] for r in region])
            rel=np.flatnonzero(array[1:]*array[:-1]<0)
            transitions=[]
            for j in rel:
                pair=(region[j]["D"],region[j+1]["D"])
                D0=max(5,pair[0]-0.05)
                D1=min(15,pair[1]+0.05)
                cross=brentq(lambda D:scaled(m,D,k,2)[-1],D0,D1,xtol=5e-13)
                transitions.append(cross)
            spectral.append(dict(m=m,kappa=k,threshold=transitions,negatives=int(np.sum(array<0)),
                                 positives=int(np.sum(array>0)),n8_min=float(array.min()),n8_max=float(array.max()),
                                 n0_n8_max_abs=float(max(abs(r["scales"][0]-r["scales"][-1]) for r in region))))
    output=dict(numpy=np.__version__,scipy=scipy.__version__,max_abs_errors=values.tolist(),samples=selected,regions=spectral)
    Path(sys.argv[2]).write_text(json.dumps(output,indent=2)+"\n")
    print("max errors",values)
    for r in spectral:
        print(r)


if __name__=="__main__":
    main()
