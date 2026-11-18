"""Full spectral-tail sums and exact intervals for bounded Hermite profiles."""
import argparse,json,math,os,time,platform
from fractions import Fraction as F
from pathlib import Path
import numpy as np
import scipy
from scipy.stats import qmc


def direct_float(delta):
    K=np.arange(math.ceil(4/delta)+3,dtype=np.float64)
    t=np.ones(len(K))
    b=np.zeros(len(K))
    for s in range(1,29):
        t*=delta/2/(s+1)*np.sqrt(np.maximum(K-2*s+2,0)*np.maximum(K-2*s+1,0))
        b+=t
    assert b[-1]>1 and np.all(b[1:]>=b[:-1])
    caps=np.minimum(1,b*b)
    tails=1-caps
    energy=np.sum(tails)
    # Independent summation by parts on the extremal spectral law.
    pmf=np.diff(np.concatenate(([0.],caps)))
    direct=np.sum(K*pmf)
    return float(energy),int(np.searchsorted(caps,1.)),float(abs(energy-direct))


def bracket(delta,denom=10**30,terms=21):
    energy_lo=energy_hi=0
    p=delta/2
    coeff=[p**s/F(math.factorial(s+1)) for s in range(terms+1)]
    for k in range(math.ceil(F(4)/delta)+10):
        falling=1
        lo=F(0); hi=F(0)
        r=min(k//2,terms)
        for s in range(1,r+1):
            falling*=(k-2*s+2)*(k-2*s+1)
            sq=math.isqrt(falling*denom**2)
            lo+=coeff[s]*F(sq,denom)
            hi+=coeff[s]*F(sq+1,denom)
        if r<k//2:
            a=delta*k/2
            tail=a**(r+1)/math.factorial(r+2)
            assert a<r+3
            hi+=tail/(1-a/(r+3))
        if lo>=1:
            return F(energy_lo,denom),F(energy_hi,denom),k
        l=max(F(0),1-hi*hi)*denom
        u=max(F(0),1-lo*lo)*denom
        energy_lo+=math.floor(l)
        energy_hi+=math.ceil(u)
    raise AssertionError("cutoff not reached")


def main():
    p=argparse.ArgumentParser();p.add_argument('shard',type=int);p.add_argument('--exponent',type=int,default=4)
    p.add_argument('--seconds',type=float,default=100);p.add_argument('--out',type=Path,required=True);a=p.parse_args()
    start=time.monotonic();cpu=time.process_time();rows=[]
    exact=[]
    d0=F(1,20+5*a.shard)
    d1=F(3,83+13*a.shard)
    for d in (d0,d1):
        x,y,k=bracket(d)
        fx,_,_=direct_float(float(d))
        assert float(x)-1e-10<=fx<=float(y)+1e-10
        exact.append(dict(delta=str(d),lower=str(x),upper=str(y),cutoff=k,computed=fx))
    while len(rows)<4096 and (not rows or time.monotonic()-start<a.seconds):
        m=qmc.Sobol(1,scramble=True,rng=np.random.default_rng(791623+a.shard+15*len(rows))).random_base2(a.exponent).ravel()
        ds=np.exp2(-.5-16*(a.shard+m)/15)
        vals=[direct_float(float(d)) for d in ds]
        rows.append(dict(count=len(ds),delta_min=float(min(ds)),delta_max=float(max(ds)),
                         max_cutoff=max(v[1] for v in vals),poly_prec=max(v[2] for v in vals),
                         scaled_energy=[min(d*v[0] for d,v in zip(ds,vals)),max(d*v[0] for d,v in zip(ds,vals))]))
    a.out.write_text(json.dumps(dict(shard=a.shard,pid=os.getpid(),elapsed=time.monotonic()-start,cpu=time.process_time()-cpu,
                   software=[platform.python_version(),np.__version__,scipy.__version__],exact=exact,results=rows),indent=2)+'\n')
if __name__=='__main__':main()
