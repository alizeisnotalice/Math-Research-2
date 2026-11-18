"""Gaussian-capped fixed-mass and fixed-mean maximum second moment."""
import argparse,json,math,os,platform,time
from pathlib import Path
import numpy as np
import scipy
from scipy.special import ndtr,ndtri
from scipy.stats import qmc

SQ=math.sqrt(2*math.pi)


def phi(x):
    return np.exp(-.5*x*x)/SQ


def evaluate(l,r,m,mu,steps=46):
    cl=ndtr(l); cr=ndtr(-r)
    pl=phi(l); pr=phi(r)
    qlo=np.zeros(len(l)); qhi=m.copy()
    for _ in range(steps):
        p=(qlo+qhi)*.5
        a=ndtri(cl+p)
        b=-ndtri(cr+(m-p))
        v=pl-phi(a)+phi(b)-pr
        choose=(v>mu)
        qlo=np.where(choose,p,qlo)
        qhi=np.where(choose,qhi,p)
    p=(qlo+qhi)*.5
    a=ndtri(cl+p)
    b=-ndtri(cr+(m-p))
    pa,pb=phi(a),phi(b)
    v=pl-pa+pb-pr
    second=m+l*pl-a*pa+b*pb-r*pr
    return a,b,p,v,second


def get_cases(shard,rep,exponent):
    u,v,w,z,h=qmc.Sobol(5,scramble=True,rng=np.random.default_rng(887452+shard+15*rep)).random_base2(exponent).T
    l=-5+7*(shard+u)/15
    width=np.exp2(-3+6*v)
    r=l+width
    reflect=l+r>0
    l,r=np.where(reflect,-r,l),np.where(reflect,-l,r)
    cl=ndtr(l);cr=ndtr(-r)
    total=1-cl-cr
    narrow=l>-2
    total[narrow]=ndtr(r[narrow])-cl[narrow]
    frac=.0005+.999*w
    frac[::23]=1e-6
    frac[1::23]=1-1e-6
    m=total*frac
    p=m*z
    a=ndtri(cl+p)
    b=-ndtri(cr+(m-p))
    mu=phi(l)-phi(a)+phi(b)-phi(r)
    return l,r,m,mu,a,b,p


def work(shard,i,exponent):
    l,r,m,mu,a0,b0,p0=get_cases(shard,i,exponent)
    a,b,p,v,second=evaluate(l,r,m,mu)
    an,bn,pn,vn,sn=evaluate(l[::37],r[::37],m[::37],mu[::37],54)
    gap=b-a
    finite=np.isfinite(second)&np.isfinite(a)&np.isfinite(b)
    assert np.all(finite) and np.all(gap>=-1e-12)
    return dict(points=len(l),nested=len(sn),max_mu_error=float(max(abs(v-mu))),
                max_mass_split_error=float(max(abs(p-p0))),
                max_second_nested=float(max(abs(second[::37]-sn))),
                max_absolute_threshold_error=float(max(np.maximum(abs(a-a0),abs(b-b0)))))


def main():
    p=argparse.ArgumentParser()
    p.add_argument('shard',type=int)
    p.add_argument('--exponent',type=int,default=11)
    p.add_argument('--seconds',type=float,default=90)
    p.add_argument('--out',type=Path,required=True)
    a=p.parse_args()
    start=time.monotonic();cpu=time.process_time();rows=[]
    while len(rows)<8192 and (not rows or time.monotonic()-start<a.seconds):
        rows.append(work(a.shard,len(rows),a.exponent))
    a.out.write_text(json.dumps(dict(shard=a.shard,pid=os.getpid(),elapsed=time.monotonic()-start,cpu=time.process_time()-cpu,
                        software=[platform.python_version(),np.__version__,scipy.__version__],results=rows),indent=2)+"\n")

if __name__=='__main__': main()
