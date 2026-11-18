"""Stable past heat-trace differences and near-quadratic transition."""
import argparse,json,os,math,time,platform
from pathlib import Path
import numpy as np
import scipy
from scipy.special import roots_legendre,erf,gammaln
from scipy.stats import qmc

rules={n:roots_legendre(n) for n in (16,32)}
LR=math.log(.25)
SQ=float(math.sqrt(2*math.pi))


def normalized(p,c,a,hlog,order):
    k=1-p/2
    s0=k*(LR-hlog)
    s1=k*(LR-(2*hlog-np.log(2*c))/p)
    f=np.array([-.5,-.25,-.125,-.06,-.025,0,.025,.06,.125,.25,.5])
    e=np.broadcast_to(np.array([0,.05,.2,.5,1,2,4,8,16,32,52]),(len(p),11))
    e=np.sort(np.clip(np.concatenate([e,s0[:,None]+f,s1[:,None]+f],axis=1),0,52),axis=1)
    dt=np.diff(e)/2;mid=(e[:,1:]+e[:,:-1])/2
    n,w=rules[order]
    s=mid[:,:,None]+dt[:,:,None]*n
    lu=LR-s/k[:,None,None]
    u=np.exp(lu)
    delta_h=lu-hlog[:,None,None]
    lm=(2*hlog)[:,None,None]-p[:,None,None]*lu-np.log(2*c)[:,None,None]
    signed=-np.expm1(np.log(2.)+delta_h)
    logabs=lm+np.log(np.maximum(abs(signed),1e-300))
    delta=np.exp(np.minimum(logabs,80.))
    A=np.exp((2-p)[:,None,None]*LR-2*s)/(2*c[:,None,None])
    low=np.where(signed>=0,A,A*np.square(-np.expm1(np.minimum(-delta_h,np.log(2.)))))
    rel=np.sign(signed)*(-np.expm1(-delta))*np.exp(-low-s+a[:,None,None]*u-u*u/2)
    return np.sum(w*dt[:,:,None]*rel,axis=(1,2))*np.exp(k*LR)


def work(shard,it,ex):
    u,v,w,q=qmc.Sobol(4,scramble=True,rng=np.random.default_rng(750363+shard+15*it)).random_base2(ex).T
    p=2-np.exp2(-11+10.9*(shard+u)/15)
    c=np.exp2(-2.5+4.5*v)
    a=-1+4*w
    hlog=-(8+155*q)*np.log(2.)
    n=0;mx=np.zeros(3)
    for j in range(0,len(p),384):
        i=slice(j,j+384);ii=slice(j,j+384,41)
        y=normalized(p[i],c[i],a[i],hlog[i],16)
        r=normalized(p[ii],c[ii],a[ii],hlog[ii],32)
        lim=np.sqrt(np.pi*c[i]/2)*erf(np.exp((2-p[i])*hlog[i]/2)/np.sqrt(2*c[i]))
        sel=p[i]>1.99
        q0=1/p[i]-.5
        cusp=np.exp(-q0*np.log(2*c[i])+gammaln(1-q0)+(2*q0)*hlog[i])
        interior=(p[i]>.95)&(p[i]<1.74)&(hlog[i]<-42)
        assert np.all(np.isfinite(y))
        mx=np.maximum(mx,[max(abs(r-y[::41])),max(abs(y[sel]-lim[sel])) if np.any(sel) else 0.,
                           max(abs(y[interior]/cusp[interior]-1)) if np.any(interior) else 0.])
        n+=len(r)
    return dict(profiles=len(p),independent=n,max_nested=float(mx[0]),
                max_near_quadratic=float(mx[1]),max_cusp_interior_relative=float(mx[2]))


def main():
    p=argparse.ArgumentParser();p.add_argument('shard',type=int);p.add_argument('--exponent',type=int,default=11)
    p.add_argument('--seconds',type=float,default=122);p.add_argument('--out',type=Path,required=True);a=p.parse_args()
    st=time.monotonic();cpu=time.process_time();rows=[]
    while len(rows)<4096 and (not rows or time.monotonic()-st<a.seconds):
        rows.append(work(a.shard,len(rows),a.exponent))
    a.out.write_text(json.dumps(dict(shard=a.shard,pid=os.getpid(),elapsed=time.monotonic()-st,cpu=time.process_time()-cpu,
        software=[platform.python_version(),np.__version__,scipy.__version__],results=rows),indent=2)+'\n')
if __name__=='__main__': main()
