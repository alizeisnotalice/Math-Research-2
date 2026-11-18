"""Stable quartic log-amplitude asymmetric Wick coefficient series and independent finite convolution."""
import argparse,json,os,time,platform
from pathlib import Path
import numpy as np
import scipy
from scipy.special import gammaln,logsumexp
from scipy.stats import qmc


def batch(shard,rep,exponent):
    grid=qmc.Sobol(7,scramble=True,rng=np.random.default_rng(701030+shard+15*rep)).random_base2(exponent)
    a,b,c,p,q,z,j4=grid.T
    m=11+np.floor(8*(shard+a)/15).astype(int)
    d=m-1
    n=2*np.floor(np.exp2(5+5*b)).astype(int)
    frac=(.10+.80*c)/(2*d)
    k=np.maximum(1,np.floor(n*frac).astype(int))
    nu=k/n
    M=n//2-d*k
    assert np.all(M>=0)
    t=2*m*k
    lam=d/(2*m-1)*(4*m/(2*m-1))**(m/d)*(1+.03*z)
    alpha=(2*p-1)*1.4
    beta=.5+1.75*q
    gamma=.025+2.975*j4
    center=-.8+1.6*z
    scale=n**(-1/(2*d))
    A=np.abs((2*k*alpha-center)*scale)
    B=(2*k*beta-.5)*scale*scale
    C=2*k*gamma*scale**4
    lc=np.log(4*C)
    assert np.all(B>=0)
    la=np.log(np.maximum(A,1e-300))
    lb=np.log(np.maximum(2*B,1e-300))
    score=np.zeros(len(n))
    term=np.zeros(len(n))
    prev=np.zeros(len(n))
    prevprev=np.full(len(n),-np.inf)
    tail1=prevprev.copy()
    tail2=prevprev.copy()
    # Exhaustive through the greatest finite index, masked at each case's M.
    for j in range(1,2*int(max(M))+1):
        new=np.logaddexp(np.logaddexp(la+prev,lb+prevprev),lc+tail2)-np.log(j)
        tail2,tail1,prevprev,prev=tail1,prevprev,prev,new
        if j%2==0:
            ell=j//2
            with np.errstate(divide='ignore',invalid='ignore'):
                mult=np.log(np.maximum(M-ell+1,1))+np.log(2/(lam*n))
            mult+=np.log(t+j-1)+np.log(t+j)-np.log(t-2*k+j-1)-np.log(t-2*k+j)
            term+=mult
            score=np.logaddexp(score,np.where(ell<=M,new+term,-np.inf))
    verify=[]
    for i in range(0,len(n),61):
        ell=np.arange(M[i]+1)
        l=ell[:,None]
        v=ell[None,:]
        u=2*(l-v)
        coef=logsumexp(np.where(u>=0,u*la[i]+v*(lb[i]-np.log(2))
                               -gammaln(u+1)-gammaln(v+1),-np.inf),axis=1)
        src=ell[:,None]-2*ell[None,:]
        coef=logsumexp(np.where(src>=0,coef[np.maximum(src,0)]+ell[None,:]*np.log(C[i])-gammaln(ell[None,:]+1),-np.inf),axis=1)
        top=t[i]+2*ell
        comb=gammaln(top+1)-gammaln(top-2*k[i]+1)
        logw=comb-comb[0]+gammaln(M[i]+1)-gammaln(M[i]-ell+1)+ell*np.log(2/(lam[i]*n[i]))
        verify.append(abs(float(logsumexp(coef+logw))-score[i]))
    assert np.all(np.isfinite(score))
    predicted=2*nu*np.abs(alpha)*(m/d)*np.sqrt((1-2*d*nu)/lam)
    normalized=score/(n*scale)
    target=(np.abs(alpha)>.8)&(m>=14)&(n>=512)
    return dict(points=len(n),independent=len(verify),max_independent_difference=max(verify),
                log_group_ratio=[float(min(score)),float(max(score))],
                scale_normalized_asymmetric_cases=int(np.count_nonzero(target)),
                max_scale_relaxation=float(max(abs(normalized[target]-predicted[target]))) if np.any(target) else 0.)


def main():
    p=argparse.ArgumentParser();p.add_argument('shard',type=int);p.add_argument('--exponent',type=int,default=9)
    p.add_argument('--seconds',type=float,default=90);p.add_argument('--out',type=Path,required=True);a=p.parse_args()
    start=time.monotonic();cpu=time.process_time();rows=[]
    while len(rows)<4096 and (not rows or time.monotonic()-start<a.seconds):
        rows.append(batch(a.shard,len(rows),a.exponent))
    a.out.write_text(json.dumps(dict(shard=a.shard,pid=os.getpid(),elapsed=time.monotonic()-start,cpu=time.process_time()-cpu,
        software=[platform.python_version(),np.__version__,scipy.__version__],results=rows),indent=2)+'\n')

if __name__=='__main__': main()
