"""Gap-sum Gibbs Fourier coefficients at Cantor-compatible harmonics."""
import argparse,json,math,os,time,platform
from pathlib import Path
import numpy as np
import scipy
from scipy.special import roots_legendre
from scipy.stats import qmc

RULES={r:tuple((q+1)/2 if i==0 else q/2 for i,q in enumerate(roots_legendre(r))) for r in (48,80)}
LEVEL_OFFSETS=np.arange(-45,83)
B_TABLE=np.arange(3,13,dtype=float)
SMALL_PHASE=np.cos(np.pi*(B_TABLE[:,None]-1)*B_TABLE[:,None]**(-np.arange(1.,129.)))
PRODUCT=np.concatenate((np.ones((10,1)),np.cumprod(SMALL_PHASE,axis=1)),axis=1)


def coefficient(b,k,a,c,order):
    logb=np.log(b)
    logn=k*logb
    logeps=-(np.log(c)+2*logn)/a
    scale=np.exp(logn+logeps)
    logH=np.log((b-2)/2)
    main=np.floor((logH-logeps)/logb).astype(int)
    ell=main[:,None]+LEVEL_OFFSETS
    weight=np.where(ell>=1,np.exp2(LEVEL_OFFSETS)[None,:],0.)
    h=logH[:,None]-ell*logb[:,None]
    rr=np.exp(np.minimum(h-logeps[:,None],math.log(120)/a[:,None]))
    s,w=RULES[order]
    power=2/a[:,None,None]
    t=rr[:,:,None]*s**power
    density=np.exp(-t**a[:,None,None])*rr[:,:,None]*power*s**(power-1)*w
    late=ell-k[:,None]
    center=np.where(late>0,-PRODUCT[b[:,None]-3,np.clip(late-1,0,128)],1.)
    # The early-gap center sign cancels its exact half-integer phase.
    offset=np.where(late>0,np.exp(logH[:,None]-np.maximum(late,1)*logb[:,None]),0.)
    phase=np.cos(2*np.pi*(scale[:,None,None]*t-offset[:,:,None]))*center[:,:,None]
    mass=weight*np.sum(density,axis=2)
    line=np.sum(weight*np.sum(density*phase,axis=2),axis=1)
    z=np.sum(mass,axis=1)
    cutoff_mass=np.sum(np.where(late<=0,mass,0.),axis=1)
    assert np.all(z>0)
    return line/z,-PRODUCT[b-3,-1],cutoff_mass/z,scale


def work(shard,rep,power):
    u,v,w,q=qmc.Sobol(4,scramble=True,rng=np.random.default_rng(682413+shard+15*rep)).random_base2(power).T
    b=3+np.floor(10*u).astype(int)
    k=2+np.floor(23*v).astype(int)
    a=.45+1.55*(shard+w)/15
    a[::5]=2
    c=np.exp2(-3+6*q)
    maxima=np.zeros(4)
    nested=0
    for lo in range(0,len(b),22):
        ii=slice(lo,lo+22); jj=slice(lo,lo+22,7)
        coeff,limit,early,rs=coefficient(b[ii],k[ii],a[ii],c[ii],48)
        ref=coefficient(b[jj],k[jj],a[jj],c[jj],80)[0]
        nested+=len(ref)
        d=np.log(2)/np.log(b[ii])
        sub=rs<.1
        dscaled=np.abs(coeff-limit)/rs**d
        ave=np.abs(coeff-limit)
        maxima=np.maximum(maxima,[max(abs(coeff[::7]-ref)),max(early/rs**d),
                                  max(dscaled[sub]) if np.any(sub) else 0.,
                                  max(ave[rs<1e-4]) if np.any(rs<1e-4) else 0.])
        assert np.all(np.isfinite(coeff))
    return dict(evaluated=len(b),nested=nested,max_nested=float(maxima[0]),
                max_exterior_scaled=float(maxima[1]),max_subcritical_scaled=float(maxima[2]),
                max_absolute_eta_small=float(maxima[3]))


def main():
    p=argparse.ArgumentParser();p.add_argument('shard',type=int);p.add_argument('--exponent',type=int,default=8)
    p.add_argument('--seconds',type=float,default=100);p.add_argument('--out',type=Path,required=True);a=p.parse_args()
    nt=time.monotonic();ct=time.process_time();rows=[]
    while len(rows)<4096 and (not rows or time.monotonic()-nt<a.seconds):
        rows.append(work(a.shard,len(rows),a.exponent))
    a.out.write_text(json.dumps(dict(shard=a.shard,pid=os.getpid(),elapsed=time.monotonic()-nt,cpu=time.process_time()-ct,
                  software=[platform.python_version(),np.__version__,scipy.__version__],results=rows),indent=2)+'\n')

if __name__=='__main__': main()
