"""Scaled circular heat Fourier coefficients for cusps and smooth minima."""
import argparse,json,math,os,time,platform
from pathlib import Path
import numpy as np
import scipy
from scipy.special import roots_legendre,gammaln
from scipy.stats import qmc

Knots=np.array([0,1e-6,1e-4,.002,.02,.10,.35,.75,1.5,3,6,12,24,48,96,160.])
Rules={i:roots_legendre(i) for i in (12,24)}


def integrate(k,j,c,a,order):
    logL=np.log(c*(k*j)**2/2)/a
    L=np.exp(logL)
    eta=j/L
    top=np.exp(np.minimum(15,np.log(np.pi**a)+a*logL))
    edges=np.minimum(Knots,top[:,None])
    width=np.diff(edges)/2
    midpoint=(edges[:,:-1]+edges[:,1:])/2
    nodes,weights=Rules[order]
    z=midpoint[:,:,None]+width[:,:,None]*nodes
    unit=(nodes+1)/2
    z[:,0,:]=edges[:,1,None]*unit**2
    jac=np.broadcast_to(width[:,:,None],z.shape).copy()
    jac[:,0,:]=edges[:,1,None]*unit
    z=np.maximum(z,np.finfo(float).tiny)
    t=np.exp(np.log(z)/a[:,None,None])
    u=t/L[:,None,None]
    sinc=np.sinc(u/(2*np.pi))
    atten=np.exp(-z*np.power(sinc,a[:,None,None]))
    base=np.exp(-z)
    cos=np.cos(eta[:,None,None]*t)
    prep=z**(1/a[:,None,None]-1)*cos*weights*jac
    integ=np.sum(prep*atten,axis=(1,2))/(np.pi*a)
    ref=np.sum(prep*base,axis=(1,2))/(np.pi*a)
    # The reference is explicitly integrated to infinity using the same finite kernel cut at z=160.
    if np.any(top<160):
        zf=(Knots[:-1]+Knots[1:])[:,None]/2+np.diff(Knots)[:,None]*nodes/2
        dzf=np.broadcast_to(np.diff(Knots)[:,None]/2,zf.shape).copy()
        zf[0,:]=Knots[1]*unit**2
        dzf[0,:]=Knots[1]*unit
        zf=zf[None,:,:]
        tf=np.exp(np.log(zf)/a[top<160,None,None])
        uf=eta[top<160,None,None]*tf
        ref[top<160]=np.sum(zf**(1/a[top<160,None,None]-1)*np.exp(-zf)*np.cos(uf)
                            *weights*dzf,axis=(1,2))/(np.pi*a[top<160])
    limit=np.exp(gammaln(1/a))/(np.pi*a)
    smooth=(a==2)
    limit[smooth]*=np.exp(-eta[smooth]**2/4)
    return integ,ref,limit,L,eta


def work(shard,rep,exponent):
    u,v,w,f=qmc.Sobol(4,scramble=True,rng=np.random.default_rng(636017+shard+15*rep)).random_base2(exponent).T
    k=np.floor(np.exp2(1+11*u)).astype(np.int64)
    j=np.floor(np.exp2(18*v)).astype(np.int64)
    k[::6]=2
    j[1::6]=1
    c=np.exp2(-2+4*w)
    a=.36+(2-.36)*(shard+f)/15
    a[::5]=2
    maxima=np.zeros(5)
    refs=0
    for base in range(0,len(k),480):
        i=slice(base,base+480)
        ii=slice(base,base+480,47)
        x,z,p,L,eta=integrate(k[i],j[i],c[i],a[i],12)
        y=integrate(k[ii],j[ii],c[ii],a[ii],24)[0]
        refs+=len(y)
        assert np.all(np.isfinite(x))
        # Use an unscaled empirical deviation. Scaled rates at very large L
        # mathematically vanish but subtraction has a finite precision floor.
        large=L>=128
        smooth=a[i]==2
        high=(a[i]<2)&(eta<1e-4)&large
        vals=[max(abs(x[::47]-y)),max(abs(x-z)),
              max(abs(x[large]-z[large])) if np.any(large) else 0,
              max(abs(x[smooth&large]-p[smooth&large])) if np.any(smooth&large) else 0,
              max(abs(x[high]-p[high])) if np.any(high) else 0]
        maxima=np.maximum(maxima,vals)
    return dict(points=len(k),nested=refs,max_nested=float(maxima[0]),max_geometry=float(maxima[1]),
                max_geometry_L128=float(maxima[2]),max_smooth_L128=float(maxima[3]),
                max_cusp_eta_small=float(maxima[4]))


def main():
    p=argparse.ArgumentParser();p.add_argument('shard',type=int);p.add_argument('--exponent',type=int,default=13)
    p.add_argument('--seconds',type=float,default=120);p.add_argument('--out',type=Path,required=True);a=p.parse_args()
    start=time.monotonic();cpu=time.process_time();rows=[]
    while len(rows)<4096 and (not rows or time.monotonic()-start<a.seconds):
        rows.append(work(a.shard,len(rows),a.exponent))
    a.out.write_text(json.dumps(dict(shard=a.shard,pid=os.getpid(),elapsed=time.monotonic()-start,cpu=time.process_time()-cpu,
        software=[platform.python_version(),np.__version__,scipy.__version__],results=rows),indent=2)+'\n')

if __name__=='__main__': main()
