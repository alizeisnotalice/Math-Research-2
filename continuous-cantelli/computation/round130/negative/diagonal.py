"""Entire integer diagonal at 948; numerical cross-check, not a certificate."""
import json
from pathlib import Path
import numpy as np
from scipy.special import roots_legendre, xlogy
from scipy.optimize import minimize_scalar


def ratios(y,n,k,order):
    a,w=roots_legendre(order)
    times=(a+1)/4+np.array([0,.5,1])[:,None]
    diffuse=np.exp(-y[None,None,:]**2/(2*times[:2,:,None]))/np.sqrt(2*np.pi*times[:2,:,None])
    source_v=times[1:,:,None]-np.array([.375,.875])[:,None,None]
    sources=.6*np.exp(-16/3)*np.exp(-(y[None,None,:]-2)**2/(2*source_v))/np.sqrt(2*np.pi*source_v)
    times=times.ravel()
    j=n-k
    peak=(xlogy(j,j/n)+xlogy(k,k/n))[:,None]
    W=np.exp(-j[:,None]*times+k[:,None]*np.log(-np.expm1(-times))-peak)*np.tile(w/4,3)[None,:]
    P=W[:,:2*order]@diffuse.reshape(-1,len(y))+W[:,2*order:]@sources[1]
    N=W[:,order:2*order]@sources[0]
    return (P-N)/(P+N)


def main():
    n=948; k=np.arange(n+1); y=np.linspace(55/24,3.5,18303)
    local_min=np.full(n+1,np.inf); arguments=np.zeros(n+1)
    for offset in range(0,len(y),1200):
        x=y[offset:offset+1200]
        Q=ratios(x,n,k,128)
        ix=np.argmin(Q,axis=1)
        r=Q[np.arange(len(k)),ix]
        improve=r<local_min
        local_min[improve]=r[improve]
        arguments[improve]=x[ix[improve]]
    checks=[]
    for i in np.argsort(local_min)[:16]:
        x=arguments[i]
        def fun(p): return float(ratios(np.array([p]),n,np.array([i]),256)[0,0])
        res=minimize_scalar(fun,bounds=(x-.001,x+.001),method='bounded',options={'xatol':1e-10})
        checks.append((n-int(i),int(i),float(res.x),float(res.fun),float(local_min[i])))
    out=dict(diagonal=n,y_count=len(y),grid_y_domain=[55/24,3.5],profiles=len(y)*len(k),
             all_grid_positive=bool(np.all(local_min>0)), weakest_refined=checks,orders=[128,256])
    Path('diagonal-948.json').write_text(json.dumps(out,indent=2)+'\n')
    print(out)

if __name__=='__main__':main()
