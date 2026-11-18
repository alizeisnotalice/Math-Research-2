"""Search for Bernstein-weight negative heat remainder in the exact round 5 model."""
import argparse,json,os,platform,time
from pathlib import Path
import numpy as np
import scipy
from scipy.special import roots_legendre


def grid(y,n,model):
    a,w=roots_legendre(n)
    intervals=[(0.,.5),(.5,1.),(1.,1.5)] if model=='three' else [(0.,.5),(.5,1.)]
    times=np.concatenate([a*(hi-lo)/2+(hi+lo)/2 for lo,hi in intervals]); t=times[:,None]
    quad=np.concatenate([w*(hi-lo)/2 for lo,hi in intervals])
    def g(v,m):
        return np.exp(-(y[None,:]-m)**2/(2*v))/np.sqrt(2*np.pi*v)
    early=g(t[:n],0.)
    if model=='three':
        e=.6*np.exp(-16/3)
        vals=np.concatenate([early,g(t[n:2*n],0.)-e*g(t[n:2*n]-3/8,2.),e*g(t[2*n:]-7/8,2.)])
    else:
        vals=np.concatenate([early,g(t[n:],0.)-.55*g(t[n:]-3/8,0.)])
    return times,quad,vals


def stable_results(j,ks,t,w,values):
    ex=-j*t[None,:]+ks[:,None]*np.log(-np.expm1(-t))[None,:]
    ref=ex.max(axis=1)
    raw=(np.exp(ex-ref[:,None])*w[None,:])@values
    return raw,ref


def profile(shard,model):
    y=np.linspace(-1.,1.,769) if model=='two' else np.linspace(1.6,4.2,1153)
    k=np.arange(0,432)
    t,w,v=grid(y,96,model)
    tc,wc,vc=grid(y,192,model)
    minima=[]; n=0
    for j in range(shard,320,15):
        d,ref=stable_results(j,k,t,w,v)
        neg=np.nonzero(np.min(d,axis=1)<0)[0]
        if not len(neg):
            n+=len(k)*len(y)
            continue
        # Recheck first negative k for this j and strongest grid candidate.
        i0,_=np.unravel_index(np.argmin(d),d.shape)
        candidates=np.unique(np.r_[neg[:12],i0])
        dc,rc=stable_results(j,k[candidates],tc,wc,vc)
        for nrow,kk in enumerate(candidates):
            m=int(np.argmin(d[kk]))
            minima.append((int(j),int(kk),float(y[m]),float(d[kk,m]),float(dc[nrow,m]),float(rc[nrow]),float(d[kk,m]*np.exp(ref[kk]))))
        n+=len(k)*len(y)
    return {'model':model,'profiles':n,'y_range':[float(y[0]),float(y[-1])],'y_count':len(y),'tested_negative':minima}


def main():
    p=argparse.ArgumentParser();p.add_argument('index',type=int);p.add_argument('out',type=Path)
    args=p.parse_args();start=time.time();cpu=time.process_time()
    results=[profile(args.index,model) for model in ('three','two')]
    ans={'shard':args.index,'pid':os.getpid(),'start_epoch':start,'end_epoch':time.time(),
         'cpu_seconds':time.process_time()-cpu,'versions':[platform.python_version(),np.__version__,scipy.__version__],
         'results':results}
    args.out.write_text(json.dumps(ans,indent=2)+'\n')

if __name__=='__main__':main()
