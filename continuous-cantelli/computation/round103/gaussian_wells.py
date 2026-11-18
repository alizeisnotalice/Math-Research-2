"""Global Gaussian test lower envelope for three separated variance wells."""
import argparse,json,math,os,time,platform
from pathlib import Path
import numpy as np
import scipy
from scipy.special import ndtr
from scipy.stats import qmc
from scipy.integrate import quad


def p(s,d):
    return np.exp(-d*d/(2*s))/np.sqrt(2*np.pi*s)


def product(y,s,l,r):
    var=s/(1+s)
    z=y/(1+s)
    root=np.sqrt(var)
    return p(1+s,y)*np.maximum(0,ndtr((r-z)/root)-ndtr((l-z)/root))


def envelope(y,e,v,l,r):
    q=np.sqrt(e*(e+v)*np.log1p(v/e)/v)
    lo=y-q;hi=y+q
    return (product(y,e,l,np.minimum(r,lo))+
            product(y,e+v,np.maximum(l,lo),np.minimum(r,hi))+
            product(y,e,np.maximum(l,hi),r))


def exactcheck(row):
    M,e,y,wells=row['M'],row['epsilon'],row['y'],row['wells']
    def piece(l,r,v):
        A=math.sqrt(e*(e+v)*math.log1p(v/e)/v)
        pts=sorted([l,r]+[u for u in (y-A,y+A) if l<u<r])
        return sum(quad(lambda x:p(1.,x)*min(p(e,y-x),p(e+v,y-x)),
                        a,b,epsabs=2e-13)[0] for a,b in zip(pts,pts[1:]))
    base=piece(-np.inf,np.inf,M)
    value=base+sum(piece(l,r,b)-piece(l,r,M) for l,r,b in wells)
    return abs(value-row['lower_density']),value/p(M+e,y)


def work(shard,rep,k):
    S=qmc.Sobol(12,scramble=True,rng=np.random.default_rng(1012887+shard+15*rep)).random_base2(k)
    M=1.04+3.2*(shard+S[:,0])/15
    e=np.exp2(-5+5.5*S[:,1])
    y=-2.2+4.4*S[:,2]
    shift=-1+2*S[:,3]
    W=.35+3.8*S[:,4]
    arcs=.035+S[:,5:10]
    arcs=arcs/arcs.sum(axis=1)[:,None]*W[:,None]
    ends=y[:,None]+shift[:,None]-W[:,None]/2+np.concatenate([np.zeros((len(M),1)),np.cumsum(arcs,axis=1)],axis=1)
    L=ends[:,::2];R=ends[:,1::2]
    depth=np.column_stack((S[:,10],S[:,11],(S[:,10]+S[:,11])%1))
    upper=(.005+.39*depth)*M[:,None]
    mass=ndtr(R)-ndtr(L)
    eligible=((M[:,None]-upper)*mass).sum(axis=1)<.985
    Y=y[:,None];E=e[:,None]
    target=p(M+e,y)
    base=envelope(y,e,M,-np.inf,np.inf)
    uplift=envelope(Y,E,upper,L,R)-envelope(Y,E,M[:,None],L,R)
    lower=base+uplift.sum(axis=1)
    most=base+uplift.max(axis=1)
    cross=eligible&(lower>target)&(most<target)
    score=np.where(cross,lower/target,0.)
    ix=int(np.argmax(score))
    assert np.all(uplift>=-2e-15)
    record=None
    if score[ix]>0:
        record=dict(M=float(M[ix]),epsilon=float(e[ix]),y=float(y[ix]),wells=np.column_stack((L[ix],R[ix],upper[ix])).tolist(),
                    lower_density=float(lower[ix]),ratio=float(score[ix]),one_well_ratio=float(most[ix]/target[ix]),
                    mean_guard=float(((M[ix]-upper[ix])*mass[ix]).sum()))
    return dict(points=len(M),feasible_budget=int(np.count_nonzero(eligible)),joint_witness_count=int(np.count_nonzero(cross)),
                maximum_joint_ratio=float(score[ix]),record=record)


def main():
    p=argparse.ArgumentParser();p.add_argument('shard',type=int);p.add_argument('--exponent',type=int,default=12)
    p.add_argument('--seconds',type=float,default=75);p.add_argument('--out',type=Path,required=True);a=p.parse_args()
    st=time.monotonic();cpu=time.process_time();res=[];top=None
    while len(res)<8192 and (not res or time.monotonic()-st<a.seconds):
        row=work(a.shard,len(res),a.exponent)
        now=row.pop("record")
        if now is not None and (top is None or now['ratio']>top['ratio']):top=now
        res.append(row)
    if top is not None:
        top['independent_quad_error'],top['independent_ratio']=exactcheck(top)
    a.out.write_text(json.dumps(dict(shard=a.shard,pid=os.getpid(),elapsed=time.monotonic()-st,
        cpu=time.process_time()-cpu,software=[platform.python_version(),np.__version__,scipy.__version__],
        results=res,top=top),indent=2)+"\n")
if __name__=='__main__':main()
