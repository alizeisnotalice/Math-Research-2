import argparse,json,math,os,time,platform
from pathlib import Path
import numpy as np
import scipy
from scipy.special import roots_legendre
from scipy.stats import qmc
k=.5*math.log(2)
sqrt2pi=math.sqrt(2*math.pi)
meshes={n:roots_legendre(n) for n in (64,96)}

def model(e,j,q,r):
    l=np.log(e); a=j*l+q; b=r*l*l
    o=.5*l+k+a/e+b/e**2
    return e+o,o,1+.5/e+(j-a)/e**2+(2*r*l-2*b)/e**3

def calc(par,second=False):
    L,j,q,r=par.T
    e=L-np.log(2*L)
    for _ in range(4):
        D,o,d=model(e,j,q,r)
        e-=(e-L+np.log(2*D))/(1+d/D)
    D,o,d=model(e,j,q,r)
    raw=np.sqrt(2*D)*np.exp(o-.5*np.exp(-2*e))/(2*(D+d))
    LL,jj,qq,rr=[p[:,None] for p in (L,j,q,r)]
    rt=np.sqrt(D)[:,None]
    n,w=meshes[96 if second else 64]
    z=(9 if second else 8)*n
    y=2*D[:,None]+2*rt*z
    if second:
        for _ in range(4):
            f,a,b=model(LL-np.log(y),jj,qq,rr); g=np.sqrt(f)
            y-=(y-2*f-2*g*z)/(1+b*(2+z/g)/y)
        f,a,b=model(LL-np.log(y),jj,qq,rr); g=np.sqrt(f)
        jac=2*g/(1+b*(2+z/g)/y)
    else:
        f,a,b=model(LL-np.log(y),jj,qq,rr)
        z=(y-2*f)/(2*np.sqrt(f)); jac=2*rt
    val=(9 if second else 8)*np.sum(w*jac*np.exp(a-np.log(y)-z*z-.5*(y*np.exp(-LL))**2),axis=1)/sqrt2pi
    return raw,val,D

def main():
    p=argparse.ArgumentParser();p.add_argument('shard',type=int);p.add_argument('--exponent',type=int,default=14)
    p.add_argument('--seconds',type=float,default=41);p.add_argument('--out',type=Path,required=True);a=p.parse_args()
    start=time.monotonic();cpu=time.process_time();rows=[]
    while len(rows)<4096 and (not rows or time.monotonic()-start<a.seconds):
        l,j,q,r=qmc.Sobol(4,scramble=True,rng=np.random.default_rng(631726+a.shard+15*len(rows))).random_base2(a.exponent).T
        par=np.column_stack((np.exp2(7+12*(a.shard+l)/15),-.3+1.1*j,-2+4*q,-1+2*r))
        mx=0.;mn=np.inf;hi=-np.inf;count=0
        for off in range(0,len(l),1024):
            b=par[off:off+1024]
            raw,val,d=calc(b); check=calc(b[::29],True)[1]
            mx=max(mx,float(max(abs(val[::29]-check))));count+=len(check)
            second=d*(2*d*(val/raw-1)-1)
            assert np.all(np.isfinite(second))
            mn=min(mn,float(min(second)));hi=max(hi,float(max(second)))
        rows.append(dict(points=len(l),crosschecks=count,maximum_abs_disagreement=mx,second_deficit_range=[mn,hi]))
    a.out.write_text(json.dumps(dict(shard=a.shard,pid=os.getpid(),elapsed=time.monotonic()-start,cpu=time.process_time()-cpu,software=[platform.python_version(),np.__version__,scipy.__version__],results=rows),indent=2)+'\n')
if __name__=='__main__': main()
