"""Exploratory continuous graph fitting. No exact Cantelli counterexample claim."""
import os
for key in ['OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS']:os.environ[key]='1'
import numpy as np
import scipy
from scipy.special import ndtr, erfc
from scipy.optimize import minimize
from numpy.polynomial.legendre import leggauss
from pathlib import Path
import json,sys,time,platform,hashlib
L=float(sys.argv[2]) if len(sys.argv)>2 else 8.
Tmax=float(sys.argv[3]) if len(sys.argv)>3 else 7.
starts=int(sys.argv[4]) if len(sys.argv)>4 else 4
maxiter=int(sys.argv[5]) if len(sys.argv)>5 else 1200
H=3.;knots=np.linspace(0,L,int(4*L)+1);fixed={0:.6,1:.6,len(knots)-1:1.}
free=np.array([i for i in range(len(knots)) if i not in fixed])
base=np.zeros(len(knots))
for i,v in fixed.items():base[i]=v

def quadrature(order):
 z,w=leggauss(order);xx=[];ww=[];BB=[]
 for i in range(len(knots)-1):
  left,right=knots[i:i+2];x=(left+right)/2+(right-left)*z/2
  weight=w*(right-left)/2*2*np.exp(-x*x/2)/np.sqrt(2*np.pi)
  b=np.zeros((order,len(knots)));b[:,i]=(right-x)/(right-left);b[:,i+1]=(x-left)/(right-left)
  xx.extend(x);ww.extend(weight);BB.extend(b)
 return np.array(xx),np.array(ww),np.array(BB)

def setup(order,ts,xis):
 x,w,B=quadrature(order);bf=B[:,free];baseq=B@base
 meanbasis=w@B;meantail=2*ndtr(-L)
 A=ts[:,None]*x;logcw=np.logaddexp(A,-A)-np.log(2)+np.log(w)[None,:]
 tail=ndtr(ts-L)+ndtr(-ts-L);st=ts*ts/2
 cosw=np.cos(xis[:,None]*x)*w;sc=xis*xis/2
 target=np.exp(-sc*H);ctail=target*np.real(erfc((L-1j*xis)/np.sqrt(2)))
 scale=5*np.exp(sc*.5)
 def evaluate(z):
  fq=baseq+bf@z
  R=np.exp(logcw-st[:,None]*fq);F=R.sum(axis=1)+tail
  real=np.log(F);jr=-st[:,None]*(R@bf)/F[:,None]
  K=cosw*np.exp(-sc[:,None]*(H-fq));cf=K.sum(axis=1)+ctail-target
  jc=sc[:,None]*(K@bf)
  return np.r_[real,cf*scale],np.vstack([jr,jc*scale[:,None]]),F,cf
 return evaluate,meanbasis,meantail

ts=np.linspace(.25,Tmax,int(6*Tmax)+3);xis=np.linspace(.25,4.,24)
evaluate,meanbasis,meantail=setup(16,ts,xis)
coef=meanbasis[free];rhs=1.-meantail-meanbasis@base
constant=rhs/coef.sum()
rng=np.random.default_rng(20261117)
trials=[];started=time.monotonic()
for seed in range(starts):
 z0=np.full(len(free),constant)
 if seed:z0+=rng.normal(0,.15,len(free));z0+=(rhs-coef@z0)/coef.sum()
 cache={}
 def get(z):
  if 'z' not in cache or not np.array_equal(z,cache['z']):
   r,J,F,cf=evaluate(z);cache.update(z=z.copy(),value=float(r@r),grad=2*J.T@r)
  return cache['value'],cache['grad']
 opt=minimize(get,z0,jac=True,method='SLSQP',bounds=[(.05,2.8)]*len(free),
              constraints=[{'type':'eq','fun':lambda z:coef@z-rhs,'jac':lambda z:coef}],
              options={'maxiter':maxiter,'ftol':1e-18,'disp':False})
 values=base.copy();values[free]=opt.x
 tr=dict(seed=seed,success=bool(opt.success),message=str(opt.message),iterations=int(opt.nit),
         objective=float(opt.fun),mean_residual=float(coef@opt.x-rhs),values=values.tolist())
 trials.append(tr);print(json.dumps({k:v for k,v in tr.items() if k!='values'}),flush=True)
best=min(trials,key=lambda r:r['objective']);values=np.array(best['values']);z=values[free]
ht=np.linspace(.013,max(12.,Tmax+4),601);hx=np.linspace(.017,6.,401)
checks=[]
for order in [32,64]:
 ev,mb,mt=setup(order,ht,hx);rr,J,F,cf=ev(z)
 checks.append(dict(order=order,mean_f=float(mb@values+mt),max_real_log_residual=float(np.max(np.abs(np.log(F)))),
                    worst_real_t=float(ht[np.argmax(np.abs(np.log(F)))]),
                    max_characteristic_residual=float(np.max(np.abs(cf))),
                    worst_characteristic_xi=float(hx[np.argmax(np.abs(cf))]),
                    real_log_residual=np.log(F).tolist(),characteristic_residual=cf.tolist()))
result=dict(parameters=dict(L=L,H=H,knots=knots.tolist(),fixed_plateau=[0,.25,.6],bounds=[.05,2.8],
                            train_t=ts.tolist(),train_xi=xis.tolist(),seed=20261117,starts=starts,
                            quadrature_order=16,maxiter=maxiter,ftol=1e-18),
            trials=trials,best=best,validation_t=ht.tolist(),validation_xi=hx.tolist(),validation=checks,
            best_max_slope=float(np.max(np.abs(np.diff(values)/np.diff(knots)))),
            runtime=dict(python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__,pid=os.getpid(),
                         elapsed_seconds=time.monotonic()-started),
            script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
            scope='finite continuous piecewise-linear fit, f=1 outside [-L,L]; not an exact solution; no interval certificate')
Path(sys.argv[1]).write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps({'best_objective':best['objective'],'best_max_slope':result['best_max_slope'],
 'checks':[{k:v for k,v in r.items() if not isinstance(v,list)} for r in checks], 'elapsed':result['runtime']['elapsed_seconds']},indent=2),flush=True)
