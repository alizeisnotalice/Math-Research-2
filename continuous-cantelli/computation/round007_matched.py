"""Floating exploration of a profile with the required polynomial tail correction."""
import os
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS']:os.environ[k]='1'
import numpy as np,scipy,json,platform,sys,hashlib
from pathlib import Path
from scipy.integrate import quad
from scipy.optimize import brentq
k=.1;d=1/(1+k)
m=quad(lambda x:np.exp(-(1+k)*x*x/2)/(1+x*x)**1.5,-np.inf,np.inf,epsabs=1e-12,epsrel=1e-12)[0]/np.sqrt(2*np.pi)
def coeff(u):
 B=d/(m*np.sqrt(1+u)-1);A=B+d;c=A*u-B*(k+3);return A,B,c
def eq(u):
 A,B,c=coeff(u)
 return np.log(B)-np.log(2*np.sqrt(2/c))-1/c+2.5*np.log1p(k)
u=brentq(eq,40,100,xtol=1e-12);A,B,c=coeff(u)
def f(x):return 1-A*np.exp(-u*x*x/2)+B*(1+x*x)**(-1.5)*np.exp(-k*x*x/2)
rows=[]
for lam in [.5,1,2,4,8,12]:
 def integrand(x):return np.exp(lam*x-lam*lam*f(x)/2-x*x/2)/np.sqrt(2*np.pi)
 cuts=[-30,-2,-.25,0,.25,2,30];v=e=0.
 for lo,hi in zip(cuts,cuts[1:]):
  val,err=quad(integrand,lo,hi,epsabs=1e-12,epsrel=1e-12,limit=250);v+=val;e+=err
 rows.append(dict(lambda_value=lam,F=v,residual=v-1,quad_error=e))
out=dict(parameters=dict(kappa=k,u=u,A=A,B=B,alpha=1-A+B,curvature=c,m=m,root_bracket=[40,100],root_residual=eq(u)),rows=rows,runtime=dict(python=platform.python_version(),scipy=scipy.__version__,numpy=np.__version__,pid=os.getpid()),scope='floating profile and integral; asymptotic and mean constraints matched only numerically; not an exact solution or certified root',script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
Path(sys.argv[1]).write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
