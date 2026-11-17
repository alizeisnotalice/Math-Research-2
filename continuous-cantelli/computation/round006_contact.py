"""Local logarithmic contact integral: numerical evidence, not an exact solution."""
import os
for key in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS']: os.environ[key]='1'
import math,json,sys,platform,hashlib
from pathlib import Path
import scipy
from scipy.integrate import quad
import sympy as S
l,a,s=S.symbols('l a s',positive=True)
assert S.simplify(l*a-s-a*a*l*l/(4*s)+(s-a*l/2)**2/s)==0
rows=[];a=1.;t0=4.
for beta in [0.,.25,.5,1.]:
 for lam in [20.,100.,500.,2000.,10000.]:
  A=a*lam/2; center=A+beta*math.log(A);scale=math.sqrt(A)
  vals=[];errs=[]
  for sign in [-1,1]:
   def integrand(y):
    t=center+scale*y
    if t<t0:return 0.
    ss=t-beta*math.log(t);u=math.exp(-t);x=a+sign*u
    exponent=-(ss-A)**2/ss-beta*math.log(t)+sign*lam*u-x*x/2
    return math.exp(exponent)*scale/math.sqrt(2*math.pi)
   lo=(t0-center)/scale
   # Split into bounded saddle window and explicit adjacent tails.
   cuts=sorted(set([lo]+[v for v in [-12.,0.,12.] if v>lo]))
   total=err=0.
   for left,right in zip(cuts,cuts[1:]):
    v,e=quad(integrand,left,right,epsabs=1e-11,epsrel=1e-11,limit=250);total+=v;err+=e
   v,e=quad(integrand,cuts[-1],math.inf,epsabs=1e-11,epsrel=1e-11,limit=250);total+=v;err+=e
   vals.append(total);errs.append(err)
  target=2*math.exp(-a*a/2)/math.sqrt(2)*A**(.5-beta)
  rows.append(dict(beta=beta,lambda_value=lam,local_integral=sum(vals),asymptotic=target,ratio=sum(vals)/target,quad_error=sum(errs)))
out=dict(parameters=dict(a=a,t0=t0,betas=[0,.25,.5,1],lambda_values=[20,100,500,2000,10000],epsabs=1e-11,epsrel=1e-11),rows=rows,
 runtime=dict(python=platform.python_version(),scipy=scipy.__version__,sympy=S.__version__,pid=os.getpid()),
 exact_check='completed square symbolic identity passed',scope='local model only; quad errors not directed certificates; no global solution',script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
Path(sys.argv[1]).write_text(json.dumps(out,indent=2)+'\n')
for row in rows:print(row)
