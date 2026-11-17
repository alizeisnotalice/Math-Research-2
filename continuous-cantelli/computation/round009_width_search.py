"""Restricted width-family fit, negative numerical evidence only."""
import os
for key in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS']:os.environ[key]='1'
from pathlib import Path
import numpy as np,scipy,json,platform,sys,hashlib
from scipy.integrate import quad
from scipy.optimize import least_squares
rows=[]
for k in [.02,.05,.1,.2,.5,1.]:
 def fun(z,detail=False):
  s=np.exp(z[1]);u=k+5*s+np.exp(z[0]);g=1+k
  m=quad(lambda x:np.exp(-g*x*x/2)/(1+s*x*x)**1.5,0,np.inf,epsabs=1e-10)[0]*np.sqrt(2/np.pi)
  B=1/g/(m*np.sqrt(1+u)-1);A=B+1/g;c=A*u-B*(k+3*s);d4=-3*A*u*u+B*(3*k*k+18*k*s+45*s*s)
  s1=-1/c-2/c**2-d4*(1/(4*c*c)+1/c**3+1/(3*c**4));q=(s1-6*g)/g**2
  eq=[np.log(B*s**-1.5)-np.log(2*np.sqrt(2/c))-1/c+2.5*np.log(g),s*q+1.5]
  return dict(s=s,u=u,A=A,B=B,c=c,residuals=eq) if detail else eq
 sol=least_squares(fun,np.log([60,.3]),bounds=([-12,-12],[15,12]),max_nfev=500,gtol=1e-12,xtol=1e-12,ftol=1e-12)
 rows.append(dict(kappa=k,nfev=sol.nfev,status=int(sol.status),message=sol.message,variables=sol.x.tolist(),result=fun(sol.x,True)))
out=dict(rows=rows,runtime=dict(python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__,pid=os.getpid()),script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),scope='single-start bounded restricted-family search; boundary saturation and nonzero residuals do not prove infeasibility')
Path(sys.argv[1]).write_text(json.dumps(out,indent=2)+'\n');print('completed',len(rows),'restricted fits')
