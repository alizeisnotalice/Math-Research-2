"""Smooth compact perturbations fitted to finite constraints; numerical roots only."""
import os
for key in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS']:os.environ[key]='1'
from pathlib import Path
import numpy as np,scipy,json,sys,platform,hashlib
from numpy.polynomial.legendre import leggauss
from scipy.optimize import root
centers=np.array([.5,1.1,1.7,2.3,2.9]);width=.2;lams=np.array([1.,2.]);rows=[]
def setup(order):
 y,w=leggauss(order);x=centers[:,None]+width*y;psi=np.exp(1-1/(1-y*y));mass=2*width*w*np.exp(-x*x/2)/np.sqrt(2*np.pi)
 return x,psi,mass
x,psi,mass=setup(128)
def evaluate(a,x=x,psi=psi,mass=mass,params=lams):
 b=a[:,None]*psi
 mean=np.sum(mass*b);fourth=np.sum(mass*(b*b-2*(x*x-1)*b))
 vals=[np.sum(mass*np.exp(-lam*lam/2)*np.cosh(lam*x)*np.expm1(-lam*lam*b/2)) for lam in params]
 return np.array([mean,fourth]+vals)
def jac(a):
 b=a[:,None]*psi
 j=[np.sum(mass*psi,axis=1),np.sum(mass*2*(b-(x*x-1))*psi,axis=1)]
 for lam in lams:j.append(np.sum(mass*np.exp(-lam*lam/2)*np.cosh(lam*x)*(-lam*lam/2)*np.exp(-lam*lam*b/2)*psi,axis=1))
 return np.array(j)[:,:4]
a=np.zeros(4)
for free in [.01,.03,.1]:
 sol=root(lambda z:evaluate(np.r_[z,free]),a,jac=lambda z:jac(np.r_[z,free]),tol=1e-11)
 a=sol.x;full=np.r_[a,free];xx,pp,mm=setup(256)
 checks=evaluate(full,xx,pp,mm,np.array([.5,1,1.5,2,2.5,3,4,5,6]))
 rows.append(dict(free=free,coefficients=full.tolist(),solver_success=bool(sol.success),solver_message=sol.message,nfev=int(sol.nfev),fitted_residuals=evaluate(full).tolist(),validation_256=dict(mean=checks[0],fourth_condition=checks[1],parameters=[.5,1,1.5,2,2.5,3,4,5,6],F_minus_one=checks[2:].tolist()),f_bounds=[1+min(0,float(min(full))),1+max(0,float(max(full)))],lipschitz_upper_formula_value=float(8/np.e/width*max(abs(full)))))
out=dict(parameters=dict(centers=centers.tolist(),width=width,fit_params=lams.tolist(),fit_order=128,validation_order=256,free_values=[.01,.03,.1]),rows=rows,linear_jacobian=jac(np.zeros(5)).tolist(),linear_condition_number=float(np.linalg.cond(jac(np.zeros(5)))),runtime=dict(python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__,pid=os.getpid()),script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),scope='floating roots and quadrature; exact finite-constraint existence proved separately by implicit-function argument')
Path(sys.argv[1]).write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
