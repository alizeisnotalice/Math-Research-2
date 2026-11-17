"""Smooth even microstructure plus finite-dimensional correction; numerical evidence only."""
import os
for key in ['OPENBLAS_NUM_THREADS','OMP_NUM_THREADS']:os.environ[key]='1'
from pathlib import Path
import numpy as np,scipy,json,sys,platform,hashlib,time
from scipy.special import expit
from scipy.optimize import root
from numpy.polynomial.legendre import leggauss
EPS=.6*np.exp(-16/3);R=14.;centers=np.array([.5,1.1,1.7,2.3]);width=.2

def segments(n):
 h=1/n;delta=1/(8*n*n);out=[]
 for k in range(int(R*n)):
  left=k/n;c=(k+.5)/n
  def rho(v):return EPS/2/np.sqrt(v)*(np.exp(-(c-2)**2/(2*v)+c*c/2)+np.exp(-(c+2)**2/(2*v)+c*c/2))
  rs,rt=rho(.125),rho(.625);rm=1-rs-rt;assert rm>0
  lens=np.array([delta,(1-4*delta)*rs,delta,(1-4*delta)*rm,delta,(1-4*delta)*rt,delta])*h/2
  ends=[(1,.5),(.5,.5),(.5,1),(1,1),(1,1.5),(1.5,1.5),(1.5,1)]
  x=left;half=[]
  for j,(length,(fa,fb)) in enumerate(zip(lens,ends)):
   b=c if j==6 else x+length
   if b>x:half.append((x,b,fa,fb))
   x=b
  out.extend(half+[(2*c-b,2*c-a,fb,fa) for a,b,fa,fb in reversed(half)])
 # Split quadrature only, retaining original smooth-ramp coordinates.
 pieces=[];pin=n**-3
 for a,b,fa,fb in out:
  cuts=[a]+[q*pin/4 for q in range(1,5) if a<q*pin/4<b]+[b]
  pieces.extend([(l,r,a,b,fa,fb) for l,r in zip(cuts,cuts[1:])])
 return np.array(pieces)

def setup(n,order):
 seg=segments(n);l,r,a,b,fa,fb=seg.T;y,w=leggauss(order)
 x=(l[:,None]+r[:,None])/2+(r-l)[:,None]*y/2
 z=np.clip((x-a[:,None])/(b-a)[:,None],1e-15,1-1e-15)
 step=expit(-1/z+1/(1-z));h=fa[:,None]+(fb-fa)[:,None]*step
 pin=x*n**3;chi=np.zeros_like(pin);take=pin<1;chi[take]=np.exp(1-1/(1-pin[take]**2))
 h=(1-chi)*h+.5*chi
 mass=(r-l)[:,None]*w*np.exp(-x*x/2)/np.sqrt(2*np.pi)
 x=x.ravel();h=h.ravel();mass=mass.ravel();tau=16*(h-.5)**2*(h-1.5)**2
 v=(x[None,:]-centers[:,None])/width;psi=np.zeros_like(v);take=abs(v)<1;psi[take]=np.exp(1-1/(1-v[take]**2));W=psi*tau
 return x,h,mass,W

def calc(a,data,params=[1.,2.]):
 x,h,mass,W=data;f=h+a@W;b=f-1
 rows=[mass@b,mass@(b*b-2*(x*x-1)*b)]
 rows.extend(mass@(np.exp(-lam*lam/2)*np.cosh(lam*x)*np.expm1(-lam*lam*b/2)) for lam in params)
 return np.array(rows)
def jac(a,data):
 x,h,mass,W=data;b=h+a@W-1
 rows=[W@mass,W@(mass*2*(b-(x*x-1)))]
 rows.extend(W@(mass*np.exp(-lam*lam/2)*np.cosh(lam*x)*(-lam*lam/2)*np.exp(-lam*lam*b/2)) for lam in [1.,2.])
 return np.array(rows)
rows=[];start=time.monotonic()
for n in [16,32,64]:
 data=setup(n,32);raw=calc(np.zeros(4),data);sol=root(lambda a:calc(a,data),np.zeros(4),jac=lambda a:jac(a,data),tol=1e-10)
 fitted=calc(sol.x,data);del data
 validation=setup(n,64);x,h,mass,W=validation;f=h+sol.x@W
 checked=calc(sol.x,validation,[.5,1,1.5,2,2.5,3,4]);chars=[float(np.exp(-xi*xi/2)*(mass@(np.cos(xi*x)*np.expm1(xi*xi*(f-1)/2)))) for xi in [.5,1,2,3]]
 row=dict(n=n,success=bool(sol.success),message=sol.message,nfev=int(sol.nfev),coefficients=sol.x.tolist(),raw_residuals=raw.tolist(),fitted_residuals=fitted.tolist(),validation_64=checked.tolist(),variance_f=float(mass@((f-1)**2)),sampled_f_bounds=[float(min(f)),float(max(f))],characteristic_residuals=chars)
 rows.append(row);print(json.dumps(row),flush=True);del validation,x,h,mass,W,f
out=dict(parameters=dict(n=[16,32,64],R=R,fit_order=32,validation_order=64,centers=centers.tolist(),width=width,pin_radius='n^-3',tau='16*(s-.5)^2*(s-1.5)^2',validation_parameters=[.5,1,1.5,2,2.5,3,4],characteristic_parameters=[.5,1,2,3]),epsilon=EPS,limit_variance=EPS/2,rows=rows,runtime=dict(python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__,pid=os.getpid(),elapsed_seconds=time.monotonic()-start),script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),scope='floating smooth microstructure on |x|<=14; quadrature/truncation not certified; exact finite-constraint correction is proved separately')
Path(sys.argv[1]).write_text(json.dumps(out,indent=2)+'\n')
