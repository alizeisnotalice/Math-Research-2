"""Independent adaptive-quadrature checks of signed intermediate heat mass.
Floating evidence only. Inputs are exploratory graphs, not exact solutions.
"""
import os
for key in ['OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS']:os.environ[key]='1'
import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.special import ndtr
from scipy.integrate import quad
import json,sys,platform,hashlib
from pathlib import Path
phi=lambda x:np.exp(-x*x/2)/np.sqrt(2*np.pi)
normal=lambda x,v:np.exp(-x*x/(2*v))/np.sqrt(2*np.pi*v)
def cdfdiff(lo,hi):
 return ndtr(-lo)-ndtr(-hi) if lo>0 else ndtr(hi)-ndtr(lo)
def check(path):
 d=json.loads(Path(path).read_text());k=np.array(d['parameters']['knots']);v=np.array(d['best']['values']);L=k[-1]
 def pieces(t):
  out=[]
  for i in range(len(k)-1):
   a,b=k[i:i+2];fa,fb=v[i:i+2];slope=(fb-fa)/(b-a)
   if fa>=t and fb>=t:continue
   if fa<t and fb<t:lo,hi=a,b
   elif slope>0:lo,hi=a,a+(t-fa)/slope
   else:lo,hi=a+(t-fa)/slope,b
   if hi>lo:out.append((lo,hi,a,fa,slope))
  return out
 def samples(t,order):
  z,w=leggauss(order);ys=[];ws=[];vs=[]
  for lo,hi,a,fa,slope in pieces(t):
   y=(lo+hi)/2+(hi-lo)*z/2
   ys.extend(y);ws.extend(w*(hi-lo)/2*phi(y));vs.extend(t-(fa+slope*(y-a)))
  return np.array(ys),np.array(ws),np.array(vs)
 def qgrid(t,x,order):
  y,w,var=samples(t,order)
  if len(y):source=((normal(x[:,None]-y,var)+normal(x[:,None]+y,var))*w).sum(axis=1)
  else:source=np.zeros_like(x)
  q=normal(x,t)-source
  if t>1:
   sd=np.sqrt((t-1)/t);mu=x/t
   tail=ndtr((mu-L)/sd)+ndtr((-L-mu)/sd)
   q-=normal(x,t)*tail
  return q
 times=np.array([.65,.8,.95,1.05,1.2,1.4,1.7]);xs=np.linspace(0,5,101)
 matrices={o:np.array([qgrid(t,xs,o) for t in times]) for o in [48,96]}
 ii,jj=np.unravel_index(np.argmin(matrices[96]),matrices[96].shape);t=float(times[ii]);x=float(xs[jj])
 def qpoint(t,x):
  value=normal(x,t);error=0.
  for lo,hi,a,fa,slope in pieces(t):
   def integrand(y):
    var=t-(fa+slope*(y-a))
    if var<=0:return 0.
    return phi(y)*(normal(x-y,var)+normal(x+y,var))
   val,err=quad(integrand,lo,hi,epsabs=1e-12,epsrel=1e-10,limit=300)
   value-=val;error+=err
  if t>1:
   sd=np.sqrt((t-1)/t);mu=x/t
   value-=normal(x,t)*(ndtr((mu-L)/sd)+ndtr((-L-mu)/sd))
  return value,error
 point,error=qpoint(t,x);left,right=x-.025,x+.025
 def intervalmass(t,left,right):
  result=cdfdiff(left/np.sqrt(t),right/np.sqrt(t));error=0.
  def kernel(y,var):
   sd=np.sqrt(var)
   return phi(y)*(cdfdiff((left-y)/sd,(right-y)/sd)+cdfdiff((left+y)/sd,(right+y)/sd))
  for lo,hi,a,fa,slope in pieces(t):
   def integrand(y):
    var=t-(fa+slope*(y-a))
    return 0. if var<=0 else kernel(y,var)
   val,err=quad(integrand,lo,hi,epsabs=1e-12,epsrel=1e-10,limit=300)
   result-=val;error+=err
  if t>1:
   val,err=quad(lambda y:kernel(y,t-1),L,np.inf,epsabs=1e-14,epsrel=1e-10,limit=100)
   result-=val;error+=err
  return result,error
 mass,masserr=intervalmass(t,left,right)
 return dict(input_file=Path(path).name,input_sha256=hashlib.sha256(Path(path).read_bytes()).hexdigest(),
             t=t,x=x,q48=float(matrices[48][ii,jj]),q96=float(matrices[96][ii,jj]),
             q_adaptive=float(point),q_adaptive_error_estimate=float(error),
             interval=[left,right],interval_mass=float(mass),interval_error_estimate=float(masserr),
             times=times.tolist(),x_grid=xs.tolist(),q96_matrix=matrices[96].tolist())
results=[check(p) for p in sys.argv[2:]]
out=dict(checks=results,runtime=dict(python=platform.python_version(),pid=os.getpid()),
         script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
         scope='floating quadrature, no interval certification; graph candidates fail exact transform identity')
Path(sys.argv[1]).write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps([{k:v for k,v in r.items() if k not in ['times','x_grid','q96_matrix']} for r in results],indent=2))
