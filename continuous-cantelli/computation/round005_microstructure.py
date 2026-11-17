"""Continuous oscillatory near-Gaussian examples; floating evaluations only.
The convergence and non-exactness are established analytically in proof 005.
"""
import os
for key in ['OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS']:os.environ[key]='1'
import numpy as np
from scipy.special import ndtr
from scipy.integrate import simpson
from numpy.polynomial.legendre import leggauss
from pathlib import Path
import json,sys,platform,time,hashlib
EPS=.6*np.exp(-16/3);R=14.;H=2.
phi=lambda x:np.exp(-x*x/2)/np.sqrt(2*np.pi)

def build(n):
 h=1/n;delta=1/(8*n*n);slope_bound=8*n**3;segments=[]
 for k in range(-int(R*n),int(R*n)):
  left=k/n;right=(k+1)/n;c=(left+right)/2
  rs=np.exp(np.log(EPS)+.5*np.log(8)-(c-2)**2/(2*.125)+c*c/2)
  rt=np.exp(np.log(EPS)+.5*np.log(8/5)-(c-2)**2/(2*.625)+c*c/2)
  rm=1-rs-rt
  assert rm>0
  lengths=np.array([delta,(1-4*delta)*rs,delta,(1-4*delta)*rm,delta,(1-4*delta)*rt,delta])*h/2
  endpoints=[(1.,.5),(.5,.5),(.5,1.),(1.,1.),(1.,1.5),(1.5,1.5),(1.5,1.)]
  x=left;half=[]
  for j,(length,(fl,fr)) in enumerate(zip(lengths,endpoints)):
   y=c if j==6 else x+length
   if y>x:half.append((x,y,fl,fr))
   x=y
  cell=half+[(2*c-b,2*c-a,fb,fa) for a,b,fa,fb in reversed(half)]
  # Pin f_n(0)=1/2 via min(g_n(x),1/2+8*n^3*|x|).
  for a,b,fa,fb in cell:
   wa=.5+slope_bound*abs(a);wb=.5+slope_bound*abs(b)
   da,db=fa-wa,fb-wb
   if da<=0 and db<=0:segments.append((a,b,fa,fb))
   elif da>=0 and db>=0:segments.append((a,b,wa,wb))
   else:
    cross=a+(b-a)*(-da)/(db-da);fc=fa+(fb-fa)*(cross-a)/(b-a)
    if da>0:segments.extend([(a,cross,wa,fc),(cross,b,fc,fb)])
    else:segments.extend([(a,cross,fa,fc),(cross,b,fc,wb)])
 return np.array(segments),slope_bound

def nodes(segments,order,shift=0.,cut=None):
 a,b,fa,fb=segments.T.copy();fa+=shift;fb+=shift
 if cut is not None:
  take=(fa<cut)|(fb<cut);a,b,fa,fb=[v[take] for v in [a,b,fa,fb]]
  cross=(fa<cut)!=(fb<cut)
  roots=np.zeros_like(a);roots[cross]=a[cross]+(b[cross]-a[cross])*(cut-fa[cross])/(fb[cross]-fa[cross])
  clipleft=cross&(fa>=cut);clipright=cross&(fb>=cut)
  a[clipleft]=roots[clipleft];fa[clipleft]=cut
  b[clipright]=roots[clipright];fb[clipright]=cut
 z,w=leggauss(order)
 x=(a[:,None]+b[:,None])/2+(b-a)[:,None]*z/2
 f=(fa[:,None]+fb[:,None])/2+(fb-fa)[:,None]*z/2
 weights=(b-a)[:,None]/2*w*phi(x)
 return x.ravel(),f.ravel(),weights.ravel()

def interval_prob(lo,hi):
 lo=np.asarray(lo);hi=np.asarray(hi)
 return np.where(lo>0,ndtr(-lo)-ndtr(-hi),ndtr(hi)-ndtr(lo))

ns=[4,8,16,32,64,128,256];params=np.array([.5,1.,2.,3.,4.]);zs=np.linspace(-10,10,801)
t=17/32;lo=41/15;hi=44/15
limit_q=float(interval_prob(lo/np.sqrt(t),hi/np.sqrt(t))-EPS*interval_prob((lo-2)/np.sqrt(t-.375),(hi-2)/np.sqrt(t-.375)))
rows=[];start=time.monotonic()
for n in ns:
 seg,M=build(n);reports=[]
 for order in [8,16]:
  x,f,w=nodes(seg,order);rawmean=float(w@f+2*ndtr(-R));shift=1-rawmean
  assert abs(shift)<.25
  fn=f+shift
  Fs=np.array([np.dot(w,np.exp(lam*x-lam*lam*fn/2)) for lam in params])
  variance=float(w@((fn-1)**2))
  xq,fq,wq=nodes(seg,order,shift,cut=t);var=t-fq
  source=float(np.dot(wq,interval_prob((lo-xq)/np.sqrt(var),(hi-xq)/np.sqrt(var))))
  q=float(interval_prob(lo/np.sqrt(t),hi/np.sqrt(t))-source)
  report=dict(order=order,raw_mean=rawmean,normalization_shift=shift,variance_f=variance,
              F=Fs.tolist(),max_F_residual=float(np.max(np.abs(Fs-1))),Q_interval=q)
  if order==16:
   variance_kernel=H-fn;density=np.empty(len(zs))
   for j in range(0,len(zs),8):
    zz=zs[j:j+8,None];kernel=np.exp(-(zz-x)**2/(2*variance_kernel))/np.sqrt(2*np.pi*variance_kernel)
    density[j:j+8]=kernel@w
   target=np.exp(-zs*zs/4)/np.sqrt(4*np.pi)
   report.update(tv_grid_approx=float(.5*simpson(np.abs(density-target),x=zs)),
                 max_density_error=float(np.max(np.abs(density-target))),density=density.tolist())
  reports.append(report)
 row=dict(n=n,structural_lipschitz_bound=M,segments=len(seg),quadrature=reports)
 rows.append(row)
 print(json.dumps(dict(n=n,shift=reports[-1]['normalization_shift'],variance_f=reports[-1]['variance_f'],
                       max_F_residual=reports[-1]['max_F_residual'],TV_grid=reports[-1]['tv_grid_approx'],
                       Q_interval=reports[-1]['Q_interval'])),flush=True)
out=dict(parameters=dict(n=ns,R=R,H=H,delta='1/(8*n^2)',symmetric_cells=True,pinned_value='.5',
                         real_parameters=params.tolist(),Q_time=t,Q_interval=[lo,hi],quadrature_orders=[8,16]),
         epsilon=EPS,limit_variance_f=EPS/2,limit_Q_interval=limit_q,z_grid=zs.tolist(),rows=rows,
         numerical_tail_bounds=dict(mean_error_bound=.5*2*ndtr(-R),
            real_F_tail_upper_at_4=float(np.exp(.375*16)*(ndtr(4-R)+ndtr(-4-R))),
            TV_outside_z10_upper=float(np.exp(-100/5.5)+np.exp(-25))),
         runtime=dict(python=platform.python_version(),numpy=np.__version__,pid=os.getpid(),elapsed_seconds=time.monotonic()-start),
         script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
         scope='floating quadrature: tail bound formulas analytic, displayed numbers not directed intervals; normalized means use numerical integration; exact normalized sequence is defined in proof')
Path(sys.argv[1]).write_text(json.dumps(out,indent=2)+'\n')
print('elapsed_seconds',out['runtime']['elapsed_seconds'],flush=True)
