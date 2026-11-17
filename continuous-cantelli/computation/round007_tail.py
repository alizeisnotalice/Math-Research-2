"""Exact normalization and tail-matching obstruction for a smooth Gaussian bump family."""
import os
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS']:os.environ[k]='1'
import sympy as S
import numpy as np
from scipy.integrate import quad
import scipy,json,sys,platform,hashlib,math
from pathlib import Path
u,v=S.Rational(7),S.Rational(1);A=S.Rational(1);B=S.Rational(1,2)
assert S.simplify(-A/S.sqrt(1+u)+B/S.sqrt(1+v))==0
alpha=1-A+B
assert alpha==S.Rational(1,2) and v==alpha/(1-alpha)
c=A*u-B*v
assert c==S.Rational(13,2)
x=S.symbols('x',real=True)
k,lam_sym,beta=S.symbols('k lam beta',real=True)
assert S.simplify(-(x-lam_sym)**2/2-k*x*x/2+beta*x+(1+k)*(x-(lam_sym+beta)/(1+k))**2/2+k*lam_sym**2/(2*(1+k))-beta*lam_sym/(1+k)-beta**2/(2*(1+k)))==0
fourth_condition=S.simplify(1/S.sqrt(15)+S.Rational(1,4)/S.sqrt(3)-S.Rational(1,3)-3/(8*S.sqrt(2)))
# f=1-exp(-7*x^2/2)+exp(-x^2/2)/2; t=exp(-x^2/2) in [0,1].
t=S.symbols('t',nonnegative=True)
assert S.expand(1-t**7+t/2-alpha)==S.expand((1-t)*(S.Rational(1,2)+sum(t**j for j in range(1,7))))
# Strict global minimum at x=0 follows from this factorization.
root=math.sqrt(math.log(2)/3)
phi=lambda x:math.exp(-x*x/2)/math.sqrt(2*math.pi)
f=lambda x:1-math.exp(-3.5*x*x)+.5*math.exp(-.5*x*x)
rows=[]
for lam in [5.,10.,20.,40.,80.]:
 # Multiply surplus and deficit by exp(alpha*lambda^2/2) to avoid underflow.
 def sp(y):
  xx=y/lam;ff=f(xx)
  exponent=lam*xx-lam*lam*(ff-.5)/2
  return phi(xx)*math.exp(exponent)*(-math.expm1(-lam*lam*(1-ff)/2))/lam
 surplus,se=quad(sp,-lam*root,lam*root,epsabs=1e-12,epsrel=1e-11,limit=250)
 # Rewrite deficit using its linear Gaussian convolution and (1-exp(-z))/z.
 def dp(xx):
  h=.5*math.exp(-xx*xx/2)*(-math.expm1(-3*xx*xx+math.log(2)))
  z=lam*lam*h/2
  corr=1. if z==0 else -math.expm1(-z)/z
  return lam*lam/4/math.sqrt(2*math.pi)*math.exp(-(xx-lam/2)**2)*(-math.expm1(-3*xx*xx+math.log(2)))*corr
 # Split around the saddle to prevent a long interval quadrature missing it.
 cuts=sorted(set([root]+[z for z in [lam/2-10,lam/2,lam/2+10] if z>root]))
 deficit=de=0.
 for l,r in zip(cuts,cuts[1:]):
  q,e=quad(dp,l,r,epsabs=1e-12,epsrel=1e-11,limit=250);deficit+=q;de+=e
 q,e=quad(dp,cuts[-1],math.inf,epsabs=1e-12,epsrel=1e-11,limit=250);deficit+=q;de+=e
 q,e=quad(dp,-math.inf,-root,epsabs=1e-12,epsrel=1e-11,limit=250);deficit+=q;de+=e
 pred_s=math.sqrt(2/float(c))*math.exp(1/float(c))/lam
 pred_d=lam*lam/(4*math.sqrt(2))
 rows.append(dict(lambda_value=lam,scaled_surplus=surplus,scaled_deficit=deficit,surplus_asymptotic=pred_s,deficit_asymptotic=pred_d,surplus_ratio=surplus/pred_s,deficit_ratio=deficit/pred_d,quadrature_errors=[se,de]))
out=dict(exact=dict(u=str(u),v=str(v),A=str(A),B=str(B),alpha=str(alpha),curvature=str(c),mean_f='1',tail_kappa='1',tail_power='0; necessary power is -3',fourth_condition_variance_minus_twice_covariance=str(fourth_condition),positive_minimum_factorization=str(S.factor(1-t**7+t/2-alpha))),rows=rows,runtime=dict(python=platform.python_version(),sympy=S.__version__,scipy=scipy.__version__,numpy=np.__version__,pid=os.getpid()),scope='floating quadrature supporting analytic asymptotics; candidate is rigorously excluded, not a solution',script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
Path(sys.argv[1]).write_text(json.dumps(out,indent=2)+'\n')
for row in rows:print(row)
