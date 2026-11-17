"""Second-order compatibility: symbolic formulas and scaled floating verification."""
import os
for key in ['OPENBLAS_NUM_THREADS','OMP_NUM_THREADS']:os.environ[key]='1'
from pathlib import Path
import sys,json,hashlib,platform,math
import sympy as S
import numpy as np,scipy
from scipy.integrate import quad
from scipy.optimize import brentq
c,d4=S.symbols('c d4',positive=True)
mu=2/c;var=2/c
m2=mu**2+var;m4=mu**4+6*mu**2*var+3*var**2
s2=S.simplify(-m2/2-d4*m4/48)
expected=-1/c-2/c**2-d4*(1/(4*c**2)+1/c**3+1/(3*c**4))
assert S.simplify(s2-expected)==0
z=S.symbols('z');assert S.expand(S.series((1+z)**-3,z,0,3).removeO()).coeff(z,2)==6
assert S.rf(3,4)/(2**2*S.factorial(2))==45
assert S.rf(5,2)/2==15
src=Path(__file__).with_name('round007-matched-results.json');p=json.loads(src.read_text())['parameters']
k,u,A,B,cc=[p[key] for key in ['kappa','u','A','B','curvature']];g=1+k;alpha=1-A+B
fourth=-3*A*u*u+B*(45+18*k+3*k*k)
ss=float(s2.subs({c:cc,d4:fourth}));dd=6*g-1.5*g*g;required_q=(ss-6*g)/g**2
Csur=math.sqrt(2/cc)*math.exp(1/cc);Ctail=B*g**2.5/2
root=brentq(lambda x:math.log(A/B)+1.5*math.log1p(x*x)-(u-k)*x*x/2,.0001,5)
rows=[]
for lam in [20.,40.,80.,160.,320.]:
 def sp(y):
  x=y/lam;delta=-A*math.expm1(-u*x*x/2)+B*math.expm1(-1.5*math.log1p(x*x)-k*x*x/2)
  ff=alpha+delta
  return math.exp(y-lam*lam*delta/2-x*x/2)/math.sqrt(2*math.pi)*(-math.expm1(-lam*lam*(1-ff)/2))
 cuts=sorted(set([-root*lam,root*lam]+[y for y in [-12.,0.,12.] if abs(y)<root*lam]))
 sv=se=0.
 for lo,hi in zip(cuts,cuts[1:]):
  val,err=quad(sp,lo,hi,epsabs=1e-12,epsrel=1e-12,limit=250);sv+=val;se+=err
 # Exact Gaussian completion for the positive part, centered at lambda/g.
 def dp(y):
  x=lam/g+y/math.sqrt(g)
  if abs(x)<root:return 0.
  ratio=math.exp(math.log(A/B)+1.5*math.log1p(x*x)-(u-k)*x*x/2)
  h=B*math.exp(-k*x*x/2)*(1+x*x)**-1.5*(1-ratio)
  zz=lam*lam*h/2;corr=1. if zz==0 else -math.expm1(-zz)/zz
  return math.exp(-y*y/2)/math.sqrt(2*math.pi)*(1+x*x)**-1.5*(1-ratio)*corr
 # Outside |y|<=14 omitted; floating check, not interval-certified.
 dv,de=quad(dp,-14,14,epsabs=1e-14,epsrel=1e-12,limit=300)
 dn=dv*lam**3*B/(2*math.sqrt(g))/Ctail
 sn=sv/Csur
 rows.append(dict(lambda_value=lam,surplus_normalized=sn,deficit_normalized=dn,surplus_second=(sn-1)*lam*lam,deficit_second=(dn-1)*lam*lam,quadrature_errors=[se,de]))
out=dict(formula=str(s2),parameters=p,derived=dict(fourth_derivative=fourth,surplus_second=ss,deficit_second=dd,required_tail_q=required_q,actual_tail_q=-1.5,Csur=Csur,Ctail=Ctail),rows=rows,runtime=dict(python=platform.python_version(),sympy=S.__version__,scipy=scipy.__version__,numpy=np.__version__,pid=os.getpid()),source_sha256=hashlib.sha256(src.read_bytes()).hexdigest(),script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),scope='exact symbolic formula; numerical profile from round7, floating quadrature and truncation not certified')
Path(sys.argv[1]).write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
