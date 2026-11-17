"""Formal jet recursion and a mean-preserving local repair; no exact Gaussian solution."""
import os
for key in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS']:os.environ[key]='1'
from pathlib import Path
import sympy as S
import numpy as np,scipy,json,sys,hashlib,platform,math
from scipy.integrate import quad
Y,z=S.symbols('Y z');jets={};records=[]
def moment(n):
 return sum(S.binomial(n,2*j)*2**(n-2*j)*S.factorial(2*j)/S.factorial(j) for j in range(n//2+1))
def expectation(poly):
 return sum(coef*moment(power[0]) for power,coef in S.Poly(S.expand(poly),Y).terms())
for k in range(1,5):
 d=S.symbols('d'+str(2*k+2));exponent=-z*Y**2/2
 for r in range(2,k+2):exponent-=(d if r==k+1 else jets[2*r])*Y**(2*r)*z**(r-1)/(2*S.factorial(2*r))
 coeff=S.expand(S.series(S.exp(exponent),z,0,k+1).removeO()).coeff(z,k)
 sk=expectation(coeff);target=S.rf(3,2*k)/S.factorial(k) # g=2
 value=S.solve(sk-target,d)[0];assert S.simplify(sk.subs(d,value)-target)==0
 jets[2*k+2]=value;records.append(dict(order=k,derivative=2*k+2,value=str(value),target=str(target)))
assert jets[4]==-S.Rational(180,19)
src=Path(__file__).with_name('round007-matched-results.json');p=json.loads(src.read_text())['parameters']
k,u,A,B,c=[p[key] for key in ['kappa','u','A','B','curvature']];g=1+k
old4=-3*A*u*u+B*(45+18*k+3*k*k);target=6*g-1.5*g*g
new4=-(target+1/c+2/c**2)/(1/(4*c*c)+1/c**3+1/(3*c**4));D=(new4-old4)/24;w=10000.
# Exact Gaussian moment cancellation for x^4*(1-(1+w)*x^2/5)*exp(-w*x^2/2).
ww=S.symbols('w',positive=True)
assert S.simplify(3/(1+ww)**S.Rational(5,2)-(1+ww)/5*15/(1+ww)**S.Rational(7,2))==0
def base(x):return 1-A*np.exp(-u*x*x/2)+B*(1+x*x)**-1.5*np.exp(-k*x*x/2)
def delta(x):return D*x**4*(1-(1+w)*x*x/5)*np.exp(-w*x*x/2)
rows=[];cuts=[-30,-2,-.25,-.1,-.03,0,.03,.1,.25,2,30]
for lam in [.5,1.,2.,4.,8.,12.]:
 vals=[]
 for repair in [False,True]:
  def integrand(x):return math.exp(lam*x-lam*lam*(base(x)+(delta(x) if repair else 0))/2-x*x/2)/math.sqrt(2*math.pi)
  val=sum(quad(integrand,a,b,epsabs=1e-12,epsrel=1e-12,limit=250)[0] for a,b in zip(cuts,cuts[1:]))
  vals.append(val)
 rows.append(dict(lambda_value=lam,before_F=vals[0],after_F=vals[1],change=vals[1]-vals[0]))
grid=np.linspace(-.2,.2,40001);err=float(np.max(np.abs(delta(grid))))
# Analytic upper bound on the uniform correction via max t^j exp(-t).
bound=abs(D)*(16/(w*w*math.e**2)+216*(1+w)/(5*w**3*math.e**3))
out=dict(formal_jets_alpha_half_c_one_q_zero=records,repair=dict(w=w,D=D,old_fourth=old4,new_fourth=new4,target_second=target,mean_change='exactly zero in the mathematical formula',sampled_sup_change=err,analytic_sup_bound_formula_value=bound),rows=rows,runtime=dict(python=platform.python_version(),sympy=S.__version__,scipy=scipy.__version__,numpy=np.__version__,pid=os.getpid()),input_sha256=hashlib.sha256(src.read_bytes()).hexdigest(),script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),scope='rational formal coefficients and exact moment cancellation; profile and quadrature are floating, sup-bound display not directed interval')
Path(sys.argv[1]).write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
