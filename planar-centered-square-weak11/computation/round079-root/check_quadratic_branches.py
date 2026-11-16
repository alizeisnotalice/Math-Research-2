from pathlib import Path
import sympy as s,json,os,platform,hashlib,time
start=time.monotonic()
n=s.symbols('n',positive=True);X,Y,c=s.symbols('X Y c',real=True);r=s.symbols('r',real=True)
h=(n+1)/n**2;eta=n**2/(n+1)**2;e=1/n**2
F=(4-eta)*r**2+eta*(X+Y-2*h)*r-eta*(X-h)*(Y-h)
assert s.simplify(F-(4*r**2-eta*(r-X+h)*(r-Y+h)))==0
err=s.factor(eta*h*(1-2*e+h));assert s.simplify(err-(n*n+n-1)/(n*n*(n+1)))==0
G=3*r*r+(X+Y)*r-X*Y
rows=[]
for var in (X,Y):
 for sign in (1,-1):
  f=s.expand(G.subs(r,sign*var+c));co=s.expand(f).coeff(var,2)
  assert co==(4 if sign==1 else 2)
  rows.append(dict(variable=str(var),sign=sign,polynomial=str(f),quadratic_coefficient=str(co)))
result=dict(scope='Symbolic identities only; no measure convergence, area certification, parameter search or C2 theorem',python=platform.python_version(),sympy=s.__version__,pid=os.getpid(),script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),finite_n_error=str(err),limit_distance_branches=rows,seconds=time.monotonic()-start)
Path(__file__).with_name('quadratic-branches.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result))
