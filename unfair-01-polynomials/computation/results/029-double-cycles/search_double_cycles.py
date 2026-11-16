import json,math,hashlib
from pathlib import Path
import sympy as sp
x=sp.symbols('x');H=1+x+x*x+x**7+x**8; C=sp.Poly(H.subs(x,x**3)*H.subs(x,x**5),x,domain=sp.ZZ)
def asc(poly):
 P=sp.Poly(poly,x);return [int(P.nth(i)) for i in range(P.degree()+1)]
def one(p):
 for ell in list(sp.primerange(2,2001)):
  if ell%p!=1:continue
  h=sp.Poly(H,x,modulus=ell); hp=sp.Poly(H.subs(x,x**p),x,modulus=ell); hq=sp.Poly(H.subs(x,x**(5 if p==3 else 3)),x,modulus=ell); cp=sp.Poly(C,x,modulus=ell)
  fh=sp.factor_list(h)[1];fp=sp.factor_list(hp)[1];fq=sp.factor_list(hq)[1];ds=[f.degree() for f,e in fp for _ in range(e)]
  if any(f.degree()==p for f,e in fh) or ds.count(p)!=2 or any(d%p==0 for d in ds if d!=p) or any(f.degree()%p==0 for f,e in fq) or sp.gcd(cp,cp.diff()).degree()!=0:continue
  rest=[d for d in ds if d!=p]+[f.degree() for f,e in fq for _ in range(e)];M=math.lcm(*(rest+[1]));assert M%p
  def rows(v):return [{'ascending':asc(f),'multiplicity':e,'degree':f.degree()} for f,e in v]
  return {'ell':ell,'p':p,'H':rows(fh),'H_xp':rows(fp),'H_xq':rows(fq),'C_squarefree':True,'M':M,'condition':'Frob^M leaves exactly two p-cycles from H(x^p)'}
 return None
out={'H_ascending':[1,1,1,0,0,0,0,1,1],'limit':2000,'p3':one(3),'p5':one(5),'sympy':sp.__version__};Path('double-cycle-results.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps({'p3':out['p3'] is not None,'p5':out['p5'] is not None,'sha':hashlib.sha256(Path('double-cycle-results.json').read_bytes()).hexdigest()}))
