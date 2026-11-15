from fractions import Fraction as Q
import sympy as S
import json
s,t=S.symbols('s t')
r=3;f=[S.Integer(1)]
for n in range(1,8):f.append(S.expand(1-s*f[-1]-(t*f[n-r] if n>=r else 0)))
p=S.Poly(20480*t*t+17728*t-8425,t)
lo,hi=Q(218,640),Q(219,640)
assert 495**2<245229<496**2
assert S.simplify(p.as_expr().subs(t,(S.sqrt(245229)-277)/640))==0
out=[]
for n in range(7):
 poly=S.Poly(f[n].subs(s,S.Rational(3,4)),t)
 lower=upper=Q(0)
 for (power,),a in poly.terms():
  a=Q(int(a.p),int(a.q));x,y=lo**power,hi**power
  lower+=a*(x if a>=0 else y);upper+=a*(y if a>=0 else x)
 assert lower>0 and upper<=1
 out.append({'index':n,'expression':str(poly.as_expr()),'lower':str(lower),'upper':str(upper)})
assert S.rem(S.Poly(f[7].subs(s,S.Rational(3,4)),t),p).is_zero
# Validate the general coefficient formula on a symbolic grid; the identity
# for all n follows separately from formal geometric/binomial expansion.
for rr in range(2,6):
 fs=[S.Integer(1)]
 for n in range(1,16):
  fs.append(S.expand(1-s*fs[-1]-(t*fs[n-rr] if n>=rr else 0)))
  formula=sum((-t)**q*sum(S.binomial(i+q,q)*(-s)**i for i in range(n-q*rr+1)) for q in range(n//rr+1))
  assert S.expand(fs[n]-formula)==0
print(json.dumps({'r':3,'s':'3/4','t':'(sqrt(245229)-277)/640','minimal_polynomial':str(p.as_expr()),'t_interval':[str(lo),str(hi)],'positive_unit_prefix':out,'F7_exactly_zero':True,'scope':'Local binary prefix through degree7 only; not a finite nonnegative factorization','general_formula_grid':{'r':[2,5],'n':[1,15]}},indent=2))
