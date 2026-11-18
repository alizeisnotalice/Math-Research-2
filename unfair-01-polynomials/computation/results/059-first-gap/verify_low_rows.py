"""Verify the sharp fourth-row bound and the degree-six obstruction."""
from pathlib import Path
import argparse
import json
import sympy as s

t,x=s.symbols('t x')
m=1-t+t*t
n=1-t*m
p=1-(1-t)*n
v=t*t*m
R=t*(1-t*m)+(t*(2*t-1)+m)**2/4
assert s.expand(4*(1-R)-(1-t)*(5*t**3+3*t*(1-t)+3))==0
F=2*t**3-4*t**2+5*t-1
H=t**6-4*t**5+6*t**4-7*t**3+4*t**2-2*t+1
assert s.expand((1-t)+(m-t*(1-t))*(1-(1-t)*m)-1-t*(t-1)*F)==0
assert s.Poly(F,t).is_irreducible
assert all(F.subs(t,w)!=0 for w in [-1,1,-s.Rational(1,2),s.Rational(1,2)])
assert s.expand(1-((1-t)*p+m*n)-v)==0
assert s.expand(t*v+p*m-1+H)==0
assert s.Poly(H,t,modulus=2).is_irreducible
assert s.Poly(H,t,modulus=2)==s.Poly(t**6+t**3+1,t,modulus=2)
assert [pow(2,i,9) for i in range(1,7)]==[2,4,8,7,5,1]
lo,hi=s.Rational(7,10),s.Rational(3,4)
assert H.subs(t,lo)>0>H.subs(t,hi)
assert s.Poly(H,t).count_roots(0,1)==1
A=1+t*x+n*x**3+p*x**4
B=1+(1-t)*x+m*x**2+v*x**5
C=s.Poly(s.expand(A*B),x)
reduced=[s.rem(C.nth(k),H,t) for k in range(10)]
assert reduced[:8]==[1,1,1,1,1,1,1,0]
assert s.expand(C.nth(8)-n*v)==0
out={'status':'verified_symbolic_identities_and_sharp_local_example','sympy':s.__version__,
     'no_first_gap_at_6_minpoly':str(F),'first_gap_at_7_minpoly':str(H),
     'p6_reduction_mod_2':str(s.Poly(H,t,modulus=2).as_expr()),
     'first_gap_at_7_root_interval':[str(lo),str(hi)],
     'root_endpoint_signs':[str(H.subs(t,q)) for q in (lo,hi)],
     'local_c0_through_c7':list(map(str,reduced[:8])),
     'local_c8':'n*t^2*m, strictly between 0 and 1'}

if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--out',type=Path,required=True)
    a=parser.parse_args()
    with a.out.open('x') as f:f.write(json.dumps(out,indent=2)+'\n')
    print(json.dumps(out))
