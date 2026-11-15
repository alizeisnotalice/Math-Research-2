"""Find a small exact 1=P13*f13+P15*f15+P17*f17 by rational linear algebra."""
import sympy as S
from itertools import combinations_with_replacement
p,q=S.symbols('p q')
b2=1-p; b9=p; b15=1-p; b17=1-q
b4=1-q-p*b2; b6=1-p*b4-q*b2; b8=1-p*b6-q*b4
b10=1-p*b8-q*b6; b11=1-b4-p*b9; b12=1-p*b10-q*b8
f=[S.expand(p*b11+q*b9+b6-1),S.expand(q*b11+b15+b8-1),S.expand(b10+p*b15+b17-1)]
def mons(d): return [p**i*q**j for i in range(d+1) for j in range(d+1-i)]
for d in range(0,11):
    M=mons(d); unknown=S.symbols('u:'+str(3*len(M)))
    identity=S.Poly(S.expand(sum(unknown[t*len(M)+h]*m*f[t] for t in range(3) for h,m in enumerate(M))-1),p,q)
    equations=[coef for _,coef in identity.terms()]
    A,b=S.linear_eq_to_matrix(equations,unknown)
    sol=S.linsolve((A,b),unknown)
    if sol is not S.EmptySet:
        sample=next(iter(sol))
        if sample != S.EmptySet:
            free=set().union(*(value.free_symbols for value in sample))
            chosen=[value.subs({symbol:0 for symbol in free}) for value in sample]
            Ps=[S.factor(sum(chosen[t*len(M)+h]*m for h,m in enumerate(M))) for t in range(3)]
            assert S.expand(sum(Ps[t]*f[t] for t in range(3))-1)==0
            print('degree',d,'unknowns',len(unknown),'equations',len(equations))
            print('f=',*[S.factor(x) for x in f],sep='\n')
            print('P=',*Ps,sep='\n')
            break
    print('no degree',d)
else: raise SystemExit('no identity through degree 10')
