"""Symbolic check of the crossed low/top coefficient recurrences at p=3,q=5."""
import sympy as s
x=s.symbols('x');p,q=3,5;a,b=13,12
f=s.symbols('f0:'+str(a+1));g=s.symbols('g0:'+str(b+1))
F=sum(f[i]*x**i for i in range(a+1));G=sum(g[i]*x**i for i in range(b+1))
A=s.expand(F.subs(x,x**p)*G.subs(x,x**q));D=p*a+q*b
B=s.expand(G.subs(x,x**p)*F.subs(x,x**q))
for n in range(10):
 lo=sum(f[n-q*t]*g[p*t] for t in range(n//q+1) if n-q*t<=a and p*t<=b)
 hi=sum(f[a-n+q*t]*g[b-p*t] for t in range(n//q+1) if a-n+q*t>=0 and b-p*t>=0)
 assert s.expand(A.coeff(x,p*n)-lo)==0
 assert s.expand(A.coeff(x,D-p*n)-hi)==0
assert s.expand(A.coeff(x,p*q)-(f[q]*g[0]+f[0]*g[p]))==0
assert s.expand(B.coeff(x,q*q)-(g[0]*f[q]+g[q]*f[q-p]))==0
# Equality edge m=n=2q for p=3,q=5: the q-2p term is absent.
fe=s.symbols('fe0:11');ge=s.symbols('ge0:11')
Fe=sum(fe[i]*x**i for i in range(11));Ge=sum(ge[i]*x**i for i in range(11))
Be=s.expand(Ge.subs(x,x**p)*Fe.subs(x,x**q))
assert s.expand(Be.coeff(x,q*q)-(ge[0]*fe[q]+ge[q]*fe[q-p]))==0
# The (3q,q) edge has disjoint A shifts for a binomial G because p does not
# divide q^2; this exact residue check is the algebraic isolation used above.
assert q*q % p != 0
# q>2p equality edge (p,q)=(3,7): B at q^2 has exactly the three
# s=0,1,2 terms, with the last one indexed by q-2p=1.
P,Q=3,7
ff=s.symbols('ff0:'+str(2*Q+1));gg=s.symbols('gg0:'+str(2*Q+1))
FF=sum(ff[i]*x**i for i in range(2*Q+1));GG=sum(gg[i]*x**i for i in range(2*Q+1))
BB=s.expand(GG.subs(x,x**P)*FF.subs(x,x**Q))
assert s.expand(BB.coeff(x,Q*Q)-(gg[0]*ff[Q]+gg[Q]*ff[Q-P]+gg[2*Q]*ff[Q-2*P]))==0
F0=1+x**(q-p)-x**q+x**(q+p)+x**(2*q)
G0=(1+x**p)*(1+x**q);H0=s.expand(F0*G0)
A0=s.Poly(s.expand(F0.subs(x,x**p)*G0.subs(x,x**q)),x)
B0=s.Poly(s.expand(G0.subs(x,x**p)*F0.subs(x,x**q)),x)
C0=s.Poly(s.expand(H0.subs(x,x**p)*H0.subs(x,x**q)),x)
assert min(A0.all_coeffs())>=0 and min(B0.all_coeffs())>=0 and set(H0.as_poly(x).all_coeffs())<={0,1}
assert C0.coeff_monomial(x**(p*q))==2
print({'p':p,'q':q,'checked_n':list(range(10)),
       'first_low_alias':str(A.coeff(x,p*q)),
       'first_top_alias':str(A.coeff(x,D-p*q)),
       'binary_4q_A_pq':str(A.coeff(x,p*q)),
       'binary_4q_B_q2':str(B.coeff(x,q*q)),
       'equality_edge_B_q2':str(Be.coeff(x,q*q)),
       'q_gt_2p_edge_B_q2':str(BB.coeff(x,Q*Q)),
       'sharp_example_crossed_nonnegative':True,'sharp_example_target_pq_coefficient':2})
