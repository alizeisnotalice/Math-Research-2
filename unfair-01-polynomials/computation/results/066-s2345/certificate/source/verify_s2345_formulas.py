#!/usr/bin/env python3
"""Check both Newton systems and all projection identities with SymPy."""

import sympy as sp

L,a,B,c,d,e,f=sp.symbols("L a B c d e f")
b=L*a+B
P=[1,L,a,b,c,d,e,f,1]


def moments(coefficients):
    S=[8]
    for k in range(1,6):
        S.append(sp.expand(-sum(coefficients[j]*S[k-j] for j in range(1,k))-k*coefficients[k]))
    return S[1:]


forward=[-L,L**2-2*a,-L**3-3*B,L**4+2*a*a+4*L*B-4*c,
         -L**5-5*L*L*B+5*B*a+5*L*c-5*d]
reverse=[-f,f*f-2*e,-f**3+3*f*e-3*d,f**4-4*f*f*e+2*e*e+4*f*d-4*c,
         -f**5+5*f**3*e-5*f*e*e+5*(e-f*f)*d+5*f*c-5*b]
for coefficients,expected in ((P,forward),(P[::-1],reverse)):
    assert moments(coefficients)==list(map(sp.expand,expected))
x,r,s,t,j=sp.symbols("x r s t j")
p,rho,z=sp.symbols("p rho z")
G=sp.symbols("g0:9")
H=sp.symbols("h")
cofactor=sum(g*x**k for k,g in enumerate(G))
irr=r*x*x+L*r*x**3+s*x**4+t*x**5
lo=sp.expand(irr*cofactor)
assert lo.coeff(x,3).subs(G[0],1)==L*r+r*G[1]
subs={G[0]:1,G[1]:-L,G[2]:L*L-p,G[3]:-B-L*(L*L-p),G[4]:rho-z}
for k,expected in ((4,s-r*p),(5,t-L*s-r*B),
                   (6,r*(rho-z+L*subs[G[3]])+s*subs[G[2]]-L*t)):
    assert sp.expand(lo.coeff(x,k).subs(subs)-expected)==0
hi=sp.expand(x**8*irr.subs(x,1/x)*cofactor)
assert hi.coeff(x,4).subs(G[0],1)==s+t*G[1]
assert hi.coeff(x,5).subs({G[0]:1,s:t*j})==L*r+t*j*G[1]+t*G[2]
assert sp.expand(hi.coeff(x,6).subs({G[0]:1,s:t*j})-(r*(1+L*G[1])+t*(j*G[2]+G[3])))==0
assert sp.expand(hi.coeff(x,6).subs({G[0]:1,s:t*j,G[3]:-H-e*G[1]-f*G[2]}).subs({f:j,e:L*r/t})-(r-t*H))==0
print(f"SymPy {sp.__version__}: both Newton systems and both endpoint projections match.")
