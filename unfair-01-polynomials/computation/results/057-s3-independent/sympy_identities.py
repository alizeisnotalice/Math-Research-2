"""Symbolic checks for the degree-eight S3/S5 and next S12 formulas."""
import json
import sympy as sp

x,a,L,b,c,d,e=sp.symbols('x a L b c d e')
P=[1,L,L**2,a,b,c,d,e,1]


def newton(poly):
    sums=[sp.Integer(8)]
    for k in range(1,9):
        sums.append(sp.expand(-sum(poly[i]*sums[k-i] for i in range(1,k))-k*poly[k]))
    return sums


s=newton(P)
assert sp.expand(s[2]+L**2)==0
assert sp.expand(s[3]-2*L**3+3*a)==0
for k,formula in {6:-6*d,7:-7*(e-L*d),8:8*L*e}.items():
    assert sp.expand(s[k]-s[k].subs({d:0,e:0})-formula)==0
G=[1,*sp.symbols('g1:8')]
G[2]=0
# c changes p5; integer low-degree cofactor correction is forced by the formal inverse.
changes={5:-c,6:c*(L-G[1]),7:c*L*G[1]}
Gp=[v+changes.get(i,0) for i,v in enumerate(G)]
A=sum(x**i*v for i,v in enumerate(G))*sum(x**i*v for i,v in enumerate([*P[:5],0,0,0,1]))
Ap=sum(x**i*v for i,v in enumerate(Gp))*sum(x**i*v for i,v in enumerate([*P[:6],0,0,1]))
difference=sp.Poly(sp.expand(Ap-A),x)
assert all(difference.nth(k)==0 for k in range(8))
assert sp.expand(difference.nth(8)-c*(G[3]+L**3-a))==0

# The forced sixth row of S12 fixes G5 and Q5 in terms of their fourth rows.
r,slope,a_bar,b2,b_bar,c3,D,rho,sigma,g1,q1,g2,q2,g3,q3=sp.symbols('r s abar b2 bbar c3 D rho sigma g1 q1 g2 q2 g3 q3',nonzero=True)
A5=a*(rho-D)+b2*g3+c3*g2+D*g1-slope*(rho-D)/r
B5=a_bar*(sigma-D)+b_bar*q3+c3*q2+D*q1-slope*(sigma-D)/r
assert sp.expand(sp.diff(A5-B5,D)-(g1-q1-a+a_bar))==0
assert sp.expand(sp.diff(A5-B5,e))==0
print(json.dumps({'sympy':sp.__version__,'status':'passed','verified':['S3 Newton-2,3,6,7,8 identities','S3 exact row8 translation','S12 row5 difference slope']}))
