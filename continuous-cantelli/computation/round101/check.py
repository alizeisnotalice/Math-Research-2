import sympy as s
x,b,d=s.symbols('x b d')
k,l,m,n,p,q=s.symbols('k l m n p q',nonzero=True)
c,r=s.symbols('c r')
J=(b-x**2/2)*(k+l*x+m*x*x/2+n*b+p*x*b+q*b*b)
J=J+d*J.diff(x,2)/2+d**2*J.diff(x,4)/8
at={x:r*d,b:d/2+c*d*d}
slope=s.solve(J.diff(x).subs(at).expand().coeff(d,1),r)[0]
curv=s.factor(s.solve(J.subs(at).expand().coeff(d,2).subs(r,slope),c)[0])
assert slope==-l/k and s.simplify(curv-(m/k-l**2/k**2)/2)==0
print(s.__version__,slope,curv)
