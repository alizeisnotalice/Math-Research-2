import sympy as s
t=s.symbols('t',real=True)
f=s.Function('f');x=s.symbols('x',real=True)
J=t*s.sin(t)/2
assert s.simplify(s.diff(J,t,2)+J-s.cos(t))==0
w=f(s.atan(x))*s.sqrt(1+x*x)
target=(f(s.atan(x))+s.Subs(s.diff(f(t),t,2),t,s.atan(x)))/(1+x*x)**s.Rational(3,2)
assert s.simplify(s.diff(w,x,2)-target)==0
print('J + Jpp = cos(t); wxx = (1+x^2)^(-3/2) * (F + Fpp)')
