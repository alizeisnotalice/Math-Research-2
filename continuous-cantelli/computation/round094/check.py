from sympy import symbols,simplify
b,h,z=symbols('b h z',positive=True)
exponent=h*z*z/2-2*b*z
assert simplify(exponent.subs(z,2*b/h)+2*b*b/h)==0
print('Laplace optimization: -2 b^2 / h')
