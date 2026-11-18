"""Check reciprocal kernel exponent and dyadic leakage bound."""
import sympy as sp

x,t,f=sp.symbols('x t f', real=True,positive=True)
original=t*x-t**2*f/2-x**2/2
reciprocal=(1/f-1)*x**2/2-f*(t-x/f)**2/2
assert sp.expand(original-reciprocal)==0
# For centers in [R,2R], target in [-T,T] with R>=2T:
# distance>=R/2, variance<=2R/a, tail exponent>=aR/16.
R,a=sp.symbols('R a',positive=True)
assert sp.simplify((R/2)**2/(2*(2*R/a))-a*R/16)==0
print(sp.__version__, 'identity and tail bound verified')
