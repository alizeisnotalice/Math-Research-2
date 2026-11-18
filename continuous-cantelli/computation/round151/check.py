import sympy as s,json
x,t,b,v=s.symbols('x t b v',positive=True)
d=x*x*(v-1)+2*b*x-4*s.log(x)
reconstructed=2*s.log(t)-b*t+b*t*t/x-2*(t/x)**2*s.log(x)-(t/x)**2*d/2
assert s.simplify(reconstructed-(2*s.log(t)-b*t+t*t*(1-v)/2))==0
print(json.dumps(dict(sympy=s.__version__,tail_exponent_residual='0')))
