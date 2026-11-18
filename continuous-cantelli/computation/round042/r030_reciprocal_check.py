"""Exact reciprocal kernel algebra and numerical saturated center-versus-delivered asymptotic."""
import json
import platform
import mpmath as m
import sympy as sp

x,t,f=sp.symbols('x t f', positive=True)
assert sp.expand((1/f-1)*x**2/2-(t-x/f)**2*f/2-(-x**2/2+t*x-t**2*f/2))==0
m.mp.dps=60
d=m.log(2)/2
def D(ell):
    return ell+m.log(ell)/2+d
def row(L):
    L=m.mpf(L)
    ell=m.findroot(lambda l: l+m.log(2*D(l))-L, L-m.log(2*L))
    D0=D(ell)
    u=m.exp(-ell)
    raw=m.exp(D0-ell-u*u/2)*m.sqrt(2*D0)/(2*(D0+1+1/(2*ell)))
    def source_y(y):
        l=L-m.log(y)
        return m.exp(y-y*y/(4*D(l))-L-(y*m.exp(-L))**2/2)/m.sqrt(2*m.pi)
    radius=2*m.sqrt(D0)
    middle=2*D0
    points=sorted(set([L/4,8*L]+[middle+j*radius for j in range(-16,17) if L/4<middle+j*radius<8*L]))
    delivered=m.quad(source_y,points)
    return {'log_t':int(L),'raw_center_density':m.nstr(raw,20),
        'delivered_core_window':m.nstr(delivered,20),'difference_times_log_t':m.nstr((raw-delivered)*L,20)}
print(json.dumps({'python':platform.python_version(),'sympy':sp.__version__,'mpmath':m.__version__,
    'mpmath_digits':m.mp.dps,'source_window':'t*u in [log(t)/4, 8log(t)]',
    'j':0,'q':0,'predicted_limit':'-1/2','rows':[row(L) for L in (20,80,320,1280)]},indent=2))
