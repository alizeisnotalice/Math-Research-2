"""Exact action checks and high-precision calibration for the all-order multiplier."""
from fractions import Fraction as Q
import json
import platform
import mpmath as mp
import sympy as sp
mp.mp.dps=95
x=sp.symbols('x')
A=-(3-2*x)*sp.log(1-x)+(2-x)*sp.log(1-x/2)-1-(1-x)*(1-2*x)/(2*x)*sp.log(1-2*x)
p=sp.series(A,x,0,26).removeO()
for j in range(2,26):
    expected=((j-3)*(1-sp.Rational(2**(j-1),j+1))+sp.Rational(1,2**(j-1)))/(j*(j-1))
    assert p.coeff(x,j)==expected
assert Q(49,18)**3<Q(84375,4096)
assert Q(11,4)**8<Q(3**12,125)
assert 2*Q(1,4)**2/(3*(1-2*Q(1,4)))==Q(1,12)
def integral_action(y):
    return -(3-2*y)*mp.log1p(-y)+(2-y)*mp.log1p(-y/2)-1-(1-y)*(1-2*y)/(2*y)*mp.log1p(-2*y)
def eval_order(m):
    d=mp.mpf(m-1); M=mp.mpf(m)
    lm=d/(2*M-1)*(4*M/(2*M-1))**(M/d)
    q=mp.mpf(2)/5 if m==2 else 1/(M*d)
    def gain(lam):
        return (q*(2*M*mp.log(M)-2*d*mp.log(d)+mp.log(2*q)-1)
                +(.5-d*q)*(1+mp.log(lam/(1-2*d*q)))-.5*(1+mp.log(lam)))
    g=(mp.log(mp.mpf(3)**12/125)-8)/10 if m==2 else integral_action(1/M)/d
    assert abs(gain(lm)-g)<mp.mpf('1e-75')
    lam=lm*mp.exp(g/(2*d*q))
    assert abs(gain(lam)-g/2)<mp.mpf('1e-75')
    actual=(M/(2*M-1))*(d/((2*M-1)*lam))**(d/M)
    assert actual<mp.mpf('.25')
    return dict(m=m,lambda_min=mp.nstr(lm,15),lambda_test=mp.nstr(lam,15),
                positivity_slack=mp.nstr(mp.mpf('.25')-actual,12),
                test_rate_gain=mp.nstr(g/2,16),order_cubed_gain=mp.nstr(M**3*g/2,14))
print(json.dumps(dict(python=platform.python_version(),sympy=sp.__version__,mpmath=mp.__version__,
                      precision=mp.mp.dps,exact_series_orders=[2,25],
                      rows=[eval_order(m) for m in [2,3,4,5,8,16,32,100,1000,10000,1000000]]),indent=2))
