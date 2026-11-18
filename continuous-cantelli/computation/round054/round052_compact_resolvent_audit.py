"""Direct heat quadrature against the resolvent and exterior transform formulas."""
import json
import platform
import mpmath as mp

mp.mp.dps=40
V=mp.mpf(2)
q=mp.mpf('.8')
s=mp.sqrt(2*q)
g=lambda x: mp.exp(-x*x/2)/mp.sqrt(2*mp.pi)
v=lambda x: mp.mpf('.2')+mp.mpf('.3')*x*x
def verify(y):
    y=mp.mpf(y)
    sign=1 if y>1 else -1
    m=mp.quad(lambda x:g(x)*mp.exp(q*v(x)+sign*s*x),[0,1])
    direct=mp.quad(lambda x:g(x)*mp.exp(-q*(V-v(x))-s*abs(y-x))/s,[0,1])
    via_mgf=mp.exp(-q*V-sign*s*y)/s*m
    assert abs(direct-via_mgf)<mp.mpf('1e-35')
    return [str(y),mp.nstr(direct,24),mp.nstr(abs(direct-via_mgf),6)]
heat=mp.quad(lambda t:mp.exp(-q*t-1/(2*t))/mp.sqrt(2*mp.pi*t),[0,1,mp.inf])
assert abs(heat-mp.exp(-s)/s)<mp.mpf('1e-35')
print(json.dumps(dict(python=platform.python_version(),mpmath=mp.__version__,
   precision=mp.mp.dps,q=str(q),V=str(V),heat_resolvent_error=mp.nstr(abs(heat-mp.exp(-s)/s),6),
   one_sided_checks=[verify('-1.7'),verify('1.8')]),indent=2))
