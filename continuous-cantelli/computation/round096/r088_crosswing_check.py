"""Check the endpoint Laplace scales for transmission across a fixed gap 2b."""
import json
import platform
import mpmath as mp

mp.mp.dps=70
b=mp.mpf('0.8')
a=2*b*b
rows=[]
for h in map(mp.mpf,('0.2','0.1','0.05','0.02','0.01','0.005','0.001')):
    # Exact integrated derivative-kernel envelope after tilting by exp(a/h).
    # Integrating r^(-3/2) exp(-a/r) once and twice yields the two scales.
    def integrand(y):
        return mp.exp(-y)/mp.sqrt(1+h*y/a)
    def integrand2(y):
        return y*mp.exp(-y)/(1+h*y/a)**mp.mpf('1.5')
    v=mp.quad(integrand,[0,1,mp.inf])/a
    x=mp.quad(integrand2,[0,1,mp.inf])/(a*a)
    direct=mp.gammainc(mp.mpf('0.5'),a/h,mp.inf)*mp.exp(a/h)/mp.sqrt(a*h)
    rows.append(dict(h=str(h),scaled_first=mp.nstr(v,16),scaled_second=mp.nstr(x,16),incomplete_gamma_difference=mp.nstr(abs(v-direct),6)))
print(json.dumps(dict(python=platform.python_version(),mpmath=mp.__version__,decimal_precision=mp.mp.dps,b=str(b),scale_exponent=str(a),speed_coefficient_limit=mp.nstr(1/a,18),radius_coefficient_limit=mp.nstr(1/a**2,18),rows=rows),indent=2))
