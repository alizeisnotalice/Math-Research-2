"""Numerical calibration of Fourier saddle normalization via shifted contours."""
import json
import platform
import mpmath as mp

mp.mp.dps = 55
x0 = mp.mpf(3)
c = mp.mpf(1)
rows = []

def gaussian(z):
    return mp.exp(-z*z/2)/mp.sqrt(2*mp.pi)


for n in (4, 6):
    m = n//2
    theta = mp.pi/(2*(n-1))
    b = mp.sin(theta)
    q = mp.cos(theta)
    decay = (n-1)*b/n
    oscillation = (n-1)*q/n
    shift = (n-2)*theta/2
    for k in (4, 16, 64):
        lam = (2*mp.pi*k+shift)/oscillation
        xi = (lam**(n-1)*m*c)**(1/(mp.mpf(n-2)))
        a = lam/xi
        def integrand(x):
            z = x + 1j*b
            return gaussian(x0+a*z)*mp.exp(lam*(1j*z-z**n/n+decay))
        width = 1/mp.sqrt(lam)
        breaks = sorted(set([-2, 0, 2]+[s*q+j*width for s in (-1,1) for j in (-2,0,2)]))
        result = mp.quad(integrand, breaks) * mp.sqrt(lam)
        leading = 2*gaussian(x0)*mp.sqrt(2*mp.pi/(n-1))
        rows.append(dict(n=n, k=k, xi=mp.nstr(xi,15),
                         real_ratio=mp.nstr(mp.re(result/leading),16),
                         imag_ratio=mp.nstr(mp.im(result/leading),16)))
print(json.dumps(dict(python=platform.python_version(), mpmath=mp.__version__,
                      dps=mp.mp.dps, x0=str(x0), c=str(c), rows=rows),indent=2))
