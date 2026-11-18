"""Shifted-contour calibration for pure inverse-power tails, including p<1."""
import json
import platform
import mpmath as mp
mp.mp.dps=50
def run(p,xi):
    p=mp.mpf(p); xi=mp.mpf(xi); d=1
    a=(p*d*xi/2)**(1/(p+1))
    lam=a*xi; delta=a/xi
    theta=mp.pi/(2*(p+1)); rot=mp.exp(1j*theta)
    phi=min(theta/2,mp.pi/8); endrot=mp.exp(1j*phi)
    def F(z):return 1j*z-z**(-p)/p-delta*z*z/2
    def Fp(z):return 1j+z**(-p-1)-delta*z
    z=mp.findroot(Fp,rot)
    S=F(z)
    C=rot*mp.sqrt(2*mp.pi)/mp.sqrt(-rot**2*(-(p+1)*z**(-p-2)-delta))
    def integrand(w):
        return mp.exp(lam*(F(w)-S)) if w else 0
    width=1/mp.sqrt(lam)
    local=[z+rot*t for t in sorted(set([mp.mpf('-.25'),mp.mpf('.25'),mp.mpf(0)]
                                +[sg*b*width for sg in (-1,1) for b in (1,3) if b*width < .25]))]
    path=[0,.4*rot]+local+[2*rot,6*rot,6*endrot,9*endrot]
    result=mp.quad(integrand,path)*mp.sqrt(lam)/C
    return dict(p=float(p),xi=int(xi),real_ratio=mp.nstr(mp.re(result),13),
                imag_ratio=mp.nstr(mp.im(result),13),real_action=mp.nstr(mp.re(S),13))
rows=[run(p,t) for p in ('.5','1','2') for t in (16,48,128)]
print(json.dumps(dict(python=platform.python_version(),mpmath=mp.__version__,
                      dps=mp.mp.dps,tail_d=1,scaled_upper_radius=9,rows=rows),indent=2))
