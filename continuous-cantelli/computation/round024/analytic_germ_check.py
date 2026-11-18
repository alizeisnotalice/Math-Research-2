"""55-digit calibration for v=u⁴(1+u/5), |u|<1, ρ(x)=p₁(x), x₀=3."""
import json
import platform
import mpmath as mp

mp.mp.dps = 55
theta = mp.pi / 6
b = mp.sin(theta)
d = mp.cos(theta)
B = 3*d/4
x0 = mp.mpf(3)
eta = mp.mpf(1)/5


def rho(z):
    return mp.exp(-z*z/2)/mp.sqrt(2*mp.pi)


def saddles(lam):
    a = 1/mp.sqrt(2*lam)
    def phase(z):
        return 1j*z-z**4*(1+eta*a*z)/4
    def derivative(z):
        return 1j-z**3-5*eta*a*z**4/4
    def curvature(z):
        return 3*z**2+5*eta*a*z**3
    z = [mp.findroot(derivative, s*d+1j*b) for s in (1,-1)]
    phases = [phase(t) for t in z]
    amplitudes = [rho(x0+a*t)*mp.sqrt(2*mp.pi)/mp.sqrt(curvature(t)) for t in z]
    offset = lam*mp.im(phases[0]-phases[1])+mp.arg(amplitudes[0])-mp.arg(amplitudes[1])
    return a, phase, phases, amplitudes, offset


def run(k):
    lam = mp.findroot(lambda l: saddles(l)[-1]-2*mp.pi*k, mp.pi*k/B)
    a, F, S, C, _ = saddles(lam)
    best = max(mp.re(s) for s in S)
    def integrand(z):
        return rho(x0+a*z)*mp.exp(lam*(F(z)-best))
    width = 1/mp.sqrt(lam)
    R = 1/a
    xcuts = {-R, -2, 0, 2, R}
    for sign in (-1,1):
        xcuts.update(sign*d+j*width for j in (-3,0,3))
    cuts = [-R] + [x+1j*b for x in sorted(xcuts)] + [R]
    observed = mp.quad(integrand, cuts)
    predicted = sum(Cj*mp.exp(lam*(sj-best)) for Cj,sj in zip(C,S))/mp.sqrt(lam)
    return dict(k=k, a=mp.nstr(a,13), xi=mp.nstr(lam/a,13),
                real_ratio=mp.nstr(mp.re(observed/predicted),16),
                imag_ratio=mp.nstr(mp.im(observed/predicted),16),
                log_amplitude_difference=mp.nstr(mp.log(abs(C[0]/C[1]))+
                                                  lam*mp.re(S[0]-S[1]),11))


print(json.dumps(dict(python=platform.python_version(), mpmath=mp.__version__,
                      dps=mp.mp.dps, n=4, eta="1/5", halfwidth=1, x0=3,
                      rows=[run(k) for k in (6, 24, 96, 384)]), indent=2))
