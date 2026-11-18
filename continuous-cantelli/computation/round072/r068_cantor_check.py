"""Quadrature by gap generations for the two-end-digit circle Cantor set."""
import json
import platform
import math
import scipy
from scipy.integrate import quad

def gap_integrals(H, delta, w, alpha):
    upper = min(H/delta, 80**(1/alpha))
    def weight(v):
        return math.exp(-v**alpha)
    opts = dict(epsabs=2e-13, epsrel=2e-13, limit=200)
    mass = quad(weight, 0, upper, **opts)[0]
    cosine = quad(weight, 0, upper, weight='cos', wvar=w, **opts)[0]
    sine = quad(weight, 0, upper, weight='sin', wvar=w, **opts)[0]
    return tuple(delta*t for t in (mass,cosine,sine))

def series(b, k, alpha, a):
    N = b**k
    delta = (a*N**2)**(-1/alpha)
    w = 2*math.pi*N*delta
    I = []
    Z = []
    phase_average = -1.
    last = k+104
    for j in range(1,last+1):
        H = (b-2)/(2*b**j)
        mass, cosine, sine = gap_integrals(H, delta, w, alpha)
        if j <= k:
            contribution = cosine
        else:
            angle = math.pi*(b-2)/b**(j-k)
            contribution = phase_average*(math.cos(angle)*cosine+math.sin(angle)*sine)
            phase_average *= math.cos(math.pi*(b-1)/b**(j-k))
        I.append(2**j*contribution)
        Z.append(2**j*mass)
    integral = math.fsum(I)
    normalization = math.fsum(Z)
    s = math.log(2,b)
    scale = N**(-2)*delta**(-1-s)
    return dict(b=b, k=k, alpha=alpha, a=a, coefficient=integral,
                mass=normalization, ratio=integral/normalization,
                principal_scale_ratio=integral/scale,
                generation_mass_tail=(2/b)**last)

a=2.
rows=[series(b,k,alpha,a) for b in (3,4) for alpha in (1.,2.,4.) for k in (4,8,12)]
for b in (3,4):
    exact = -math.prod(math.cos(math.pi*(b-1)/b**j) for j in range(1,40))
    tail = next(row for row in rows if row['b']==b and row['alpha']==1 and row['k']==12)
    assert abs(tail['ratio']-exact) < 0.001
    for row in rows:
        if row['b']==b and row['alpha']==4 and row['k']>=8:
            assert math.copysign(1,row['coefficient']) == (1 if b==3 else -1)
    prev=series(b,7,2.,a)
    curr=series(b,8,2.,a)
    inc=(b/2)**8*curr['coefficient']-(b/2)**7*prev['coefficient']
    expected=2**(1-8)*math.sqrt(math.pi)/(2*math.sqrt(a))*math.exp(-math.pi**2/a)
    assert abs(inc-expected)<1e-10
print(json.dumps(dict(python=platform.python_version(),scipy=scipy.__version__,
                      arithmetic='double-precision numerical consistency check',randomness='none',
                      gaussian_cutoff_exponent=80, absolute_tolerance=2e-13,relative_tolerance=2e-13,
                      additional_fine_generations=104, rows=rows),indent=2))
