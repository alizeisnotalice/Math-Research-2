"""Direct finite-window complex quadrature for the critical, normalized linear response."""
import json
import platform
import numpy as np
import scipy
from scipy.integrate import quad_vec

a = 0.0002
B = 32*a
z = 1 + 1j*a
g0 = 1/np.sqrt(2*np.pi)
limit = -(z*z/2)*g0*np.exp(-z*z)/2 / np.sqrt(1+2*B*(a-1j))

def pair(w):
    return [float(w.real), float(w.imag)]

rows = []
for n in (20000, 80000, 320000):
    d = 1/(2*np.pi*n)
    L = np.sqrt(B/d)
    def integrand(w):
        t = L*w
        factor = 1+d*t
        u = d*factor
        T = 1+factor**4
        A = L**2+T
        alpha = -T/A
        log_derivative = -1/(2*A)+t**2/(2*A**2)-1j*z*t*L**2/A**2-z*z*L**4/(2*A**2)
        response = 0.5*(np.exp(1j*d*t*t/factor)+np.exp(1j*(t+t/factor)))
        gaussian = L*g0*np.exp(-t*t/(2*A))/np.sqrt(A)
        correction = np.exp(-a*d*t*t/factor + 1j*z*t*alpha-z*z*L**2*T/(2*A))
        return g0*np.exp(-u*u/2)*factor**4*gaussian*correction*log_derivative*response
    result, err = quad_vec(integrand, -10, 10, epsabs=1e-11, epsrel=1e-11, limit=10000)
    rows.append(dict(n=n, L=float(L), response=pair(result), absolute_error=float(abs(result-limit)),
                     estimated_quadrature_error=float(err)))
assert rows[-1]["absolute_error"] < rows[0]["absolute_error"] / 8
print(json.dumps(dict(python=platform.python_version(), scipy=scipy.__version__, numpy=np.__version__,
                      a=a, B=B, q=4, p=1, omega=1, kappa=1, theta=0, H="cosine", randomness="none",
                      quantity="linear response only divided by eta, in complex modulation",
                      window=[-10, 10], epsabs=1e-11, epsrel=1e-11,
                      predicted_limit=pair(limit), rows=rows), indent=2))
