"""Level-by-level annular transform with Legendre integration in scaled gap distance."""
import json
import math
import platform
import numpy as np
import scipy
from scipy.special import roots_legendre

def harmonics(rho, alpha, max_n, order):
    n = np.arange(1, max_n+1, dtype=float)
    sig = n**(-2/alpha)
    phase = 2*math.pi*n
    centers = (-1.)**n
    I = np.zeros(n.size)
    Z = np.zeros(n.size)
    nodes, weights = roots_legendre(order)
    nodes = (1+nodes)/2
    weights = weights/2
    for level in range(175):
        H = (1-2*rho)*rho**level/2
        upper = np.minimum(H/sig, 70**(1/alpha))
        u = upper[:,None]*nodes
        mw = np.exp(-u**alpha)*weights
        factor = 2**(level+1)*upper*sig
        Z += factor*mw.sum(axis=1)
        I += factor*centers*(np.cos(phase[:,None]*(H-sig[:,None]*u))*mw).sum(axis=1)
        centers *= np.cos(math.pi*n*(1-rho)*rho**level)
    return I, Z, centers

rho = 0.4
dim = math.log(2)/math.log(1/rho)
rows = []
errors = []
for alpha in (1.4,1.8):
    I,Z,mu = harmonics(rho,alpha,1024,76)
    I2,Z2,_ = harmonics(rho,alpha,1024,108)
    errors.append(float(np.max(np.abs(I/Z-I2/Z2))))
    for low in (32,64,128,256):
        sel = slice(low-1,4*low)
        scale = low**(1-dim)
        per_freq_sig = np.arange(1,1025,dtype=float)**(-2/alpha)
        rows.append(dict(alpha=alpha, low=low, high=4*low,
                         natural_energy=float(np.sum(mu[sel]**2)/scale),
                         tilted_energy=float(np.sum((I2[sel]/Z2[sel])**2)/scale),
                         weighted_raw_energy=float(np.sum((I2[sel]/per_freq_sig[sel]**(1-dim))**2)/scale)))
assert max(errors) < 1e-5
print(json.dumps(dict(python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__,
                      rho=rho, dimension=dim, cutoff_exponent=70, levels=175, orders=[76,108],
                      normalized_max_quadrature_discrepancies=errors,
                      deterministic=True, rows=rows),indent=2))
