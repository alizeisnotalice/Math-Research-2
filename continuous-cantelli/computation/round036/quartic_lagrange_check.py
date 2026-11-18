"""Independent full-quartic coefficient check by moment-tilted Laguerre integration."""
import json
import math
import platform

import scipy
from scipy.special import roots_genlaguerre


def gaussian_taylor(x, order):
    data = [math.exp(-x*x/2)/math.sqrt(2*math.pi)]
    prev = 0.
    for j in range(order):
        nxt = (-x*data[-1]-prev)/(j+1)
        prev = data[-1]
        data.append(nxt)
    return data


def coefficient(order, d, x0, quadrature_order):
    data0 = gaussian_taylor(x0, order)
    nodes,weights = roots_genlaguerre(quadrature_order, (order-1)/2)
    weights = weights / sum(weights)
    terms = []
    for node,weight in zip(nodes,weights):
        for sign in (-1,1):
            b = sign*math.sqrt(2*d*node)
            data1 = gaussian_taylor(x0-1/b, order)
            prefactor = 1.
            value = 0.
            for j in range(order+1):
                value += prefactor*(data0[j]+(-1)**j*data1[j])
                prefactor *= -(order-j)/((2*order-j)*b)
            terms.append(weight*value/2)
    return math.fsum(terms)


rows = []
for x0,d in ((0,1),(3,.4)):
    g0 = gaussian_taylor(x0, 2)
    for n in (16,32,64,128,256):
        obs = coefficient(n,d,x0,56)
        independent = coefficient(n,d,x0,32)
        pred = 2*g0[0]+g0[2]/(d*(2*n-1))
        rows.append(dict(x0=x0,d=d,order=n,ratio=obs,leading=2*g0[0],
                         two_term_prediction=pred,scaled_error=n*n*(obs-pred),
                         quadrature_difference=abs(independent-obs)))
print(json.dumps(dict(python=platform.python_version(),scipy=scipy.__version__,
                      degrees=[32,56],rows=rows),indent=2))
