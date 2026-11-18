"""Gauss-Jacobi nodes for q^(k-1/2)(1-q)^(-k) dq on [0,1]."""
import numpy as np
from scipy.special import roots_jacobi

ORDER=64
with open("nodes.bin","wb") as f:
    for k in (.25,.5):
        z,w=roots_jacobi(ORDER,-k,k-.5)
        np.column_stack(((z+1)/2,w/np.sqrt(2))).astype("<f4").tofile(f)
