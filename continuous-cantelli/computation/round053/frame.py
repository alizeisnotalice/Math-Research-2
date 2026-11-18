"""Deterministic finite-section check of the triangular Fourier frame."""
import numpy as np
import json
from fractions import Fraction
from math import pi
centers=np.cumsum(1+0.4*np.sin(np.arange(250))**2)
T=8
matrix=T*np.sinc(T*(centers[:,None]-centers[None,:])/(2*pi))**2
vals=np.linalg.eigvalsh(matrix)
bound=Fraction(40,3*T)
assert vals[0]>T-float(bound) and vals[-1]<T+float(bound)
print(json.dumps(dict(numpy=np.__version__,min_spacing=1,half_length=T,row_bound=str(bound),smallest=vals[0],largest=vals[-1]),indent=2))
