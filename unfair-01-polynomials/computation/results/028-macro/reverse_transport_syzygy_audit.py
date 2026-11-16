#!/usr/bin/env python3
"""Exact symbolic audit of the reverse 2a -> b transport row identity."""
import sympy as sp

La,Lb,Ua,Ub,D,z=sp.symbols("La Lb Ua Ub D z")
V=La*Lb
R=Ua*Lb
S=La*Ub-D
T=z*D                 # z represents u^{-h}
W=Ua*Ub
lhs=sp.expand(R*S*W-S**2*T-W**2*V)
assert sp.expand(lhs + D*(R*W+z*S**2)) == 0
print("RSW-S^2T-W^2V = -D*(RW+z*S^2) verified exactly")
