"""Finite independent calibration of the first two even feedback coordinates.

The proof of existence to all orders uses compactness and is independent of this solver.
"""
import json
import math
import platform
import mpmath as mp
from fractions import Fraction

mp.mp.dps=65
R=mp.mpf('.4')
g0=1/mp.sqrt(2*mp.pi)
e2,e4=mp.mpf('.08'),mp.mpf('.04')

def bump(u,e):
    q=abs(u)/e
    if q<=mp.mpf('.5'): return mp.mpf(1)
    if q>=1: return mp.mpf(0)
    t=2*q-1
    a=mp.exp(-1/t)
    b=mp.exp(-1/(1-t))
    return b/(a+b)

def he(k,z):
    x,y=mp.mpf(1),z
    for j in range(1,k):
        x,y=y,z*y-j*x
    return y if k else 1

def moments(a2,a4):
    def I(n):
        def integrand(u):
            if u==0: return mp.mpf(0)
            w=1+a2*u*u*bump(u,e2)+a4*u**4*bump(u,e4)
            assert w>mp.mpf('.5')
            v=u**4*w
            z=u/mp.sqrt(v)
            return (2*g0**2*mp.exp(-(u*u+z*z)/2)*he(n-2,z)
                    /(mp.factorial(n)*v**(mp.mpf(n-1)/2)))
        return mp.quad(integrand,[0,e4/2,e4,e2,R])
    return [I(6),I(8)]

def desired(I6,I8):
    a2=mp.mpf('.5')-14-2*I6/g0
    a4=a2/2-mp.mpf('.125')-mp.mpf('22.5')*(2*a2-mp.mpf('.5'))-495-2*I8/g0
    return a2,a4

a2=a4=mp.mpf(0)
for n in range(3):
    a2,a4=desired(*moments(a2,a4))
I6,I8=moments(a2,a4)
B6=g0*((a2-mp.mpf('.5'))/2+7)
B8=g0*((a4-a2/2+mp.mpf('.125'))/2+mp.mpf('11.25')*(2*a2-mp.mpf('.5'))+mp.mpf('247.5'))
residuals=[abs(I6+B6),abs(I8+B8)]
assert max(residuals)<mp.mpf('1e-25')
for n,j,pred in [(6,2,Fraction(7)),(8,2,Fraction(45,4)),(8,3,Fraction(495,2))]:
    actual=Fraction(math.factorial(n+2*j-2),2**j*math.factorial(j)*math.factorial(n))
    assert actual==pred
print(json.dumps(dict(python=platform.python_version(),mpmath=mp.__version__,decimal_places=mp.mp.dps,
    model="even two-coordinate bump example",r=str(R),eps=[str(e2),str(e4)],
    a2=mp.nstr(a2,32),a4=mp.nstr(a4,32),
    exterior=[mp.nstr(I6,32),mp.nstr(I8,32)],
    diagonal=[mp.nstr(B6,32),mp.nstr(B8,32)],
    residual=[mp.nstr(x,8) for x in residuals],iters=3,
    note="numerical finite two-coordinate calibration only"),indent=2))
