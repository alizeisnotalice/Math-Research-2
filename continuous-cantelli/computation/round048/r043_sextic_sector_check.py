"""Numerical audit of a nonanalytic angular-flat sextic far fold and phase locking."""
import cmath
import json
import math
import platform
import scipy
from scipy.optimize import brentq,newton
from scipy.special import gammaln
import numpy as np

rho=2/(3*math.sqrt(3))
def far_fold():
    e=.2
    def pieces(U):
        l=cmath.log(1/(e*U))
        A=cmath.exp(-l*l)
        H=1+A
        Hp=2*l*A/U
        Hpp=(4*l*l-2*l-2)*A/U**2
        F=U-U**3*H
        Fp=1-3*U**2*H-U**3*Hp
        Fpp=-6*U*H-6*U**2*Hp-U**3*Hpp
        return F,Fp,Fpp
    h=brentq(lambda U:pieces(U)[1].real,.4,.8)
    rh=pieces(h)[0].real
    def g(x):
        return cmath.exp(-(2+x)**2/2)/math.sqrt(2*math.pi)
    L=g(e*h).real/(2*math.sqrt(rh*(-pieces(h)[2].real/2)))
    rows=[]
    for N in [12, 32, 64,128]:
        M=16*N
        scaled_rad=1-2/N
        r=rh*scaled_rad
        U=brentq(lambda x:pieces(x)[0].real-r,h,1.1)
        obs=0j
        for j in range(M):
            angle=2*math.pi*j/M
            W=r*cmath.exp(1j*angle)
            for _ in range(12):
                f,p,_=pieces(U)
                du=(f-W)/p
                U-=du
                if abs(du)<2e-15:break
            assert abs(pieces(U)[0]-W)<2e-13
            obs+=(-g(e*U)/pieces(U)[1])*cmath.exp(-1j*N*angle)
        C=math.exp(gammaln(2*N+1)-2*gammaln(N+1)-N*math.log(4))
        ratio=obs/M/scaled_rad**N/(C*L)
        rows.append(dict(N=N,ratio=ratio.real,imag_err=ratio.imag))
    return dict(e=e,fold_U=h,fold_W=rh,L=L,rows=rows)

def phase():
    c=1.
    a3=8.
    p3=a3/(27*rho)
    coeff=p3*(2/c)**.75
    rows=[]
    for j in [1,2,4,8]:
        N=4*round((2*math.pi*j/coeff)**4/4)
        center=math.sqrt(c*N/2)
        phases=[]
        for offset in (-N**.125,0,N**.125):
            e=(center+offset)**(-.5)
            def vals(U):
                A=cmath.exp(a3*(e*U)**3)
                Ap=3*a3*e**3*U**2*A
                return U+U**3*A,1+3*U**2*A+U**3*Ap
            h=newton(lambda u: vals(u)[1],1j/math.sqrt(3),tol=2e-14)
            W=vals(h)[0]
            phases.append(N*cmath.phase(W/(1j*rho))-2*math.pi*j)
        rows.append(dict(j=j,N=N,center_phase_error=phases[1],
                         window_phase_error=max(abs(u) for u in phases)))
    return rows
print(json.dumps(dict(python=platform.python_version(),scipy=scipy.__version__,
                      far=far_fold(),phase=phase()),indent=2))
