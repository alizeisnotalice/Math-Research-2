"""Binary64 Bessel cross-check and independent logarithmic cusp-to-smooth limit."""
import json,math
from pathlib import Path
import numpy as np
from scipy.integrate import quad
from scipy.special import ive,roots_legendre
from circle_mapped import integrate
from circle import integrate as original

k,j=np.meshgrid(np.arange(2,33),np.arange(1,65),indexing='ij')
k,j=k.ravel(),j.ravel()
p=np.linspace(.3,2,len(k))
a=np.full(len(k),2.)
y,_,_,L,_=integrate(k,j,p,a,24)
old,_,_,_,_=original(k,j,p,a,24)
bessel=L*ive(j,p*(k*j)**2)
assert np.all(np.isfinite(bessel))
intervals=np.minimum(np.array([0,.125,.5,1,2,3,5,8,14.])[None,:],np.pi*L[:,None])
n,w=roots_legendre(24)
width=np.diff(intervals)/2
midpoint=(intervals[:,:-1]+intervals[:,1:])/2
t=midpoint[:,:,None]+width[:,:,None]*n
b=np.sinc(t/(2*np.pi*L[:,None,None]))
scaled_direct=np.sum(w*width[:,:,None]*np.exp(-(t*b)**2)*np.cos(j[:,None,None]*t/L[:,None,None]),axis=(1,2))/np.pi
bessel_audit=dict(points=len(k),old_max_abs=float(max(abs(old-bessel))),
                 mapped_max_abs=float(max(abs(y-bessel))),smooth_coordinate_max_abs=float(max(abs(scaled_direct-bessel))))
rows=[]
for beta in (0.,1.,2.,4.):
    k,c=2,.875
    limit=math.exp(-math.exp(-beta)/(2*c*k*k))/(2*math.sqrt(math.pi))
    for logj in (8.,32.,128.,512.,2048.):
        alpha=2-beta/logj
        eta=math.exp((1-2/alpha)*logj-(2/alpha)*math.log(k)+math.log(2/c)/alpha)
        got=quad(lambda t:math.exp(-t**alpha)*math.cos(eta*t),0,32,epsabs=2e-13)[0]/math.pi
        rows.append(dict(beta=beta,log_harmonic=logj,alpha=alpha,eta=eta,reference=got,limit=limit,error=abs(got-limit)))
out=dict(bessel=bessel_audit,transition=rows)
Path('transition_bessel.json').write_text(json.dumps(out,indent=2)+"\n")
print(bessel_audit)
print([r["error"] for r in rows if r["log_harmonic"]==2048])
