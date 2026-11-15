import math,cmath,json,sys
from fractions import Fraction as Q
pi=math.pi
def ratio(h,mode,direction):
 x,y=h*direction[0],h*direction[1];psi=abs(x)+abs(y)
 if mode=="compact":v=math.cos(3*x)
 elif mode=="time":v=math.exp(-psi)
 elif mode=="shift":v=cmath.exp(-psi-1j*x)
 elif mode=="centered_shift":v=math.exp(-psi)*math.cos(3*x)
 elif mode=="anisotropic":v=math.exp(-abs(x)-2*abs(y))
 return 2*pi*(1-v)/psi
rows=[]
for mode in ["compact","time","shift","centered_shift","anisotropic"]:
 for d in [(1,0),(0,1),(-1,0),(1,1)]:
  v=ratio(1e-5,mode,d)
  rows.append(dict(mode=mode,direction=d,ratio_real=float(v.real),ratio_imag=float(v.imag)))
weights=[Q(1,4)]*4;times=[Q(1,2),Q(1,2),Q(3,2),Q(3,2)];locs=[(3,0),(-3,0),(0,3),(0,-3)]
assert sum(weights)==1 and all(sum(w*z[j] for w,z in zip(weights,locs))==0 for j in [0,1])
T=sum(w*t for w,t in zip(weights,times));assert T==1
ms=[]
for t in [.5,1.,2.]:
 m=(math.atan(t)-t/(1+t*t))/pi
 assert m>0
 ms.append(dict(t=t,proven_uniform_D_lower=m,exact_mass_float=2*pi*t))
print(json.dumps(dict(python=sys.version,scope="scalar Fourier direction sanity, exact location/time moments, no kernel optimization",finite_heavy_mix_mean_time=str(T),finite_heavy_mix_total_mass="2*pi",unshifted_time_lower_bounds=ms,low_frequency_h=1e-5,Fourier_direction_samples=rows),indent=2))
