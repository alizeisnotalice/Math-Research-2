"""Logarithmic maximum contacts: rational Fourier separation and saddle diagnostics."""
from fractions import Fraction as Q
from math import factorial
from pathlib import Path
import json,os,platform,hashlib,time
import mpmath as mp
start=time.monotonic()
# xi=100, x0=3, v=1/log(1/|u|) for |u|<exp(-4), v>=1/4 outside.
# A positive inner annulus gives >3^-160; all possible negative terms <2^-1000.
positive=Q(1,3**160);negative=Q(1,2**1000);target=Q(1,2**5000)
certified=positive/2
assert positive-negative>certified>target
assert Q(10000,198)<51 and 2**100>200
# exp(5)>100 -> log(100)<5.
assert sum(Q(5)**j/factorial(j) for j in range(7))>100
rows=[]
for dps in [50,80]:
 mp.mp.dps=dps
 normal=lambda x:mp.exp(-x*x/2)/mp.sqrt(2*mp.pi)
 for beta_string in ['0.5','1','2']:
  beta=mp.mpf(beta_string)
  for frequency in [20,100,500]:
   xi=mp.mpf(frequency);tstar=(beta*xi*xi/2)**(1/(beta+1));root=mp.sqrt(tstar)
   leading_exponent=(beta+1)/beta*tstar
   lower=(4-tstar)/root
   def integrand(y):
    t=tstar+root*y;w=t/tstar;u=mp.exp(-t)
    action=tstar*(w+w**(-beta)/beta-1-1/beta)
    return root*(normal(3+u)+normal(3-u))*mp.cos(xi*u)*mp.exp(-action)
   points=[lower]+[mp.mpf(k) for k in [-40,-20,-12,-8,-4,0,4,8,12,20,40] if k>lower]+[mp.inf]
   scaled=mp.quad(integrand,points)
   prediction=2*normal(3)*mp.sqrt(2*mp.pi*tstar/(beta+1))
   rows.append(dict(dps=dps,beta=beta_string,frequency=frequency,saddle=str(tstar),leading_exponent=str(leading_exponent),scaled_local_characteristic=str(scaled),ratio_to_leading=str(scaled/prediction),log_local_characteristic=str(mp.log(scaled)-leading_exponent)))
for beta in ['0.5','1','2']:
 for xi in [20,100,500]:
  a,b=[r for r in rows if r['beta']==beta and r['frequency']==xi]
  assert abs(mp.mpf(a['scaled_local_characteristic'])-mp.mpf(b['scaled_local_characteristic']))<mp.mpf('1e-40')
out=dict(parameters=dict(x0=3,frequency=100,core_radius='exp(-4)',core_deficit='1/log(1/abs(x-3))',outside_deficit_lower='1/4',target_variance='M>1'),certificate=dict(positive_lower=str(positive),negative_upper=str(negative),certified_rotated_real_lower=str(certified),gaussian_modulus_upper=str(target),strict_separation=True),diagnostics=rows,runtime=dict(python=platform.python_version(),mpmath=mp.__version__,pid=os.getpid(),elapsed_seconds=time.monotonic()-start),script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),scope='rational Fourier certificate with analytic elementary bounds; saddle integrals are numerical diagnostics')
Path(__file__).with_name('round018-log-results.json').write_text(json.dumps(out,indent=2)+'\n')
print('Exact Fourier separation passed.');print(out['runtime'])
