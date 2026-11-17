"""Quadratic maximum local density: rational witness and transformed quadrature."""
from fractions import Fraction as Q
from pathlib import Path
import json,platform,os,time,hashlib
import mpmath as mp
start=time.monotonic()
# Elementary constants used below. Polynomial division gives
# integral_0^1 x^4(1-x)^4/(1+x^2) dx = 22/7 - pi > 0.
assert Q(1,7)-Q(4,6)+Q(5,5)-Q(4,3)+4==Q(22,7)
assert 4*sum(Q((-1)**j,2*j+1) for j in range(8))>3
assert Q(49,25)<2
assert sum(Q(2)**j / __import__('math').factorial(j) for j in range(5))==7
# M=2, x0=3, r=1/4, epsilon=2^-600. On sqrt(eps)<=|u|<=r,
# source density >= sqrt(2)/(2pi)*exp(-185/32)*log(r/sqrt(eps)).
# sqrt2>7/5, 2pi<44/7, exp(185/32)<3^6, log2>2/3.
annulus_lower=Q(49,220)*Q(1,729)*298*Q(2,3)
target_upper=Q(1,21);core_floor=Q(1,20)
assert annulus_lower>core_floor>target_upper
assert Q(2)+Q(1,2**600)<Q(9,4)
rows=[]
for dps in [80,220]:
 mp.mp.dps=dps
 normal=lambda x:mp.exp(-x*x/2)/mp.sqrt(2*mp.pi)
 coeff=normal(mp.mpf(3))*mp.exp(-mp.mpf(1)/2)/mp.sqrt(2*mp.pi)
 for power in [20,100,300,600]:
  eps=mp.mpf(2)**(-power);root=mp.sqrt(eps);T=mp.asinh(mp.mpf(1)/(4*root))
  # u=sqrt(eps)*sinh(t) removes the singular Jacobian exactly.
  def both(t):
   u=root*mp.sinh(t)
   return (normal(3+u)+normal(3-u))*mp.exp(-mp.tanh(t)**2/2)/mp.sqrt(2*mp.pi)
  knots=sorted(set([mp.mpf(0),mp.mpf(1),mp.mpf(5),T-5,T]))
  knots=[x for x in knots if 0<=x<=T]
  core=mp.quad(both,knots);H=2+eps
  target=mp.exp(-9/(2*H))/mp.sqrt(2*mp.pi*H)
  rows.append(dict(dps=dps,power=power,core_density=str(core),gaussian_density=str(target),core_over_target=str(core/target),asymptotic_coefficient=str(coeff),core_over_log=str(core/mp.log(1/eps))))
for n in [20,100,300,600]:
 a,b=[r for r in rows if r['power']==n]
 assert abs(mp.mpf(a['core_density'])-mp.mpf(b['core_density']))<mp.mpf('1e-70')
out=dict(parameters=dict(M=2,x0=3,r='1/4',local_deficit='M-f(x)=(x-3)^2',epsilon='2^-600'),certificate=dict(annulus_lower=str(annulus_lower),core_floor=str(core_floor),target_upper=str(target_upper),gap_floor=str(core_floor-target_upper),certified=True),diagnostics=rows,runtime=dict(python=platform.python_version(),mpmath=mp.__version__,pid=os.getpid(),elapsed_seconds=time.monotonic()-start),script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),scope='exact rational witness plus analytic inequalities; quadrature is diagnostic only')
Path(__file__).with_name('round016-contact-results.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out,indent=2))
