"""Exact rational cosine core certificate; independent high-precision diagnostics."""
from fractions import Fraction as Q
from pathlib import Path
import json,hashlib,os,platform,time
import mpmath as mp
start=time.monotonic()
xi=Q(1,1000);rho=Q(1,10**12);R=10;H=2
# e_*=(3/5) exp(-16/3) > (3/5)/3**6; D >= 7 e_*/160.
e_lower=Q(3,5)/3**6;D_lower=Q(7,160)*e_lower
sixth_bound=Q(15)+45*Q(3,2)+45*Q(3,2)**2+15*Q(3,2)**3
assert sixth_bound==Q(1875,8)
full_gap=D_lower*xi**4/8-(sixth_bound+120)*xi**6/720
# 2*P(|G|>10) < (1/5)*(3/8)**50 using sqrt(2pi)>2 and e>8/3.
tail_loss=Q(1,5)*Q(3,8)**50
perturbation_loss=xi**2*rho/2
core_gap=full_gap-tail_loss-perturbation_loss
margin=Q(3,10**18)
assert core_gap>margin
exact={k:str(v) for k,v in dict(xi=xi,rho=rho,e_lower=e_lower,D_lower=D_lower,sixth_bound=sixth_bound,full_gap=full_gap,tail_loss=tail_loss,perturbation_loss=perturbation_loss,core_gap=core_gap,claimed_margin=margin).items()}
# Integrals are diagnostics, not used by the rational certificate.
rows=[]
for dps in [50,80]:
 mp.mp.dps=dps
 e=mp.mpf(3)/5*mp.exp(-mp.mpf(16)/3);xx=mp.mpf(1)/1000
 normal=lambda x:mp.exp(-x*x/2)/mp.sqrt(2*mp.pi)
 def ratio(v,x):
  return mp.exp(x*x/2-(x-2)**2/(2*v))/mp.sqrt(v)
 def f(x):return 1+e*(ratio(mp.mpf(5)/8,x)-ratio(mp.mpf(1)/8,x))/2
 points=[-10,-6,-3,0,1,2,3,4,5,6,8,10]
 D=e/2-mp.quad(lambda x:(f(x)-1)**2*normal(x),points)
 # Stable evaluation of 1-exp(-a)*cos(b).
 def test(x):
  a=(2-f(x))*xx*xx/2;b=xx*x
  return -mp.expm1(-a)+mp.exp(-a)*2*mp.sin(b/2)**2
 core=mp.quad(lambda x:test(x)*normal(x),points)
 target=-mp.expm1(-xx*xx)
 gap=core-target
 rows.append(dict(dps=dps,D_truncated_display=str(D),core_test=str(core),gaussian_test=str(target),core_excess=str(gap)))
assert abs(mp.mpf(rows[0]['core_excess'])-mp.mpf(rows[1]['core_excess']))<mp.mpf('1e-45')
out=dict(parameters=dict(H=H,R=R,core='fbar from round013 section5',bound='0<=f<=2 globally; abs(f-fbar)<=1e-12 on [-10,10]'),exact_certificate=exact,diagnostics=rows,certified=True,runtime=dict(python=platform.python_version(),mpmath=mp.__version__,pid=os.getpid(),elapsed_seconds=time.monotonic()-start),script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),scope='certificate uses Fraction and analytic moment/tail estimates; quadrature is numerical only')
Path(__file__).with_name('round015-fourier-results.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out,indent=2))
