"""High derivative cusp exclusion: exact tail budget and numerical asymptotics."""
from fractions import Fraction as Q
from pathlib import Path
import json,os,platform,hashlib,time
import mpmath as mp
start=time.monotonic()
epsilon=Q(1,2**60)
# On |u|<=eps: gamma density>1/2187 and negative kernel curvature
# integrates to magnitude > eps^-1/2 /34992.
core_magnitude=Q(2**30,34992);tail_positive=Q(30)*4**3
assert epsilon<Q(1,4)
assert core_magnitude-tail_positive>20000
# Positive kernel curvature outside |u|>=1/4 is bounded by 30/|u|^3:
# max w^(5/2) exp(-w/2) occurs at w=5; sqrt(3125)<60.
assert 5**5<60**2
rows=[]
for dps in [60,100]:
 mp.mp.dps=dps
 normal=lambda x:mp.exp(-x*x/2)/mp.sqrt(2*mp.pi)
 for alpha_numer,alpha_denom,n in [(1,1,1),(1,3,3)]:
  alpha=mp.mpf(alpha_numer)/alpha_denom;p=mp.mpf(n)+mp.mpf('0.5')
  he=(lambda z:z*z-1) if n==1 else (lambda z:z**6-15*z**4+45*z*z-15)
  h0=he(mp.mpf(0));pred=2*normal(3)/mp.sqrt(2*mp.pi)*h0/alpha*mp.beta(1/alpha,p-1/alpha)
  for power in [20,60,120]:
   eps=mp.mpf(2)**(-power);r=mp.mpf(1)/4;S=mp.log1p(r**alpha/eps)
   def integrand(s):
    e1=mp.expm1(s);u=(eps*e1)**(1/alpha);z=u/mp.sqrt(eps*mp.exp(s))
    jac=mp.exp((1-p)*s)*e1**(1/alpha-1)/alpha
    return (normal(3+u)+normal(3-u))*he(z)*mp.exp(-z*z/2)*jac/mp.sqrt(2*mp.pi)
   knots=sorted(set([mp.mpf(0),mp.mpf(1),mp.mpf(5),S/2,S-5,S]))
   knots=[x for x in knots if 0<=x<=S]
   scaled=mp.quad(integrand,knots)
   actual=scaled*eps**(1/alpha-p)
   rows.append(dict(dps=dps,alpha=f'{alpha_numer}/{alpha_denom}',derivative_order=2*n,power=power,core_derivative=str(actual),scaled_derivative=str(scaled),limit_coefficient=str(pred),ratio_to_limit=str(scaled/pred)))
for alpha in ['1/1','1/3']:
 for power in [20,60,120]:
  a,b=[r for r in rows if r['alpha']==alpha and r['power']==power]
  assert abs(mp.mpf(a['scaled_derivative'])-mp.mpf(b['scaled_derivative']))<mp.mpf('1e-50')
out=dict(parameters=dict(M=2,x0=3,r='1/4',certificate_deficit='M-f(x)=abs(x-3)',epsilon=str(epsilon)),certificate=dict(core_negative_magnitude_lower=str(core_magnitude),outside_positive_upper=str(tail_positive),full_second_derivative_upper=str(tail_positive-core_magnitude),certified_upper='-20000',target_second_derivative_positive=True),diagnostics=rows,runtime=dict(python=platform.python_version(),mpmath=mp.__version__,pid=os.getpid(),elapsed_seconds=time.monotonic()-start),script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),scope='certificate uses analytic kernel bounds and exact Fraction; high-precision quadrature only checks asymptotics')
Path(__file__).with_name('round017-cusp-results.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out['certificate'],indent=2));print(out['runtime'])
