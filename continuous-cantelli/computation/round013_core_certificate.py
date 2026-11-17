"""Directed rational enclosures of a piecewise-linear core, with arbitrary nonnegative tails."""
from fractions import Fraction as Q
from math import factorial,ceil,isqrt
from pathlib import Path
import json,sys,os,platform,hashlib,time
import mpmath as mp
mp.mp.dps=70;BITS=180;SCALE=1<<BITS;ORDER=80

def down(x):return (x.numerator*SCALE)//x.denominator
def up(x):return -((-x.numerator*SCALE)//x.denominator)
def exp_iv(x):
 k=0;s=x
 while abs(s)>Q(1,2):s/=2;k+=1
 value=sum(s**j/factorial(j) for j in range(41));err=2*abs(s)**41/factorial(41)
 lo,hi=down(value-err),up(value+err);assert lo>0
 for _ in range(k):lo=lo*lo//SCALE;hi=(hi*hi+SCALE-1)//SCALE
 return lo,hi

def atan_iv(x):
 n=40;v=sum((-1)**j*x**(2*j+1)/Q(2*j+1) for j in range(n));e=x**(2*n+1)/Q(2*n+1)
 return v-e,v+e
# Machin identity checked through tangent addition with its usual first-quadrant branch.
tan2=lambda t:2*t/(1-t*t)
a=tan2(tan2(Q(1,5)));b=Q(1,239);assert (a-b)/(1+a*b)==1
al,ah=atan_iv(Q(1,5));bl,bh=atan_iv(Q(1,239))
pi_lo,pi_hi=down(16*al-4*bh),up(16*ah-4*bl)
sqrt_lo=isqrt(2*pi_lo*SCALE);sqrt_hi=isqrt(2*pi_hi*SCALE)+1
assert sqrt_lo*sqrt_lo<=2*pi_lo*SCALE and sqrt_hi*sqrt_hi>=2*pi_hi*SCALE

def segment(l,r,fl,fr,lam):
 mid=(l+r)/2;radius=(r-l)/2;assert 0<radius<1
 slope=(fr-fl)/(r-l);fm=(fl+fr)/2
 q0=-mid*mid/2+lam*mid-lam*lam*fm/2;a=lam-mid-lam*lam*slope/2
 cs=[Q(1),a]
 for j in range(2,ORDER+1):cs.append((a*cs[-1]-cs[-2])/j)
 integral=sum(cs[j]*2*radius**(j+1)/(j+1) for j in range(0,ORDER+1,2))
 remainder=2*radius*3**ceil(abs(a)+Q(1,2))*radius**(ORDER+1)/(1-radius)
 il=max(0,down(integral-remainder));ih=up(integral+remainder);el,eh=exp_iv(q0)
 return il*el//SCALE,(ih*eh+SCALE-1)//SCALE

rows=[];start=time.monotonic()
for filename,lam,rho,margin in [('round004-L8-T12-long.json',Q(2),Q(1,2000000),Q(1,10000000)),('round004-L12-T12.json',Q(5,2),Q(1,1000000),Q(3,1000000))]:
 src=Path(__file__).with_name(filename);d=json.loads(src.read_text(),parse_float=str)
 knots=list(map(Q,d['parameters']['knots']));vals=list(map(Q,d['best']['values']));lo=hi=0
 for l,r,fl,fr in zip(knots,knots[1:],vals,vals[1:]):
  for aa,bb,fa,fb in [(l,r,fl,fr),(-r,-l,fr,fl)]:
   il,ih=segment(aa,bb,fa,fb,lam);lo+=il;hi+=ih
 flo=lo*SCALE//sqrt_hi;fhi=(hi*SCALE+sqrt_lo-1)//sqrt_lo
 el,eh=exp_iv(-lam*lam*rho/2);robust=flo*el//SCALE
 assert Q(robust,SCALE)>1+margin
 display=lambda v:str(mp.mpf(v)/SCALE)
 rows.append(dict(input_file=filename,input_sha256=hashlib.sha256(src.read_bytes()).hexdigest(),R=str(knots[-1]),lambda_value=str(lam),core_tolerance=str(rho),certified_margin=str(margin),core_lower_dyadic=str(Q(flo,SCALE)),core_upper_dyadic=str(Q(fhi,SCALE)),robust_lower_dyadic=str(Q(robust,SCALE)),core_lower_display=display(flo),core_upper_display=display(fhi),robust_lower_display=display(robust),tail_scope='arbitrary measurable f>=0 outside the core; no upper bound or continuity needed for the lower certificate'))
 print(filename,display(flo),display(robust),flush=True)
out=dict(parameters=dict(dyadic_bits=BITS,local_Taylor_order=ORDER,exp_Taylor_order=40,atan_terms=40),certificates=rows,pi_interval=[str(Q(pi_lo,SCALE)),str(Q(pi_hi,SCALE))],sqrt_2pi_interval=[str(Q(sqrt_lo,SCALE)),str(Q(sqrt_hi,SCALE))],runtime=dict(python=platform.python_version(),mpmath=mp.__version__,pid=os.getpid(),elapsed_seconds=time.monotonic()-start),script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),scope='all asserted margins certified by directed rational arithmetic; decimal displays only are high precision approximations')
Path(sys.argv[1]).write_text(json.dumps(out,indent=2)+'\n')
