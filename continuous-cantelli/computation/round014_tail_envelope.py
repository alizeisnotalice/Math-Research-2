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


# Core upper bound plus maximal possible nonnegative tail contribution.
filename='round004-L12-T12.json';lam=Q(2);rho=Q(1,10000000);margin=Q(8,10000000)
src=Path(__file__).with_name(filename);d=json.loads(src.read_text(),parse_float=str)
knots=list(map(Q,d['parameters']['knots']));vals=list(map(Q,d['best']['values']));R=knots[-1];lo=hi=0
for l,r,fl,fr in zip(knots,knots[1:],vals,vals[1:]):
 for aa,bb,fa,fb in [(l,r,fl,fr),(-r,-l,fr,fl)]:
  il,ih=segment(aa,bb,fa,fb,lam);lo+=il;hi+=ih
flo=lo*SCALE//sqrt_hi;fhi=(hi*SCALE+sqrt_lo-1)//sqrt_lo
# Mills bound Phi(-z)<=exp(-z^2/2)/(z*sqrt(2*pi)); z>0.
tail=0
for z in [R-lam,R+lam]:
 assert z>0
 el,eh=exp_iv(lam*lam/2-z*z/2)
 tail+=up(Q(eh,sqrt_lo)/z)
el,eh=exp_iv(lam*lam*rho/2);robust=(fhi*eh+SCALE-1)//SCALE+tail
assert Q(robust,SCALE)<1-margin
show=lambda x:str(mp.mpf(x)/SCALE)
certificate=dict(input_file=filename,input_sha256=hashlib.sha256(src.read_bytes()).hexdigest(),R=str(R),lambda_value=str(lam),rho=str(rho),margin=str(margin),core_lower=str(Q(flo,SCALE)),core_upper=str(Q(fhi,SCALE)),tail_upper=str(Q(tail,SCALE)),robust_total_upper=str(Q(robust,SCALE)),core_upper_display=show(fhi),tail_upper_display=show(tail),robust_upper_display=show(robust))

# Independent high-precision sharp Lipschitz tail envelopes for the L8 core.
mp.mp.dps=90
num=lambda x:mp.mpf(Q(x).numerator)/Q(x).denominator
def integral(a,b,fa,fb,l):
 mid=(a+b)/2;rad=(b-a)/2;v=(fb-fa)/(b-a);fm=(fa+fb)/2
 q=-mid*mid/2+l*mid-l*l*fm/2;aa=l-mid-l*l*v/2
 u=(-rad-aa)/mp.sqrt(2);w=(rad-aa)/mp.sqrt(2)
 diff=(mp.erfc(u)-mp.erfc(w)) if u>0 else ((mp.erfc(-w)-mp.erfc(-u)) if w<0 else mp.erf(w)-mp.erf(u))
 return mp.exp(q+aa*aa/2)*diff/2
Phi=lambda x:mp.erfc(-x/mp.sqrt(2))/2
source=Path(__file__).with_name('round004-L8-T12-long.json');dd=json.loads(source.read_text(),parse_float=str)
k=list(map(num,dd['parameters']['knots']));v=list(map(num,dd['best']['values']));RR=k[-1];ll=mp.mpf('2.73079');L=mp.mpf(11);H=mp.mpf(3);endpoint=v[-1]
core=mp.mpf(0)
for a,b,fa,fb in zip(k,k[1:],v,v[1:]):
 core+=integral(a,b,fa,fb,ll)+integral(-b,-a,fb,fa,ll)
def right_tail(l,maximum):
 if maximum:
  end=RR+endpoint/L
  return integral(RR,end,endpoint,mp.mpf(0),l)+mp.exp(l*l/2)*Phi(l-end)
 end=RR+(H-endpoint)/L
 return integral(RR,end,endpoint,H,l)+mp.exp(l*l*(1-H)/2)*Phi(l-end)
minimum=core+right_tail(ll,False)+right_tail(-ll,False)
maximum=core+right_tail(ll,True)+right_tail(-ll,True)
sharp=dict(input=source.name,input_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),lambda_value=str(ll),L=str(L),H=str(H),core=str(core),minimum=str(minimum),maximum=str(maximum),scope='90-digit displays of analytically sharp single-parameter envelope; no mean constraint imposed, no interval classification')
# Rational endpoint bounds certify the strict straddling, without relying on displays.
kk=list(map(Q,dd['parameters']['knots']));vv=list(map(Q,dd['best']['values']));qq=Q(273079,100000)
assert all(0<=v<=3 for v in vv) and vv[-1]==1
assert max(abs((vv[j+1]-vv[j])/(kk[j+1]-kk[j])) for j in range(len(kk)-1))<=11
cl=ch=0
for a,b,fa,fb in zip(kk,kk[1:],vv,vv[1:]):
 for aa,bb,va,vb in [(a,b,fa,fb),(-b,-a,fb,fa)]:
  il,ih=segment(aa,bb,va,vb,qq);cl+=il;ch+=ih
def normalized(il,ih):return il*SCALE//sqrt_hi,(ih*SCALE+sqrt_lo-1)//sqrt_lo
cl,ch=normalized(cl,ch)
def tail_iv(l,maximum):
 end=Q(8)+(Q(1,11) if maximum else Q(2,11));level=Q(0) if maximum else Q(3)
 il,ih=segment(Q(8),end,Q(1),level,l);il,ih=normalized(il,ih)
 z=end-l;assert z>0
 el,eh=exp_iv(l*l*(1-level)/2-z*z/2)
 il+=down(Q(el,sqrt_hi)*z/(1+z*z));ih+=up(Q(eh,sqrt_lo)/z)
 return il,ih
ml,mh=cl,ch;xl,xh=cl,ch
for l in [qq,-qq]:
 il,ih=tail_iv(l,False);ml+=il;mh+=ih
 il,ih=tail_iv(l,True);xl+=il;xh+=ih
assert mh<SCALE and xl>SCALE
straddling=dict(minimum_lower=str(Q(ml,SCALE)),minimum_upper=str(Q(mh,SCALE)),maximum_lower=str(Q(xl,SCALE)),maximum_upper=str(Q(xh,SCALE)),minimum_upper_display=show(mh),maximum_lower_display=show(xl),strict_straddling_certified=True)
out=dict(certificate=certificate,sharp_envelope=sharp,straddling_certificate=straddling,runtime=dict(python=platform.python_version(),mpmath=mp.__version__,pid=os.getpid()),script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),scope='core rejection is directed rational certificate; sharp-envelope decimals are numerical')
Path(sys.argv[1]).write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
