"""Rational Taylor interval quadrature for the fixed-angle AH Jordan region."""
from pathlib import Path
from sigma_parameter import SIGMA_OUT
from fractions import Fraction as F
import runpy,io,contextlib,json,sys,hashlib
with contextlib.redirect_stdout(io.StringIO()):d=runpy.run_path(str(Path(__file__).with_name('certify_ah_path.py')))
I=d['I'];order=3
class J:
 def __init__(self,a):
  self.c=a.c[:] if isinstance(a,J) else [I(a)]+[I(0) for _ in range(order)]
 @classmethod
 def make(cls,a):
  q=cls(0);q.c=a;return q
 def __add__(self,b):
  b=J(b);return J.make([x+y for x,y in zip(self.c,b.c)])
 __radd__=__add__
 def __neg__(self):return J.make([-x for x in self.c])
 def __sub__(self,b):return self+-J(b)
 def __rsub__(self,b):return J(b)+-self
 def __mul__(self,b):
  b=J(b);return J.make([sum((self.c[k]*b.c[n-k] for k in range(n+1)),I(0)) for n in range(order+1)])
 __rmul__=__mul__
 def inv(self):
  out=[I(1)/self.c[0]]
  for n in range(1,order+1):out.append(-sum((self.c[k]*out[n-k] for k in range(1,n+1)),I(0))/self.c[0])
  return J.make(out)
 def __truediv__(self,b):return self*J(b).inv()
 def sqrt(self):
  out=[self.c[0].sqrt()]
  for n in range(1,order+1):out.append((self.c[n]-sum((out[k]*out[n-k] for k in range(1,n)),I(0)))/(2*out[0]))
  return J.make(out)
def vec(p):return tuple(J(x) for x in p)
def add(a,b):return tuple(x+y for x,y in zip(a,b))
def sub(a,b):return tuple(x-y for x,y in zip(a,b))
def mul(t,a):return tuple(t*x for x in a)
def dot(a,b):return sum((x*y for x,y in zip(a,b)),J(0))
def perp(a):return (-a[1],a[0])
r3=J(d['r3']);half=J(F(1,2));M=vec(d['M']);Z=vec(d['Z']);T=vec(d['T'])
nL=(half,-r3/2);nN=(-half,-r3/2);eL=(r3/2,half);eN=(-r3/2,half)
def paired(X):
 def inverse(n,e):
  dist=dot(n,X)-half;foot=sub(X,mul(dist,n));return sub(foot,mul((1-dist*dist).sqrt(),e))
 L=inverse(nL,eL);N=inverse(nN,eN)
 c=L[0]+half;lx=half*c+r3/2*(1-c*c).sqrt()-half
 Lp=(lx,(lx-1)/r3);Np=sub(mul(J(2),M),N)
 q=sub(Np,Lp);dsq=dot(q,q)
 return sub(mul(half,add(Lp,Np)),mul((1-dsq/4).sqrt()/dsq.sqrt(),perp(q)))
def coefficients(interval):
 w=J(interval);w.c[1]=I(1)
 P=paired(add(Z,mul(w,sub(T,Z))))
 # Green density x y' - y x', Taylor coefficients through degree 2.
 return [sum(((j+1)*(P[0].c[n-j]*P[1].c[j+1]-P[1].c[n-j]*P[0].c[j+1]) for j in range(n+1)),I(0)) for n in range(3)]
def crossI(a,b):return a[0]*b[1]-a[1]*b[0]
def atanI(z):
 assert max(abs(z.lo),abs(z.hi))<F(1,2)
 result=I(0);power=z
 for k in range(32):
  result=result+((-1)**k)*power/(2*k+1);power=power*z*z
 rem=max(abs(z.lo),abs(z.hi))**65/65
 return result+I(-rem,rem)
def arc(a,b,c):
 va=d['sub'](a,c);vb=d['sub'](b,c);dp=d['dot'](va,vb);assert dp.lo>0
 return crossI(c,d['sub'](b,a))+atanI(crossI(va,vb)/dp)
Tp=tuple(x.c[0] for x in paired(T));Xas=d['upper_cross'](d['G'],d['E'])
fixed=arc(Tp,Xas,d['G'])+arc(Xas,d['Y'],d['E'])+arc(d['Y'],d['T'],d['C'])+crossI(d['T'],d['Z'])
N=int(sys.argv[1]) if len(sys.argv)>1 else 128
assert N>=1
integral=I(0);max_f2=F(0)
for j in range(N):
 cell=I(F(j,N),F(j+1,N));mid=I(F(2*j+1,2*N))
 f0=coefficients(mid)[0];f2=coefficients(cell)[2]
 integral=integral+f0/F(N)+f2/F(12*N**3)
 max_f2=max(max_f2,abs(f2.lo),abs(f2.hi))
signed=(fixed+integral)/2
assert signed.hi<0 or signed.lo>0
area=-signed if signed.hi<0 else signed
print(json.dumps({'python':sys.version,'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),"sigma_degrees_interval": SIGMA_OUT,'cells':N,'jet_order':order,'midpoint_remainder_formula':'width*f(mid) + width^3/12 * enclosure(f_second_derivative/2)','area_interval':area.out(),'area_display_only':[float(area.lo),float(area.hi)],'width':str(area.hi-area.lo),'maximum_abs_density_quadratic_coefficient':str(max_f2),'paired_tangent_point':[x.out() for x in Tp],'scope':'area of the fixed-angle candidate AH Jordan region only, not an independently proved universal cover','passed':True},indent=2))
