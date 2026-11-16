from fractions import Fraction as F
from math import comb,isqrt
def pf(s):a,b=s.split('/',1);return F(int(a),int(b))
def cj(x):return pf(x['re']),pf(x['im'])
def sub(a,b):return a[0]-b[0],a[1]-b[1]
def mul(a,b):return a[0]*b[0]-a[1]*b[1],a[0]*b[1]+a[1]*b[0]
def pw(a,n):
 r=(F(1),F(0))
 while n:
  if n&1:r=mul(r,a)
  a=mul(a,a);n//=2
 return r
def ab(a,s):
 q=a[0]*a[0]+a[1]*a[1];t=(q.numerator*s*s)//q.denominator;n=isqrt(t);return F(n,s),F(n+1,s)
def tay(c,z):
 out=[]
 for k in range(len(c)):
  a=(F(0),F(0))
  for n in range(k,len(c)):
   w=pw(z,n-k);a=(a[0]+c[n]*comb(n,k)*w[0],a[1]+c[n]*comb(n,k)*w[1])
  out.append(a)
 return out
def mother(c,z,r,s):
 t=tay(c,z);_,L=ab(t[0],s)
 for k in range(2,len(t)):_,u=ab(t[k],s);L+=u*r**k
 l,_=ab(t[1],s);return L,l*r
def lift(a,ar,z,r,s,q):
 _,L=ab(sub(pw(z,q),a),s);L+=ar
 for k in range(2,q+1):_,u=ab((comb(q,k)*pw(z,q-k)[0],comb(q,k)*pw(z,q-k)[1]),s);L+=u*r**k
 l,_=ab((q*pw(z,q-1)[0],q*pw(z,q-1)[1]),s);return L,l*r
def dis(a,r,b,s):return (a[0]-b[0])**2+(a[1]-b[1])**2>(r+s)**2
