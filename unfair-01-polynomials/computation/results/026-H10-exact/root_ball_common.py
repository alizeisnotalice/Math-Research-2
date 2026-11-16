"""Exact rational complex and Rouché bounds for the independent H10 checker."""
from fractions import Fraction as F
from math import comb,isqrt
def add(a,b):return a[0]+b[0],a[1]+b[1]
def sub(a,b):return a[0]-b[0],a[1]-b[1]
def mul(a,b):return a[0]*b[0]-a[1]*b[1],a[0]*b[1]+a[1]*b[0]
def scale(a,q):return a[0]*q,a[1]*q
def powc(a,n):
 r=(F(1),F(0))
 while n:
  if n&1:r=mul(r,a)
  a=mul(a,a);n//=2
 return r
def sq(a):return a[0]*a[0]+a[1]*a[1]
def parse_fraction(s):
 a,b=s.split('/',1);return F(int(a),int(b))
def complex_from_json(x):return parse_fraction(x['re']),parse_fraction(x['im'])
def complex_from_strings(x):return F(x[0]),F(x[1])
def abs_bounds(a,scale_n):
 q=sq(a);t=(q.numerator*scale_n*scale_n)//q.denominator;n=isqrt(t);return F(n,scale_n),F(n+1,scale_n)
def taylor(c,z):
 out=[]
 for k in range(len(c)):
  s=(F(0),F(0))
  for n in range(k,len(c)):
   if c[n]:s=add(s,scale(powc(z,n-k),c[n]*comb(n,k)))
  out.append(s)
 return out
def polynomial_rouche_bounds(c,z,r,scale_n):
 t=taylor(c,z);_,L=abs_bounds(t[0],scale_n)
 for k in range(2,len(t)):_,u=abs_bounds(t[k],scale_n);L+=u*r**k
 l,_=abs_bounds(t[1],scale_n);return L,l*r
def lift_rouche_bounds(ac,ar,z,r,scale_n,q):
 _,L=abs_bounds(sub(powc(z,q),ac),scale_n);L+=ar
 for k in range(2,q+1):_,u=abs_bounds(scale(powc(z,q-k),comb(q,k)),scale_n);L+=u*r**k
 l,_=abs_bounds(scale(powc(z,q-1),q),scale_n);return L,l*r
def disk_disjoint(a,r,b,s):return sq(sub(a,b))>(r+s)**2
