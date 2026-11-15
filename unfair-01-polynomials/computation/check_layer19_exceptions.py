"""Independent integer-polynomial reconstruction of three support contradictions."""
from pathlib import Path
import json
N=3
Z=(0,)*N

def const(c):return {Z:c} if c else {}
def var(i):m=list(Z);m[i]=1;return {tuple(m):1}
def add(a,b):
 o=dict(a)
 for m,c in b.items():o[m]=o.get(m,0)+c
 return {m:c for m,c in o.items() if c}
def scale(a,c):return {m:v*c for m,v in a.items() if v*c}
def sub(a,b):return add(a,scale(b,-1))
def mul(a,b):
 o={}
 for i,c in a.items():
  for j,d in b.items():
   m=tuple(x+y for x,y in zip(i,j));o[m]=o.get(m,0)+c*d
 return {m:c for m,c in o.items() if c}
def power(a,n):
 o=const(1)
 for _ in range(n):o=mul(o,a)
 return o
def subst(a,values):
 o={}
 for m,c in a.items():
  t=const(c)
  for i,e in enumerate(m):t=mul(t,power(values.get(i,var(i)),e))
  o=add(o,t)
 return o
def up(cs):
 return {(i,0,0):v for i,v in enumerate(cs) if v}
def vector(a):
 assert all(m[1:]==(0,0) for m in a)
 return [a.get((i,0,0),0) for i in range(max((m[0] for m in a),default=0)+1)]
def prefix(I,J,A,end):
 sums={i+j for i in I for j in J};bs=[const(1)]
 for k in range(1,end+1):
  b=const(int(k in sums))
  for i in I:
   if 0<i<=k:b=sub(b,mul(A[i],bs[k-i]))
  bs.append(b)
 return bs
p,q,r=var(0),var(1),var(2);one=const(1)
I1=[0,1,2,5,9,13];J1=[0,1,2,3,6,7,10,14,18]
b=prefix(I1,J1,{0:one,1:p,2:q},4)
lin=up([-1,2,-3]);tail=up([1,-1,1,-1,1]);assert b[4]==add(add(power(q,2),mul(lin,q)),tail)
disc=sub(power(lin,2),scale(tail,4));assert disc==mul(up([-1,1]),up([3,3,-3,5]))
# p in (0,1): second discriminant factor =3+3p(1-p)+5p^3>0.
I2=[0,1,3,8,13];J2=[0,1,2,4,7,9,14,19]
b=prefix(I2,J2,{0:one,1:p,3:q,8:r},10)
qof=up([1,-1,1,-1]);assert b[3]==sub(qof,q)
post=[subst(x,{1:qof}) for x in b];rof=add(post[8],r);assert all(m[1:]==(0,0) for m in rof)
assert subst(b[8],{1:qof,2:rof})=={}
K=up([-2,3,-2,2]);H=up([1,-6,18,-31,44,-44,32,-15,5,1])
assert post[5]==mul(mul(up([-1,1]),p),K)
assert subst(b[10],{1:qof,2:rof})==mul(up([-1,1]),H)
U=up([-172,-962,1003,-2320,2664,-2922,1609,-991,-165]);V=up([235,2,330])
assert add(mul(U,K),mul(V,H))==const(579)
I3=[0,1,3,6,8,13];J3=[0,1,2,4,9,14,19]
def diag(k):return [(i,k-i) for i in I3 if k-i in J3]
assert diag(6)==[(6,0)]
assert diag(8)==[(6,2),(8,0)]
assert diag(27)==[(8,19),(13,14)]
assert diag(14)==[(0,14),(13,1)]
assert diag(1)==[(0,1),(1,0)]
assert diag(2)==[(0,2),(1,1)]
# Forced a6=1,b2=1-p,b14=1-p,b1=p,a1=1-p where p=a8>0.
assert sub(add(up([1,-1]),mul(up([1,-1]),p)),one)==scale(power(p,2),-1)
first_three={'cases':[{'I':I1,'J':J1,'certificate':'negative quadratic discriminant','discriminant_coefficients':vector(disc)},{'I':I2,'J':J2,'certificate':'integer Bezout after division-free substitutions','A3':vector(qof),'A8':vector(rof),'K':vector(K),'H':vector(H),'U':vector(U),'V':vector(V),'Bezout_constant':579},{'I':I3,'J':J3,'certificate':'six exact diagonals force a8 squared zero','diagonals':{str(k):diag(k) for k in (6,8,27,14,1,2)}}],'arithmetic':'pure Python integers; no CAS imports'}

# Case 4: derive missing B3,B4,B5 without fixing unneeded endpoints.
I4=[0,1,3,4,6,9,14];J4=[0,1,2,6,7,9,12,17]
b4=prefix(I4,J4,{0:one,1:p,3:q,4:r},5)
q4=up([1,-1,1,-1]);r4=up([0,2,-2,2,-1])
assert subst(b4[3],{1:q4})=={}
assert subst(b4[4],{1:q4,2:r4})=={}
assert subst(b4[5],{1:q4,2:r4})==up([0,0,1,-1,1])
assert all(k not in J4 for k in (3,4,5))

# Case 5: independent rational-function arithmetic over Z[p], common d=2p-1.
# No symbolic engine, solving routine, or generator import is used.
d=up([-1,2])
class R:
 def __init__(self,n,e=0):self.n,self.e=n,e
 def __add__(a,b):
  e=max(a.e,b.e)
  return R(add(mul(a.n,power(d,e-a.e)),mul(b.n,power(d,e-b.e))),e)
 def __neg__(a):return R(scale(a.n,-1),a.e)
 def __sub__(a,b):return a+-b
 def __mul__(a,b):return R(mul(a.n,b.n),a.e+b.e)
 def zero(a):return not a.n
 def equals(a,n,e=0):return sub(mul(a.n,power(d,e)),mul(n,power(d,a.e)))=={}
O=R(one);X=R(p)
A={0:O,1:X,14:O};B={0:O,17:O}
B[1]=O-X;B[2]=O-X*B[1];A[3]=O-X*B[2]
# c5 = g - (2p-1) a4; g is derived at a4=0.
g=X*(O-A[3]*B[1])+A[3]*B[2]-O
assert g.e==0 and g.n==up([0,-2,5,-5,4,-2])
from fractions import Fraction
assert sum(Fraction(v)*Fraction(1,2)**m[0] for m,v in g.n.items())==Fraction(-3,16)
A[4]=R(g.n,1);B[4]=O-A[4]-A[3]*B[1]
B[6]=O-A[3];A[6]=O-B[6]-A[4]*B[2]
B[9]=O-A[6];A[9]=O-B[9]-A[3]*B[6]
B[12]=O-A[9];B[7]=O-A[4]-A[9]*B[12]
I5=sorted(A);J5=sorted(B)
cert=json.loads((Path(__file__).resolve().parent/'results/019-case5-bezout.json').read_text())
assert (I5,J5)==(cert['I'],cert['J'])
def residual(k):
 s=-O
 for i in I5:
  if k-i in B:s=s+A[i]*B[k-i]
 return s
# All substitutions must actually satisfy their defining support diagonals.
used=[1,2,3,4,5,20,6,23,9,26,21]
assert all(residual(k).zero() for k in used)
P,Q,U,V=[up(cert[k+'_ascending']) for k in ('P','Q','U','V')]
assert residual(7).equals(P,2)
assert residual(18).equals(scale(Q,-1),2)
assert add(mul(U,P),mul(V,Q))==const(cert['rhs']) and cert['rhs']!=0

def check():
 return {'verified':True,'arithmetic':'pure Python integer polynomials and Fraction; no CAS',
 'cases':first_three['cases']+[
 {'I':I4,'J':J4,'certificate':'B5=p^2(p^2-p+1)>0 contradicts support','A3':vector(q4),'A4':vector(r4)},
 {'I':I5,'J':J5,'certificate':'support-derived substitutions and integer Bezout',
  'used_diagonals':used+[7,18],'zero_denominator_branch':'g(1/2)=-3/16',
  'bezout_constant':cert['rhs']} ]}
if __name__=='__main__':print(json.dumps(check(),indent=2))
