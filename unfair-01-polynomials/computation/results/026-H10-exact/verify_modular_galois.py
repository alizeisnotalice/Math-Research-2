"""Integer finite-field verifier for the fixed H10 S10 certificate."""
import hashlib,json
from math import isqrt,lcm
def trim(a):
 a=list(a)
 while len(a)>1 and a[-1]==0:a.pop()
 return a
def mul(a,b,p):
 z=[0]*(len(a)+len(b)-1)
 for i,x in enumerate(a):
  for j,y in enumerate(b):z[i+j]=(z[i+j]+x*y)%p
 return trim(z)
def rem(a,b,p):
 a=trim([x%p for x in a]);b=trim([x%p for x in b]);iv=pow(b[-1],-1,p)
 while a!=[0] and len(a)>=len(b):
  s=len(a)-len(b);q=a[-1]*iv%p
  for j,x in enumerate(b):a[s+j]=(a[s+j]-q*x)%p
  a=trim(a)
 return a
def gcd(a,b,p):
 while b!=[0]:a,b=b,rem(a,b,p)
 iv=pow(a[-1],-1,p);return [(x*iv)%p for x in a]
def power(a,n,f,p):
 r=[1];a=rem(a,f,p)
 while n:
  if n&1:r=rem(mul(r,a,p),f,p)
  a=rem(mul(a,a,p),f,p);n//=2
 return r
def irr(f,p):
 n=len(f)-1;x=rem([0,1],f,p);v=x
 for j in range(1,n+1):
  v=power(v,p,f,p)
  if n%j==0 and j<n and gcd(f,[(v[i]-(x[i] if i<len(x) else 0))%p for i in range(max(len(v),len(x)))],p)!=[1]:return False
 return v==x
def verify(r):
 C=r['coefficients_ascending'];assert len(C)==11 and C[0]==C[-1]==1 and hashlib.sha256(json.dumps(C,separators=(',',':')).encode()).hexdigest()==r['sha256']
 kinds=[]
 for row in r['reductions']:
  p=row['prime'];assert p>1 and all(p%d for d in range(2,isqrt(p)+1));pr=[1];ds=[]
  for f in row['factors']:
   a=f['coefficients_ascending'];assert f.get('multiplicity',1)==1 and a[-1]==1 and irr(a,p);pr=mul(pr,a,p);ds.append(len(a)-1)
  assert pr==[x%p for x in C] and gcd(pr,[(i*C[i])%p for i in range(1,len(C))],p)==[1]
  kinds.append(sorted(ds))
 assert sorted(kinds)==[[1,9],[2,3,5],[10]]
 return {'status':'S10_by_exact_modular_certificates','reductions':kinds}
