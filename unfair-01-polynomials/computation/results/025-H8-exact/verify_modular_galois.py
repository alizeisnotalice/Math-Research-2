"""Independent finite-field factor certificates; Python integer arithmetic only."""
from math import isqrt,lcm
import hashlib,json

def trim(a):
 a=list(a)
 while len(a)>1 and a[-1]==0:a.pop()
 return a

def monic(a,p):
 a=trim([v%p for v in a]);assert a!=[0]
 z=pow(a[-1],-1,p);return [(v*z)%p for v in a]

def add(a,b,p,sign=1):
 return trim([((a[i] if i<len(a) else 0)+sign*(b[i] if i<len(b) else 0))%p for i in range(max(len(a),len(b)))])

def mul(a,b,p):
 c=[0]*(len(a)+len(b)-1)
 for i,v in enumerate(a):
  for j,w in enumerate(b):c[i+j]=(c[i+j]+v*w)%p
 return trim(c)

def rem(a,b,p):
 a=trim([v%p for v in a]);b=trim([v%p for v in b]);assert b!=[0];iv=pow(b[-1],-1,p)
 while a!=[0] and len(a)>=len(b):
  shift=len(a)-len(b);v=a[-1]*iv%p
  for j,w in enumerate(b):a[shift+j]=(a[shift+j]-v*w)%p
  a=trim(a)
 return a

def gcd(a,b,p):
 while b!=[0]:a,b=b,rem(a,b,p)
 return monic(a,p)

def power(a,n,f,p):
 r=[1];a=rem(a,f,p)
 while n:
  if n&1:r=rem(mul(r,a,p),f,p)
  a=rem(mul(a,a,p),f,p);n//=2
 return r

def prime(p):return p>=2 and all(p%d for d in range(2,isqrt(p)+1))
def prime_divisors(n):return [q for q in range(2,n+1) if n%q==0 and prime(q)]

def irreducible(f,p):
 assert prime(p) and len(f)>=2 and f[-1]%p==1
 n=len(f)-1;x=rem([0,1],f,p);v=x;targets={n//q for q in prime_divisors(n)}
 for j in range(1,n+1):
  v=power(v,p,f,p)
  if j in targets and gcd(f,add(v,x,p,-1),p)!=[1]:return False
 return v==x

def reduction(C,row):
 p=row['prime'];assert prime(p)
 prod=[1];degrees=[]
 for item in row['factors']:
  f=item['coefficients_ascending'];assert item.get('multiplicity',1)==1
  assert all(type(v)==int and 0<=v<p for v in f) and f[-1]==1
  assert irreducible(f,p),'factor not irreducible'
  prod=mul(prod,f,p);degrees.append(len(f)-1)
 assert prod==trim([v%p for v in C]),'factor product mismatch'
 derivative=trim([(i*C[i])%p for i in range(1,len(C))])
 assert gcd(prod,derivative,p)==[1],'not squarefree'
 return sorted(degrees)

def verify(rec):
 C=rec['coefficients_ascending'];N=len(C)-1;assert N>=3 and C[0]==C[-1]==1 and set(C)<={0,1}
 assert rec['degree']==N and all(type(x)==int for x in C)
 assert hashlib.sha256(json.dumps(C,separators=(',',':')).encode()).hexdigest()==rec['sha256']
 witnessed=set();details=[]
 for row in rec['reductions']:
  ds=reduction(C,row);roles=[]
  if ds==[N]:witnessed.add('transitive');roles.append('transitive')
  if ds==[1,N-1]:witnessed.add('long_cycle');roles.append('long_cycle')
  if ds.count(2)==1 and all(d==2 or d%2 for d in ds):
   exponent=lcm(*(d for d in ds if d!=2));assert exponent%2
   witnessed.add('transposition');roles.append('transposition')
  details.append({'prime':row['prime'],'degrees':ds,'verified_roles':roles})
 assert witnessed=={'transitive','long_cycle','transposition'},witnessed
 return {'name':rec.get('name'),'degree':N,'status':'full_symmetric_by_exact_modular_certificates','reductions':details,'mathematical_dependency':'Frobenius factor-degree cycle theorem at squarefree reduction; transitivity plus (N-1)-cycle implies 2-transitivity; conjugates of one transposition give all transpositions.'}

if __name__=='__main__':
 import json,sys
 raw=json.load(open(sys.argv[1]));records=raw['certificates'] if isinstance(raw,dict) else raw
 print(json.dumps([verify(r) for r in records],indent=2))
