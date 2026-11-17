from pathlib import Path
from fractions import Fraction
import json,time,hashlib,platform
start=time.monotonic()
def fwht(a):
 a=list(a);h=1
 while h<len(a):
  for i in range(0,len(a),2*h):
   for j in range(i,i+h):a[j],a[j+h]=a[j]+a[j+h],a[j]-a[j+h]
  h*=2
 return a
rows=[]
for mu in range(3,7):
 M=1<<(2*mu);Q=lambda x:sum(((x>>(2*j))&1)*((x>>(2*j+1))&1) for j in range(mu))%2
 B=[x for x in range(M) if Q(x)];C=[int(not Q(x)) for x in range(M)];W=fwht(C)
 assert W[0]==(M+(1<<mu))//2 and {abs(x) for x in W[1:]}=={1<<(mu-1)}
 T=1<<(2*mu-2);bound=Fraction(M*len(B),T*T);assert bound<8
 floor=bound.numerator//bound.denominator;lower=W[0]-(floor-1)*(1<<(mu-1));assert lower>0
 N=1<<mu
 extrema=[r for r in range(1,N) if r*(N-r)<=N-1];assert extrema==[1,N-1]
 a=[0]*M
 for x in range(1,N):a[x]=1;a[x*N]=-1
 W=fwht(a);assert set(W)<={0,N,-N};assert sum(x!=0 for x in W)==2*(N-1)
 rows.append({'mu':mu,'type1_core_size':len(B),'type1_spectral_count_bound':str(bound),'type1_gram_lower_bound':lower,'type2_core_size':2*(N-1),'type2_reference_spectrum_count':sum(x!=0 for x in W),'type2_extremal_r':extrema})
out={'rows':rows,'arithmetic':'exact integers and fractions','python':platform.python_version(),'elapsed_seconds':time.monotonic()-start,'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
Path(__file__).with_suffix('.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
