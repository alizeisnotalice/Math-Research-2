import itertools,json,hashlib,time
from pathlib import Path
start=time.monotonic()
def fwht(mask):
 a=[1-2*((mask>>x)&1) for x in range(256)]
 h=1
 while h<256:
  for p in range(0,256,2*h):
   for j in range(p,p+h):a[j],a[j+h]=a[j]+a[j+h],a[j]-a[j+h]
  h*=2
 return a
f=sum((((x&15)&(x>>4)).bit_count()%2)<<x for x in range(256))
assert set(map(abs,fwht(f)))=={16}
rows=[]
def add(sigma,h,kind):
 g=sum((((((x&15)&sigma[x>>4]).bit_count()%2)^((h>>(x>>4))&1)))<<x for x in range(256))
 A=f^g
 assert A.bit_count()==32 and set(map(abs,fwht(g)))=={16}
 rows.append({'kind':kind,'sigma':sigma,'h':h,'support_hex':hex(A)})
for ys in [(0,1,2,3),(0,1,2,4)]:
 for p in itertools.permutations(ys):
  if any(a==b for a,b in zip(ys,p)):continue
  sigma=list(range(16))
  for a,b in zip(ys,p):sigma[a]=b
  for bits in range(16):add(sigma,sum(((bits>>i)&1)<<y for i,y in enumerate(ys)),str(ys))
sigma=list(range(16));sigma[0],sigma[1]=1,0
for bits in range(4):add(sigma,bits|4,'transposition_plus_full')
add(list(range(16)),3,'two_full')
assert len(rows)==293
out={'schema':'a28-r020-standard-mm-weight32-cases-v1','scope':'standard base x dot y only','cases':rows,'count':len(rows),'bent_checks':294,'elapsed_seconds':time.monotonic()-start,'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
Path(__file__).with_name('mm-cases.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps({k:out[k] for k in ['count','bent_checks','elapsed_seconds']}))
