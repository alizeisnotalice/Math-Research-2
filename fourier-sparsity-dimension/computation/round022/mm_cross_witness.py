from pathlib import Path
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
for mu,r in [(3,0),(3,1),(3,2),(4,0),(4,1),(5,0)]:
 k=mu+r;N=1<<(2*k);q=1<<k;mask=(1<<mu)-1;rmask=(1<<r)-1
 f=[];g=[];A=[]
 for t in range(N):
  x=t&mask;y=(t>>mu)&mask;z=(t>>(2*mu))&rmask;w=t>>(2*mu+r)
  base=((x&y).bit_count()+(z&w).bit_count())%2
  a=base^int(x==0 and z==0);b=base^int(y==0 and z==0)
  f.append(1-2*a);g.append(1-2*b)
  if a!=b:A.append(t)
 assert all(abs(v)==q for v in fwht(f)) and all(abs(v)==q for v in fwht(g))
 assert len(A)==2*((1<<mu)-1)*(1<<r)
 rows.append({'mu':mu,'free_dimension':r,'k':k,'support_size':len(A),'q':q,'both_bent':True})
out={'rows':rows,'python':platform.python_version(),'arithmetic':'integer FWHT; exact','elapsed_seconds':time.monotonic()-start,'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
Path(__file__).with_suffix('.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
