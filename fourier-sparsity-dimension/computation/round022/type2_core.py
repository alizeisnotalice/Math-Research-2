from pathlib import Path
import json,time,hashlib,platform
start=time.monotonic()
B=[x for x in range(64) if ((x&7)==7)^((x>>3)==7)]
assert len(B)==14
def fwht(a):
 a=list(a);h=1
 while h<len(a):
  for i in range(0,len(a),h*2):
   for j in range(i,i+h):a[j],a[j+h]=a[j]+a[j+h],a[j]-a[j+h]
  h*=2
 return a
allowed=[]
for m in range(1<<14):
 a=[0]*64
 for j,x in enumerate(B):a[x]=1-2*((m>>j)&1)
 if set(fwht(a))<={0,-8,8}:allowed.append(m)
aff={sum((c^((u&x).bit_count()&1))<<j for j,x in enumerate(B)) for c in range(2) for u in range(64)}
remaining=set(allowed);orbits=[]
while remaining:
 m=min(remaining);o={m^v for v in aff};assert o<=remaining;orbits.append(m);remaining-=o
out={'core':B,'sign_patterns':1<<14,'allowed_count':len(allowed),'affine_restrictions':len(aff),'orbit_representatives':orbits,'allowed_masks':allowed,'elapsed_seconds':time.monotonic()-start,'python':platform.python_version(),'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
Path(__file__).with_suffix('.json').write_text(json.dumps(out,indent=2)+'\n')
print({k:v for k,v in out.items() if k not in ['allowed_masks','core']})
