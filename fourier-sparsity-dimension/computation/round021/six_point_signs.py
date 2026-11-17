from pathlib import Path
import json,hashlib,time
start=time.monotonic()
B=[x for x in range(16) if (((x&1)*((x>>1)&1))^(((x>>2)&1)*((x>>3)&1)))]
def fwht(a):
 a=list(a);h=1
 while h<len(a):
  for i in range(0,len(a),2*h):
   for j in range(i,i+h):a[j],a[j+h]=a[j]+a[j+h],a[j]-a[j+h]
  h*=2
 return a
allowed=[]
for mask in range(64):
 h=[0]*16
 for j,x in enumerate(B):h[x]=1-2*((mask>>j)&1)
 w=fwht(h)
 if set(w)<={0,4,-4}:allowed.append({'mask':mask,'spectrum':w,'support':[i for i,z in enumerate(w) if z]})
A=set(r['mask'] for r in allowed)
translations={sum((((a&x).bit_count()&1)^c)<<j for j,x in enumerate(B)) for a in range(16) for c in range(2)}
assert len(B)==6 and len(translations)==32
out={'B':B,'sign_patterns_tested':64,'allowed_count':len(allowed),'affine_restriction_count':len(translations),'allowed':allowed,'single_affine_orbit':bool(A) and {min(A)^v for v in translations}==A,'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'elapsed_seconds':time.monotonic()-start}
Path(__file__).with_name('six-point-signs.json').write_text(json.dumps(out,indent=2)+'\n')
print({k:out[k] for k in ['B','allowed_count','affine_restriction_count','single_affine_orbit']})
