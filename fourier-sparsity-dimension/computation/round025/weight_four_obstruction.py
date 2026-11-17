from pathlib import Path
import time,json,hashlib,platform
start=time.monotonic()
def fwht(a):
 a=list(a);h=1
 while h<len(a):
  for i in range(0,len(a),2*h):
   for j in range(i,i+h):a[j],a[j+h]=a[j]+a[j+h],a[j]-a[j+h]
  h*=2
 return a
rows=[]
for k in [2,3,4]:
 q=1<<k;N=q*q;plane=list(range(4))
 basis=[(1<<w)^(1 if w in plane else 0) for w in range(1,q)]
 signs=[]
 for point in range(N):
  z=point&(q-1);w=point>>k
  g=((z&w).bit_count()%2)^int((w>>2)==0 and (z&3)==3)
  signs.append(1-2*g)
 count=1<<(q-1);h=0
 for j in range(count):
  if j:h^=basis[(j&-j).bit_length()-1]
  assert sum((h>>w)&1 for w in plane)%2==0
  vals=[s*(1-2*((h>>(x>>k))&1)) for x,s in enumerate(signs)]
  assert all(abs(v)==q for v in fwht(vals))
 fiber=fwht(signs[:q]);assert sum(x!=0 for x in fiber)==4
 rows.append({'k':k,'q':q,'seed_dimension':q-1,'plane':plane,'compatible_phases_checked':count,'all_bent':True,'z_fiber_w0_spectrum_size':4,'fixed_split_MM':False,'base_phase_hex':hex(sum((s<0)<<i for i,s in enumerate(signs)))})
out={'rows':rows,'arithmetic':'exact F2 and integer Walsh','python':platform.python_version(),'elapsed_seconds':time.monotonic()-start,'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
Path(__file__).with_suffix('.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
