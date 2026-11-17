from pathlib import Path
import json,time,hashlib,platform,random
start=time.monotonic()
def rank(vs):
 p={}
 for x in vs:
  while x:
   b=x.bit_length()-1
   if b in p:x^=p[b]
   else:p[b]=x;break
 return len(p)
def fwht(a):
 a=list(a);h=1
 while h<len(a):
  for i in range(0,len(a),2*h):
   for j in range(i,i+h):a[j],a[j+h]=a[j]+a[j+h],a[j]-a[j+h]
  h*=2
 return a
rows=[]
for k in [2,3,4,5]:
 q=1<<k;N=q*q;D=2*q-1
 truth=lambda fn:sum(int(fn(x&(q-1),x>>k))<<x for x in range(N))
 gens=[truth(lambda z,w,u=u:w==u) for u in range(q)]
 for j in range(k):
  for u in range(1<<j):gens.append(truth(lambda z,w,j=j,u=u:((z>>j)&1) and (w&((1<<j)-1))==u))
 assert rank(gens)==D
 hcheck=(1<<q)-1;offset=q+q//2-1;acheck=((1<<(q//2))-1)<<offset
 checks=[hcheck,acheck]
 pivots=[0,offset];kernel=[]
 for j in range(D):
  if j in pivots:continue
  c=1<<j
  for check,pivot in zip(checks,pivots):
   if (check&c).bit_count()%2:c^=1<<pivot
  kernel.append(c)
 assert rank(kernel)==D-2 and all((c&hcheck).bit_count()%2==0 and (c&acheck).bit_count()%2==0 for c in kernel)
 def word(c):
  out=0
  while c:
   b=(c&-c).bit_length()-1;out^=gens[b];c&=c-1
  return out
 base=truth(lambda z,w:(z&w).bit_count()%2)
 low=[];minimum=N
 if k<=3:
  for c in range(1<<D):
   weight=word(c).bit_count()
   if 0<weight<2*q:low.append(c)
   if weight and all((c&h).bit_count()%2==0 for h in checks):minimum=min(minimum,weight)
  assert len(low)==3*q and minimum==2*q
  assert all(any((c&h).bit_count()%2 for h in checks) for c in low)
 rng=random.Random(20261124+k);sample_count=512
 for _ in range(sample_count):
  t=rng.randrange(1<<len(kernel));c=0
  for j,b in enumerate(kernel):
   if t>>j&1:c^=b
  w=word(c);assert not w or w.bit_count()>=2*q
  assert all(abs(x)==q for x in fwht([1-2*(((base^w)>>z)&1) for z in range(N)]))
 assert word(3).bit_count()==2*q and all((3&h).bit_count()%2==0 for h in checks)
 rows.append({'k':k,'q':q,'original_dimension':D,'retained_dimension':D-2,'check_rows_hex':[hex(x) for x in checks],'full_original_weight_enumeration':k<=3,'low_weight_count':len(low) if k<=3 else None,'exact_minimum_weight_by_full_enumeration':minimum if k<=3 else None,'explicit_weight_2q_witness':True,'sampled_bent_phases':sample_count,'seed':20261124+k})
out={'rows':rows,'scope':'full weight enumeration k2,3; 512 fixed-seed bent-phase samples per k2..5; general result separately proved','python':platform.python_version(),'elapsed_seconds':time.monotonic()-start,'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
Path(__file__).with_suffix('.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
