from pathlib import Path
import hashlib,json,platform,time
start=time.monotonic()
def fwht(a):
 a=list(a);h=1
 while h<len(a):
  for i in range(0,len(a),2*h):
   for j in range(i,i+h):a[j],a[j+h]=a[j]+a[j+h],a[j]-a[j+h]
  h*=2
 return a
def rank(vs):
 p={}
 for x in vs:
  while x:
   b=x.bit_length()-1
   if b in p:x^=p[b]
   else:p[b]=x;break
 return len(p)
rows=[]
for mu,r in [(3,0),(3,1),(4,0)]:
 k=mu+r;n=2*k;N=1<<n;q=1<<k;mask=(1<<mu)-1;rmask=(1<<r)-1
 f=[];A=[]
 for t in range(N):
  x=t&mask;y=(t>>mu)&mask;z=(t>>(2*mu))&rmask;w=t>>(2*mu+r)
  f.append(((x&y).bit_count()+(z&w).bit_count()+int(x==0 and z==0))%2)
  A.append(int(z==0 and ((x==0) != (y==0))))
 patches=[set(),set()]
 for xi in range(N):
  x=xi&mask;y=(xi>>mu)&mask;w=xi>>(2*mu+r)
  if w==0 and x==0 and (y&1):patches[0].add(xi)
  if w==0 and y==0 and (x&1):patches[1].add(xi)
 assert not patches[0]&patches[1] and all(len(S)==q//2 for S in patches)
 normal_rows=[]
 for S in patches:
  b=min(S);directions=[x^b for x in S];assert rank(directions)==k-1
  normals=[z for z in range(N) if all((z&d).bit_count()%2==0 for d in directions)]
  assert rank(normals)==k+1;normal_rows+=normals
 # Parameter order: n linear coefficients, constant, indicator A. Each
 # spectral patch translates by the linear coefficient only.
 pr=rank(normal_rows);assert pr==k+mu
 K=[z for z in range(1<<(n+2)) if all((z&u).bit_count()%2==0 for u in normal_rows)]
 restrictions=[];VR=[]
 for par in K:
  lin=par&(N-1);c=(par>>n)&1;t=par>>(n+1)
  restrictions.append(sum((((lin&x).bit_count()%2)^c^(t*A[x]))<<x for x in range(N) if A[x]))
 assert rank(restrictions)==1 and len(K)==1<<(r+2)
 phases=0
 for lin in range(N):
  for c in range(2):
   for t in range(2):
    signs=[1-2*(f[x]^((lin&x).bit_count()%2)^c^(t*A[x])) for x in range(N)]
    assert all(abs(z)==q for z in fwht(signs))
    partial=fwht([signs[x]*A[x] for x in range(N)])
    assert set(partial)<={0,-q,q}
    actual={xi for xi,z in enumerate(partial) if z}
    for S in patches:assert {xi^lin for xi in S}<=actual
    phases+=1
 rows.append({'mu':mu,'free_dimension':r,'k':k,'phase_count':phases,'all_bent':True,'patch_sizes':[len(S) for S in patches],'joint_patch_codimensions':[k+1,k+1],'parameter_projection_rank':pr,'K_dimension':r+2,'K_restriction_rank':rank(restrictions)})
out={'rows':rows,'scope':'real MM common family RM1+one cross indicator; partial spectral patches','python':platform.python_version(),'elapsed_seconds':time.monotonic()-start,'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
Path(__file__).with_suffix('.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
