from pathlib import Path
import json,hashlib,platform,time
start=time.monotonic();k=3;q=8;N=64;C={0,1,2,3}
def fwht(a):
 a=list(a);h=1
 while h<len(a):
  for i in range(0,len(a),2*h):
   for j in range(i,i+h):a[j],a[j+h]=a[j]+a[j+h],a[j]-a[j+h]
  h*=2
 return a
pi=list(range(q));pi[3],pi[4]=pi[4],pi[3]
signs=[1-2*(((x&(q-1))&pi[x>>k]).bit_count()%2) for x in range(N)]
assert all(abs(x)==q for x in fwht(signs))
partial=fwht([s if (x>>k) in C else 0 for x,s in enumerate(signs)])
D={pi[w] for w in C};actual={x for x,v in enumerate(partial) if v};expected={a|(b<<k) for a in D for b in range(q)}
assert actual==expected and set(partial)<={0,-q,q} and D=={0,1,2,4}
H=[h for h in range(1<<q) if sum((h>>w)&1 for w in C)%2==0];assert len(H)==128
for h in H:assert all(abs(x)==q for x in fwht([s*(1-2*((h>>(x>>k))&1)) for x,s in enumerate(signs)]))
# A four-point affine flat containing zero must be closed under XOR.
assert 1^2 not in D
out={'k':k,'q':q,'component':sorted(C),'permutation':pi,'frequency_component':sorted(D),'frequency_component_affine':False,'all_H_phases_checked':len(H),'all_bent':True,'partial_support_size':len(actual),'partial_values':sorted(set(partial)),'python':platform.python_version(),'elapsed_seconds':time.monotonic()-start,'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
Path(__file__).with_suffix('.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
