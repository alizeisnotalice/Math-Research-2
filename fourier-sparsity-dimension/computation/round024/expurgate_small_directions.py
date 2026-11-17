from pathlib import Path
import json,time,hashlib,platform
start=time.monotonic();k=3;q=1<<k;N=q*q;FULL=(1<<N)-1
def truth(fn):return sum(int(fn(x&(q-1),x>>k))<<x for x in range(N))
def basis(vs):
 p={};out=[]
 for original in vs:
  x=original
  while x:
   b=x.bit_length()-1
   if b in p:x^=p[b]
   else:p[b]=x;out.append(original);break
 return out
def span(bs):
 S=[0]
 for b in bs:S += [x^b for x in S]
 return S
def fwht(a):
 a=list(a);h=1
 while h<len(a):
  for i in range(0,len(a),h*2):
   for j in range(i,i+h):a[j],a[j+h]=a[j]+a[j+h],a[j]-a[j+h]
  h*=2
 return a
C=[FULL]+[truth(lambda z,w,j=j:((z|(w<<k))>>j)&1) for j in range(2*k)]
gens=[truth(lambda z,w,a=a:w==a) for a in range(q)]
for j in range(k):
 for a in range(1<<j):gens.append(truth(lambda z,w,j=j,a=a:((z>>j)&1) and (w&((1<<j)-1))==a))
bs=basis(C+gens);assert len(bs)==2*q-1 and bs[:len(C)]==C
quot=span(bs[len(C):]);aff=span(C);bad=[];coset_min=[];base=truth(lambda z,w:(z&w).bit_count()%2)
for c,v in enumerate(quot):
 weights=[(v^a).bit_count() for a in aff];coset_min.append(min(weights))
 if any(0<t<2*q for t in weights):bad.append(c)
 for a in aff:assert all(abs(x)==q for x in fwht([1-2*(((base^v^a)>>z)&1) for z in range(N)]))
assert 0 not in bad
remaining=bad[:];checks=[];steps=[]
while remaining:
 h=min(range(1,len(quot)),key=lambda h:sum((h&c).bit_count()%2==0 for c in remaining))
 new=[c for c in remaining if (h&c).bit_count()%2==0]
 assert len(new)<len(remaining)
 checks.append(h);steps.append([len(remaining),len(new)]);remaining=new
safe=[c for c in range(len(quot)) if all((h&c).bit_count()%2==0 for h in checks)]
nonzero_weights=[(quot[c]^a).bit_count() for c in safe for a in aff if quot[c]^a]
assert min(nonzero_weights)>=2*q
rank=len(basis(checks));assert len(safe)==1<<(len(bs)-len(C)-rank)
out={'k':k,'q':q,'original_dimension':len(bs),'affine_dimension':len(C),'quotient_dimension':len(bs)-len(C),'all_original_bent_phases_checked':len(quot)*len(aff),'bad_coset_count':len(bad),'bad_cosets':bad,'check_rows':checks,'steps':steps,'removed_rank':rank,'retained_dimension':len(bs)-rank,'retained_phase_count':len(safe)*len(aff),'retained_minimum_weight':min(nonzero_weights),'retained_quotient_elements':safe,'original_basis_hex':[hex(x) for x in bs],'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'python':platform.python_version(),'elapsed_seconds':time.monotonic()-start}
Path(__file__).with_suffix('.json').write_text(json.dumps(out,indent=2)+'\n');print({key:v for key,v in out.items() if key not in ['bad_cosets','retained_quotient_elements','original_basis_hex']})
