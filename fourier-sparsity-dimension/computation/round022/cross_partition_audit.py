from pathlib import Path
import hashlib,json,time,platform
start=time.monotonic()
def rank(v):
 p={}
 for a in v:
  while a:
   k=a.bit_length()-1
   if k in p:a^=p[k]
   else:p[k]=a;break
 return len(p)
rows=[]
for mu in range(3,7):
 N=1<<mu
 pieces=[]
 for arm in range(2):
  for j in range(mu):
   pieces.append({((1<<j)|x)<<(arm*mu) for x in range(1<<j)})
 target=set(range(1,N))|{x<<mu for x in range(1,N)}
 union=set();dims=[];point_pieces=0
 for S in pieces:
  assert not union&S;union|=S
  b=min(S);d=rank([x^b for x in S]);assert len(S)==1<<d;dims.append(d)
  point_pieces+=d==0
 assert union==target and point_pieces==2
 # Every translation ell in the intrinsic 2mu-space is checked. Each
 # translated partition is obtained independently from its affine masks.
 checks=0
 for ell in range(N*N):
  actual={x^ell for x in target}
  predicted=set()
  for arm in range(2):
   for j in range(mu):
    base=ell^((1<<j)<<(arm*mu))
    predicted.update(base^(x<<(arm*mu)) for x in range(1<<j))
  assert actual==predicted;checks+=1
 projected_normals=[]
 for arm in range(2):
  for j in range(mu):
   free={arm*mu+t for t in range(j)}
   projected_normals.extend(1<<t for t in range(2*mu) if t not in free)
 projection_rank=rank(projected_normals);assert projection_rank==2*mu
 rows.append({'mu':mu,'pieces':len(pieces),'dimensions':dims,'point_pieces':point_pieces,'translations_checked':checks,'normal_projection_rank_for_full_translation_parameter':projection_rank,'constant_parameter_kernel_dimension':2*mu+1-projection_rank})
out={'scope':'local cross partition; not a global bent family','rows':rows,'elapsed_seconds':time.monotonic()-start,'python':platform.python_version(),'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
Path(__file__).with_suffix('.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
