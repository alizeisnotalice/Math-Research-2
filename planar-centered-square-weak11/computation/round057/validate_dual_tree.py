#!/usr/bin/env python3
"""Direct-merge check of R57's frozen dual dynamic union tree."""
import json,random,hashlib,sys,os
from pathlib import Path
from full_area_sided_sweep import DualUnionTree
p=Path(__file__).resolve().parent; rng=random.Random(570021);N=31;D=1<<32
ys=[i*i*D for i in range(N)];dy=[rng.randrange(-100,101)*D for i in range(N)]
tree=DualUnionTree([(v,v) for v in ys],[(v,v) for v in dy]);iv=[(i,j) for i in range(N) for j in range(i+1,N)];active=set()
for step in range(5000):
 k=rng.randrange(len(iv));a,b=iv[k]
 if k in active:active.remove(k);tree.update(a,b,-1)
 else:active.add(k);tree.update(a,b,1)
 unions=[]
 for a,b in sorted(iv[j] for j in active):
  if unions and a<=unions[-1][1]:unions[-1]=(unions[-1][0],max(unions[-1][1],b))
  else:unions.append((a,b))
 length=sum(ys[b]-ys[a] for a,b in unions);derivative=sum(dy[b]-dy[a] for a,b in unions)
 assert tree.result()==((length,length),(derivative,derivative))
out={'status':'passed','scope':'5000 dual-tree insert/delete states match independent sorted direct interval merge for both length and derivative','seed':570021,'coordinate_count':N,'update_count':5000,'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'sweep_source_sha256':hashlib.sha256((p/'full_area_sided_sweep.py').read_bytes()).hexdigest(),'python':sys.version}
assert os.path.ismount('/Volumes/data');(p/'dual-tree-validation.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
