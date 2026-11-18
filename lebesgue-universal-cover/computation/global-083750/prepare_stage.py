import json, hashlib, shutil
from pathlib import Path
from zipfile import ZipFile
import sys
c=Path(sys.argv[1]);d=Path(sys.argv[2]);d.mkdir(parents=True,exist_ok=False)
index=json.loads((c/'tomar-083731-score-archive-index.json').read_bytes());sha=lambda b:hashlib.sha256(b).hexdigest();seen=set()
for a in index['archives']:
 raw=(c/a['name']).read_bytes();assert sha(raw)==a['sha256']
 with ZipFile(c/a['name']) as z:
  assert len(z.namelist())==a['members']
  for name in z.namelist():
   p=Path(name);assert len(p.parts)==2 and p.parts[0] in {f'stage{i}' for i in range(1,5)} and p.name.startswith('part') and p.suffix in ('.npz','.json') and name not in seen
   seen.add(name);(d/p.parent).mkdir(exist_ok=True);(d/p).write_bytes(z.read(name))
assert len(seen)==8192
for src,dest,k in [('tomar-083731-stage-ledger.json','stage-ledger.json','stage_ledger_sha256'),('tomar-083731-retained.json','retained.json','retained_sha256')]:
 b=(c/src).read_bytes();assert sha(b)==index[k];(d/dest).write_bytes(b)
print('verified stage members',len(seen))
