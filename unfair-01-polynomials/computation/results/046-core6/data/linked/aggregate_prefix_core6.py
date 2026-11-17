#!/Users/zzh/.venvs/main/bin/python
import gzip, hashlib, json
from collections import Counter
from pathlib import Path
H=Path(__file__).resolve().parent
parts=(H/'prefix_tree_0_3708.jsonl.gz',H/'prefix_tree_3708_7417.jsonl.gz')
c=Counter(); depth=Counter(); total=0
for p in parts:
 with gzip.open(p,'rt') as f:
  for line in f:
   r=json.loads(line); total+=1; c[r['status']]+=1; depth[len(r['levels'])]+=1
assert total==7417 and c==Counter({'DEAD':7417})
out={'status':'COMPLETE_PREFIX_EXHAUSTION','cores':total,'statuses':dict(c),'death_depth_distribution':dict(sorted(depth.items())), 'tree_sha256':{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in parts}, 'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(), 'scope':'Exact bounded-prefix-tree exhaustion conditional on the prior finite core screens and endpoint filters.'}
q=H/'core6_prefix_complete_summary.json';q.write_text(json.dumps(out,indent=2,sort_keys=True)+'\n');print(json.dumps(out,indent=2));print('summary_sha256='+hashlib.sha256(q.read_bytes()).hexdigest())
