import hashlib,json,sys
from pathlib import Path
import numpy as np
B=Path(sys.argv[1]);H=lambda b:hashlib.sha256(b).hexdigest();seed=json.loads((B/'seed.json').read_bytes()); ledger=json.loads((B/'screen/stage-ledger.json').read_bytes());assert ledger['seed_receipt_sha256']==H((B/'seed.json').read_bytes()) if 'seed_receipt_sha256' in ledger else ledger['seed_receipt_sha']==H((B/'seed.json').read_bytes())
prev=seed['remaining'];out=[]
for k, layer in enumerate(ledger['audits'],1):
 assert len(layer['work_files'])==1024; expected=dict(layer['work_files']);tested=accepted=remaining=0
 for part,digest in expected.items():
  file=B/'screen'/f'stage{k}'/f'part{part:04}.npz';rec=json.loads(file.with_suffix('.json').read_text());src=B/('seeds' if k==1 else f'screen/stage{k-1}')/file.name
  assert H(file.read_bytes())==digest==rec['data_sha256'] and H(src.read_bytes())==rec['config']['source_sha256'] and rec['config']['goal']==120696470013529293
  with np.load(file,allow_pickle=False) as x:a=x['accepted'];b=x['remaining'];v=x['full_floor']
  with np.load(src,allow_pickle=False) as x:p=x['pos' if k==1 else 'remaining']
  assert np.array_equal(np.sort(np.r_[a,b]),p) and np.all(v>=120696470013529293) and [len(p),len(a),len(b)]==[rec['tested'],rec['accepted'],rec['remaining']];tested+=len(p);accepted+=len(a);remaining+=len(b)
 assert tested==prev==layer['tested'] and accepted==layer['accepted'] and remaining==layer['remaining'];prev=remaining;out.append(dict(layer=k,parts=1024,tested=tested,accepted=accepted,remaining=remaining,calls=layer['calls']))
(B/'layers-audit.json').write_text(json.dumps(dict(seed=seed['remaining'],layers=out,complete_four_layers=len(out)==4,scope='independent ledger/packet partition and threshold checks, not geometric rerun'),indent=2)+'\n');print(json.dumps(out))
