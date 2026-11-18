"""Raise .83750 whole-leaf floors using archived, pre-existing whole-leaf certificates."""
import hashlib,io,json,sys
from pathlib import Path
from zipfile import ZipFile
import numpy as np
C,B,O=map(Path,sys.argv[1:]); O.mkdir(exist_ok=False);h=lambda b:hashlib.sha256(b).hexdigest()
report=json.loads((B/'base-report.json').read_bytes());assert report['canonical_leaves']==6710314 and report['max_index_target']=='67/80'; goal=120696470013529293
positions={};floors={}
for row in report['rows']:
 part=row[0];z=np.load(B/'base-index'/f'part{part:04}.npz',allow_pickle=False);p=z['pos'];v=z['full_floor'];assert [h(p.tobytes()),h(v.tobytes())]==row[-2:];positions[part]=p;floors[part]=v.copy()
counts=[]; inputs={}
for code in ('32','35','40','41','42','43','44'):
 path=C/f'tomar-0837{code}-stage-ledger.json';blob=path.read_bytes();ledger=json.loads(blob);assert ledger['complete'] and len(ledger['audits'])==4 and ledger['parts']==1024; inputs[path.name]=h(blob)
 wanted={(k,part):digest for k,layer in enumerate(ledger['audits'],1) for part,digest in layer['work_files']}; assert len(wanted)==4096
 seen=set(); updates=0; before=sum(int(np.count_nonzero(v<goal)) for v in floors.values());total_accepted=0
 for archive in sorted(C.glob(f'tomar-scored-0837{code}*.zip')):
  inputs[archive.name]=h(archive.read_bytes())
  with ZipFile(archive) as z:
   for name in z.namelist():
    if not name.endswith('.npz'): continue
    stem=Path(name);k=int(stem.parts[0].removeprefix('stage'));part=int(stem.stem.removeprefix('part'));key=(k,part);data=z.read(name);assert key in wanted and key not in seen and h(data)==wanted[key];seen.add(key)
    with np.load(io.BytesIO(data),allow_pickle=False) as f: accepted=f['accepted']; lower=f['full_floor']
    assert accepted.dtype==np.int32 and lower.dtype==np.int64 and accepted.shape==lower.shape and np.all(lower>=ledger['integer_goal']);total_accepted+=len(accepted)
    p=positions[part];ix=np.searchsorted(p,accepted);mask=ix<len(p);safe=ix[mask];a=accepted[mask];eq=p[safe]==a;safe=safe[eq];vals=lower[mask][eq];old=floors[part][safe];updates+=int(np.count_nonzero(vals>old));floors[part][safe]=np.maximum(old,vals)
 assert seen==set(wanted);after=sum(int(np.count_nonzero(v<goal)) for v in floors.values());counts.append(dict(stage=code,accepted=total_accepted,improvements=updates,before=before,after=after,excluded=before-after));print(counts[-1],flush=True)
rows=[]
for part,p in sorted(positions.items()):
 f=floors[part]; keep=f<goal;p=p[keep];f=f[keep];np.savez_compressed(O/f'part{part:04}.npz',pos=p,full_floor=f);rows.append([part,len(p),h(p.tobytes()),h(f.tobytes())])
out=dict(goal=goal,canonical_leaves=6710314,base_report_sha256=h((B/'base-report.json').read_bytes()),base=1202260,stages=counts,remaining=sum(x[1] for x in rows),rows=rows,inputs_sha256=inputs,scope='maxima of prior whole-original-leaf integer lower certificates; old geometry inherited, not re-executed')
(B/'merged-report.json').write_text(json.dumps(out,indent=2)+'\n')
