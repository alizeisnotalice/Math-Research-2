import hashlib,io,json,sys,tarfile
from pathlib import Path
import numpy as np
B,S,O=map(Path,sys.argv[1:]);O.mkdir(exist_ok=False);h=lambda b:hashlib.sha256(b).hexdigest();r=json.loads((B/'merged-report.json').read_text());assert r['remaining']==753913;goal=120696470013529293
ps={};fs={}
for part,count,hp,hf in r['rows']:
 with np.load(B/'merged-index-v2'/f'part{part:04}.npz',allow_pickle=False) as z:p=z['pos'];f=z['full_floor']
 assert len(p)==count and [h(p.tobytes()),h(f.tobytes())]==[hp,hf];ps[part]=p;fs[part]=f.copy()
inputs={};stages=[]
for k in range(1,5):
 path=S/f'native45-screen-{k}.tar.gz';inputs[path.name]=h(path.read_bytes());before=sum(int(np.count_nonzero(f<goal)) for f in fs.values());seen=set();acc=0
 with tarfile.open(path) as t:
  mf=json.load(t.extractfile('MANIFEST.json'));d={name:(size,sha) for name,size,sha in mf};assert len(d)==2048 and len(mf)==2048
  for m in t:
   if not m.name.endswith('.npz'):continue
   assert m.name in d and m.name.startswith(f'screen/stage{k}/part'); data=t.extractfile(m).read();assert (len(data),h(data))==d[m.name];part=int(Path(m.name).stem[4:]);assert part not in seen;seen.add(part)
   with np.load(io.BytesIO(data),allow_pickle=False) as z:a=z['accepted'];lower=z['full_floor']
   assert a.dtype==np.int32 and lower.dtype==np.int64 and a.shape==lower.shape and np.all(lower>=120689264254125501);acc+=len(a);p=ps[part];ix=np.searchsorted(p,a);valid=ix<len(p);ix=ix[valid];aa=a[valid];sel=p[ix]==aa;ix=ix[sel];fs[part][ix]=np.maximum(fs[part][ix],lower[valid][sel])
 assert len(seen)==1024;after=sum(int(np.count_nonzero(f<goal)) for f in fs.values());stages.append(dict(layer=k,accepted=acc,before=before,after=after,excluded=before-after));print(stages[-1],flush=True)
rows=[]
for part,p in sorted(ps.items()):
 f=fs[part];keep=f<goal;p=p[keep];f=f[keep];np.savez_compressed(O/f'part{part:04}.npz',pos=p,full_floor=f);rows.append([part,len(p),h(p.tobytes()),h(f.tobytes())])
out=dict(goal=goal,canonical_leaves=6710314,prior_remaining=r['remaining'],remaining=sum(a[1] for a in rows),stages=stages,rows=rows,volumes_sha256=inputs,prior_report_sha256=h((B/'merged-report.json').read_bytes()),scope='inherited .83745 whole-leaf integer scores; geometry is not rerun by this merger')
(B/'merged45-report.json').write_text(json.dumps(out,indent=2)+'\n')
