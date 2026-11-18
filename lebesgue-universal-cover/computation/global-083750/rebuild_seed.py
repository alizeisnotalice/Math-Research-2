"""Rebuild the seed from the original tree, Git inputs and frozen 0.83745 layers."""
import hashlib,json,subprocess,sys,zipfile
from pathlib import Path
P=Path(__file__).resolve().parent;C,T,CACHE,W=map(lambda x:Path(x).resolve(),sys.argv[1:]);W.mkdir(parents=True,exist_ok=False);py=sys.executable
run=lambda script,args:subprocess.run([py,'-B',str(script),*map(str,args)],check=True)
run(P/'prepare_stage.py',[C,W/'stages31'])
with (W/'base-report.json').open('wb') as o:subprocess.run([py,'-B',str(C/'audit_tomar_multi_target_scores.py'),str(C),str(T),str(W/'stages31'),str(W/'build'),str(CACHE),'--emit',str(W/'base-index'),'--targets','.83745','.83750'],check=True,stdout=o)
run(P/'merge_old_scores.py',[C,W,W/'merged-index-v2']);run(P/'merge_45.py',[W,C/'six-input-083745/volumes',W/'merged45-index']);run(P/'prepare_seed.py',[W,C])
with zipfile.ZipFile(P/'volumes/seed50.zip') as z:
 rows=json.loads(z.read('MANIFEST.json'));assert len(rows)==1026
 for name,size,digest in rows:
  old=z.read(name);new=(W/name).read_bytes();assert len(old)==size and hashlib.sha256(old).hexdigest()==digest and new==old,name
print(json.dumps(dict(complete=True,parts=1024,remaining=699060,seed_sha256=hashlib.sha256((W/'seed.json').read_bytes()).hexdigest(),scope='historical whole-leaf area kernels are inherited; not rerun here')))
