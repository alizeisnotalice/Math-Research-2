"""Build self-contained lossless total53 semantic receipts; no raw ledgers copied."""
import gzip,hashlib,json,os,sys,time
from pathlib import Path
T=Path('/Volumes/data/CodexTemp/T2/unfair-01-polynomials');B=Path(__file__).parent;R=T/'sync/unfair-01-polynomials/computation/results/053-support53';archive=R/'archive';old=T/'round043-macbook-support54-corrected/run';transfer=T/'round050-parent/support53';remote=B/'support53-remote'
def sha(b):return hashlib.sha256(b).hexdigest()
def main():
 assert os.path.ismount('/Volumes/data');start=time.monotonic();assert not R.exists();(archive/'input').mkdir(parents=True);(archive/'cache').mkdir();(R/'audits').mkdir()
 sources={'run__events.jsonl':old/'events.jsonl','run__launch-manifest.json':old/'launch-manifest.json','run__payload.py':old/'payload.py','run__exit-receipt.json':old/'exit-receipt.json','plan__support54-plan.json':T/'agents/support54_macbook43/plan/support54-plan.json','plan__frozen_support.py':T/'agents/support54_macbook43/plan/frozen_support.py','code__runner.py':transfer/'support53_batch_runner.py','code__audited_chunk_verifier.py':transfer/'code__audited_chunk_verifier.py',**{n:remote/n for n in ['progress53.json','exit53.json','launch53.json']}}
 inputs=[]
 for name,p in sources.items():
  raw=p.read_bytes();data=gzip.compress(raw,mtime=0);q=archive/'input'/f'{name}.gz';q.write_bytes(data);assert len(data)<10*1024**2;inputs.append({'file':str(q.relative_to(archive)),'bytes':len(data),'sha256':sha(data),'source_name':name,'source_sha256':sha(raw)})
 caches=[];total=0
 for m in range(1,27):
  folder=remote/f'cache/t53/m{m:02}';files=sorted(folder.glob('i*.json'));assert files;q=archive/'cache'/f'cache-m{m:02}.jsonl.gz'
  with q.open('wb') as f:
   with gzip.GzipFile(fileobj=f,mode='wb',mtime=0,filename='') as z:
    for p in files:
     raw=p.read_bytes();z.write((json.dumps({'relativepath':str(p.relative_to(remote/'cache/t53')),'sha256':sha(raw),'utf8':raw.decode()},separators=(',',':'))+'\n').encode())
  data=q.read_bytes();assert len(data)<10*1024**2;caches.append({'file':str(q.relative_to(archive)),'bytes':len(data),'sha256':sha(data),'records':len(files)});total+=len(files)
 assert total==83681
 files=[{'path':str(p.relative_to(archive)),'bytes':p.stat().st_size,'sha256':sha(p.read_bytes())} for folder in ['input','cache'] for p in sorted((archive/folder).iterdir())]
 manifest={'schema':'support53-semantic-lossless-v1','scope':'Actual independent semantic receipts plus provenance. Archive audit is not another semantic replay.','shared_package':'.','shared_input_archives':[],'input_archives':inputs,'files':files,'coverage':{'total':53,'cap':6,'semantic_cache_records':total,'cache_archives':caches},'builder_sha256':sha(Path(__file__).read_bytes()),'python':sys.version,'seconds':time.monotonic()-start};(archive/'package-manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
 for src,name in [(B/'raw53-audit/summary.json','raw53.json'),(B/'support53-complete-audit.json','semantic53-parent.json'),(remote/'exit53.json','exit53.json')]: (R/'audits'/name).write_bytes(src.read_bytes())
 for src in [Path(__file__),B/'audit_raw53.py',B/'audit_semantic53.py']:(R/src.name).write_bytes(src.read_bytes())
 print(json.dumps({'records':total,'archive_bytes':sum(f['bytes'] for f in files),'seconds':manifest['seconds'],'package':str(R)}))
if __name__=='__main__':main()
