"""Freeze full returned total54 semantic receipts and disjoint provenance losslessly."""
import gzip,hashlib,json,os,sys,time,shutil
from pathlib import Path
B=Path(__file__).resolve().parent;T=B.parent;S=T/'sync';R=S/'unfair-01-polynomials/computation/results/056-support54';remote=B/'support54-remote'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
 assert os.path.ismount('/Volumes/data');assert not R.exists();start=time.monotonic();a=R/'archive';(a/'input').mkdir(parents=True);(a/'cache').mkdir();(R/'audits').mkdir()
 inputs=[]
 sources={p.name:p for p in (remote/'input').iterdir() if p.is_file()}
 sources.update({n:remote/n for n in ['progress54.json','exit54.json','launch54.json']})
 sources.update({'runner.py':T/'round053-parent/support54_runner.py','audited_chunk_verifier.py':remote/'code/audited_chunk_verifier.py'})
 for name,p in sorted(sources.items()):
  raw=p.read_bytes();q=a/'input'/f'{name}.gz';q.write_bytes(gzip.compress(raw,mtime=0));assert q.stat().st_size<10*1024**2
  inputs.append({'name':name,'path':str(q.relative_to(a)),'bytes':q.stat().st_size,'sha256':sha(q),'raw_sha256':hashlib.sha256(raw).hexdigest()})
 caches=[];total=0
 for m in range(1,28):
  files=sorted((remote/f'cache/t54/m{m:02}').glob('i*.json'));assert files;q=a/'cache'/f'm{m:02}.jsonl.gz'
  with q.open('wb') as f:
   with gzip.GzipFile(fileobj=f,mode='wb',mtime=0,filename='') as z:
    for p in files:
     raw=p.read_bytes();z.write((json.dumps({'relativepath':str(p.relative_to(remote/'cache/t54')),'sha256':hashlib.sha256(raw).hexdigest(),'utf8':raw.decode()},separators=(',',':'))+'\n').encode())
  assert q.stat().st_size<10*1024**2;caches.append({'path':str(q.relative_to(a)),'bytes':q.stat().st_size,'sha256':sha(q),'records':len(files)});total+=len(files)
 assert total==101583
 (a/'package-manifest.json').write_text(json.dumps({'schema':'support54-semantic-lossless-v1','inputs':inputs,'caches':caches,'coverage':{'total':54,'cap':6,'chunks':total,'old':10784,'new':90799},'builder_sha256':sha(Path(__file__)),'python':sys.version,'seconds':time.monotonic()-start},indent=2)+'\n')
 for src,name in [(B/'raw54-audit/summary.json','raw54.json'),(B/'semantic54-audit.json','semantic54.json'),(B/'semantic54-audit-exit.json','semantic54-exit.json')]:shutil.copyfile(src,R/'audits'/name)
 for src in [Path(__file__),B/'audit_raw54.py',B/'audit_semantic54.py']:shutil.copyfile(src,R/src.name)
 print(json.dumps({'chunks':total,'bytes':sum(x['bytes'] for x in inputs+caches),'seconds':time.monotonic()-start}))
if __name__=='__main__':main()
