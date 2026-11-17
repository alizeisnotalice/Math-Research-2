"""Lossless semantic-cache archive; reuse already public identical inputs."""
import argparse,gzip,hashlib,json,os,sys,time
from pathlib import Path
def sha(b):return hashlib.sha256(b).hexdigest()
def main():
    a=argparse.ArgumentParser();a.add_argument('--cache',type=Path,required=True);a.add_argument('--receipts',type=Path,required=True);a.add_argument('--shared',type=Path,required=True);a.add_argument('--out',type=Path,required=True);q=a.parse_args();r=q.out;r.mkdir(exist_ok=True,parents=True);(r/'cache').mkdir(exist_ok=True);(r/'input').mkdir(exist_ok=True);begin=time.monotonic()
    m=json.loads((q.shared/'package-manifest.json').read_text());shared=[]
    for x in m['input_archives']:
        if x['source_name'] in ('exit51.json','progress51.json','launch51.json'):continue
        b=(q.shared/x['file']).read_bytes();assert sha(b)==x['sha256'] and sha(gzip.decompress(b))==x['source_sha256'];shared.append(x)
    inputs=[]
    for name in ['exit52.json','progress52.json']:
        raw=(q.receipts/name).read_bytes();packed=gzip.compress(raw,mtime=0);p=r/'input'/f'{name}.gz';p.write_bytes(packed);inputs.append({'file':str(p.relative_to(r)),'bytes':len(packed),'sha256':sha(packed),'source_name':name,'source_sha256':sha(raw)})
    caches=[];total=0
    for degree in range(1,27):
        folder=q.cache/f'm{degree:02}';files=sorted(folder.glob('i*.json'));assert files
        p=r/'cache'/f'cache-m{degree:02}.jsonl.gz'
        with p.open('wb') as target:
            with gzip.GzipFile(fileobj=target,mode='wb',mtime=0,filename='') as z:
                for f in files:
                    raw=f.read_bytes();record={'relativepath':str(f.relative_to(q.cache)),'sha256':sha(raw),'utf8':raw.decode('utf8')};z.write((json.dumps(record,separators=(',',':'))+'\n').encode())
        b=p.read_bytes();assert len(b)<10*1024**2;caches.append({'file':str(p.relative_to(r)),'bytes':len(b),'sha256':sha(b),'records':len(files)});total+=len(files)
    assert total==83681
    rows=[{'path':str(f.relative_to(r)),'bytes':f.stat().st_size,'sha256':sha(f.read_bytes())} for d in ('input','cache') for f in sorted((r/d).iterdir())]
    doc={'schema':'support52-semantic-lossless-v1','scope':'Archive of actual independent semantic receipts, not another semantic replay. Parent verifier independently rebuilds all I supports and every binding.','shared_package':'../../048-support51/archive','shared_input_archives':shared,'input_archives':inputs,'files':rows,'coverage':{'total':52,'cap':6,'semantic_cache_records':total,'cache_archives':caches},'builder_sha256':sha(Path(__file__).read_bytes()),'python':sys.version,'elapsed_seconds':time.monotonic()-begin}
    (r/'package-manifest.json').write_text(json.dumps(doc,indent=2)+'\n');print(json.dumps({'records':total,'archive_bytes':sum(x['bytes'] for x in rows),'seconds':doc['elapsed_seconds']}))
if __name__=='__main__':main()
