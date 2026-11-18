"""Streaming complete local raw total54 audit; independent of semantic success flags."""
import argparse,concurrent.futures as cf,gzip,hashlib,itertools,json,os,time
from pathlib import Path
from collections import Counter

def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(1<<20),b''):h.update(b)
 return h.hexdigest()
def atom(p,o):
 t=p.with_suffix('.partial');t.write_text(json.dumps(o,sort_keys=True)+'\n');os.replace(t,p)
def field(m,events,run,out,binding):
 p=out/f'm{m:02}.json'
 if p.exists():
  o=json.loads(p.read_text());assert o['binding']==binding and o['status']=='passed';return o
 n=54-m;counts=Counter();gzbytes=rawbytes=0;rows=[]
 for i,e in sorted(events.items()):
  src=Path(e['_run'])/'t54'/f'm{m:02}n{n:02}'/f'i{i:06}.jsonl.gz';assert sha(src)==e['gzip_sha256']
  h=hashlib.sha256();size=lines=0;last=b''
  with gzip.open(src,'rb') as f:
   for b in iter(lambda:f.read(1<<20),b''):h.update(b);size+=len(b);lines+=b.count(b'\n');last=b[-1:]
  assert size==0 or last==b'\n';assert h.hexdigest()==e['raw_sha256'];assert lines==e['counts']['canonical']
  counts.update(e['counts']);gzbytes+=src.stat().st_size;rawbytes+=size;rows.append([i,src.stat().st_size,size,h.hexdigest()])
 o={'status':'passed','m':m,'chunks':len(events),'counts':dict(counts),'gzip_bytes':gzbytes,'raw_bytes':rawbytes,'binding':binding,'chunk_sizes_raw_hashes':rows};atom(p,o);return o

def main():
 ap=argparse.ArgumentParser();ap.add_argument('--run',type=Path,required=True);ap.add_argument('--out',type=Path,required=True);ap.add_argument('--workers',type=int,default=4);a=ap.parse_args();a.out.mkdir(parents=True,exist_ok=True);start=time.monotonic()
 events={m:{} for m in range(1,28)};expected={}
 for m in events:
  expected[m]=[list((0,*inside,m)) for size in range(2,min(m+1,6)+1) for inside in itertools.combinations(range(1,m),size-2)]
 for l in (a.run/'union-events.jsonl').open():
  e=json.loads(l)
  if e.get('event')!='chunk' or e.get('total')!=54:continue
  m,i=e['m'],e['index'];assert m in events and i not in events[m];assert e['n']==54-m and e['I']==expected[m][i];events[m][i]=e
 assert sum(map(len,events.values()))==101583
 assert all(set(events[m])==set(range(len(expected[m]))) for m in events)
 binding={'events':sha(a.run/'union-events.jsonl'),'source':sha(Path(__file__))};done=[]
 with cf.ThreadPoolExecutor(max_workers=a.workers) as ex:
  jobs=[ex.submit(field,m,events[m],a.run,a.out,binding) for m in events]
  for job in cf.as_completed(jobs):
   x=job.result();done.append(x);print(json.dumps({k:x[k] for k in ('m','chunks','status','raw_bytes')}),flush=True)
 counts=Counter()
 for x in done:counts.update(x['counts'])
 assert counts['unknown']==0 and len(done)==27
 o={'status':'passed','chunks':sum(x['chunks'] for x in done),'counts':dict(counts),'gzip_bytes':sum(x['gzip_bytes'] for x in done),'raw_bytes':sum(x['raw_bytes'] for x in done),'binding':binding,'seconds':time.monotonic()-start,'scope':'All planned total54 I supports, every gzip/raw content hash and canonical line count; semantic classification separately verified remotely.'};atom(a.out/'summary.json',o);print(json.dumps(o),flush=True)
if __name__=='__main__':main()
