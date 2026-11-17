import concurrent.futures,gzip,hashlib,json,os,sys
from pathlib import Path
run=Path(sys.argv[1]);out=Path(sys.argv[2]);workers=int(sys.argv[3]);out.mkdir(parents=True,exist_ok=True)
ev={}
for l in (run/'events.jsonl').open():
 x=json.loads(l)
 if x.get('event')=='chunk' and x.get('total')==52:ev[(x['m'],x['index'])]=x
def one(k):
 e=ev[k];p=run/'t52'/f'm{e["m"]:02}n{e["n"]:02}'/f'i{e["index"]:06}.jsonl.gz';b=p.read_bytes();raw=gzip.decompress(b);ok=(hashlib.sha256(b).hexdigest()==e['gzip_sha256'] and hashlib.sha256(raw).hexdigest()==e['raw_sha256'] and (0 if not raw else len(raw.splitlines()))==e['counts']['canonical']);return k,ok,len(b),len(raw)
with concurrent.futures.ThreadPoolExecutor(max_workers=workers) as ex:r=list(ex.map(one,sorted(ev)))
bad=[k for k,ok,_,_ in r if not ok];o={'status':'passed' if not bad else 'rejected','chunks':len(r),'bad':bad[:20],'gzip_bytes':sum(x[2] for x in r),'raw_bytes':sum(x[3] for x in r),'events_sha256':hashlib.sha256((run/'events.jsonl').read_bytes()).hexdigest()};q=out/'raw52.json.partial';q.write_text(json.dumps(o,sort_keys=True)+'\n');os.replace(q,out/'raw52.json');print(json.dumps(o))
