"""Integrity/coverage binder for completed nine34 receipts; not a reclassification."""
import argparse,gzip,hashlib,json,math,time
from collections import Counter
from pathlib import Path
def sha(p):
 h=hashlib.sha256()
 with Path(p).open('rb') as f:
  for b in iter(lambda:f.read(1<<20),b''):h.update(b)
 return h.hexdigest()
def main():
 p=argparse.ArgumentParser();p.add_argument('--ledger-root',type=Path,required=True);p.add_argument('--mapping',type=Path,required=True);p.add_argument('--record',type=Path,required=True);a=p.parse_args();t=time.monotonic();mapping=json.loads(a.mapping.read_text());mapped=mapping['minimal_unknown_inputs'];by={}
 for x in mapped:by[(Path(x['source_unknown']['path']).parent.name,x['source_unknown']['line'])]=(x['I'],x['J'])
 seen=set();tot=Counter();unknown=[];rows=[]
 for m in range(1,18):
  d=a.ledger_root/f'nine34-m{m:02d}';s=json.loads((d/'summary.json').read_text());v=json.loads((d/'verification.json').read_text());j=json.loads((d/'job.json').read_text());c=s['config'];assert c['total']==34 and c['min_total']==33 and c['max_terms']==9 and c['min_degree']==c['max_first_degree']==m and c['version']==3 and s['complete'] and v['verified'] and j['status']=='verified' and j['process_exited'] and s['counts']==v['counts']
  pairs=[(x['m'],x['n']) for x in s['shards']];assert all(x not in seen for x in pairs);seen.update(pairs)
  for x in s['shards']:
   # Both endpoints are mandatory, so at most seven interior support points.
   number=lambda q:sum(math.comb(q-1,k) for k in range(min(7,q-1)+1));assert x['raw']==number(x['m'])*number(x['n']) and x['complete']
  h=hashlib.sha256();n=0
  with gzip.open(d/'ledger.jsonl.gz','rb') as f:
   for b in iter(lambda:f.read(1<<20),b''):h.update(b);n+=len(b)
  assert h.hexdigest()==s['ledger_sha256_uncompressed']==v['ledger_sha256_uncompressed']
  for line_no,line in enumerate((d/'unknown.jsonl').read_text().splitlines(),1):
   z=json.loads(line);unknown.append((d.name,line_no,z['I'],z['J']))
  tot.update(s['counts']);rows.append({'shard':d.name,'gzip_sha256':sha(d/'ledger.jsonl.gz'),'uncompressed_sha256':h.hexdigest(),'uncompressed_bytes':n,'counts':s['counts']})
 expected={(m,33-m) for m in range(1,17)}|{(m,34-m) for m in range(1,18)};assert seen==expected and len(seen)==33 and tot['unknown']==len(unknown)==46
 assert set((d,l) for d,l,_,_ in unknown)==set(by) and all(by[d,l]==(I,J) for d,l,I,J in unknown)
 out={'status':'nine34_completed_receipts_bound_to_ledger_bytes','degree_pairs':33,'shards':17,'counts':dict(tot),'unknown_rows':46,'unknown_mapping_bound':True,'ledgers':rows,'mapping_sha256':sha(a.mapping),'elapsed_seconds':time.monotonic()-t,'scope':'One sequential decompression/hash pass plus receipt/coverage binding; no per-record mathematical rerun and no claim about pending external certificate.'};a.record.write_text(json.dumps(out,indent=2)+'\n');print(json.dumps({k:out[k] for k in ('status','counts','elapsed_seconds')}))
if __name__=='__main__':main()
