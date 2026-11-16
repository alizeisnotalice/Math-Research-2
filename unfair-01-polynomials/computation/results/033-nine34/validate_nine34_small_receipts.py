"""Lightweight supplement binding the already-executed nine34 byte-hash receipt."""
import argparse,hashlib,json
from collections import Counter
from pathlib import Path
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def main():
 p=argparse.ArgumentParser();p.add_argument('--ledger-root',type=Path,required=True);p.add_argument('--mapping',type=Path,required=True);p.add_argument('--old-receipt',type=Path,required=True);p.add_argument('--source-root',type=Path,required=True);p.add_argument('--record',type=Path,required=True);a=p.parse_args()
 old=json.loads(a.old_receipt.read_text());assert sha(a.old_receipt)=='ba8ee877c781ae9fc878c14a186cb83fa88a42a9bed1b93f24bbb02fdb025ccf' and old['unknown_rows']==46
 source={'support_search_layer.py':sha(a.source_root/'support_search_layer.py'),'verify_support_layer.py':sha(a.source_root/'verify_support_layer.py')}
 mapping=json.loads(a.mapping.read_text())['minimal_unknown_inputs'];keys=[(Path(x['source_unknown']['path']).parent.name,x['source_unknown']['line']) for x in mapping];assert len(keys)==len(set(keys))==46
 total=Counter();small=[];actual=[]
 for m in range(1,18):
  d=a.ledger_root/f'nine34-m{m:02d}';s=json.loads((d/'summary.json').read_text());v=json.loads((d/'verification.json').read_text());j=json.loads((d/'job.json').read_text());assert s['script_sha256']==source['support_search_layer.py'] and v['verifier_sha256']==source['verify_support_layer.py'] and j['returncode']==0 and j['status']=='verified' and j['process_exited'] and s['counts']==v['counts']
  c=s['counts'];assert Counter(c)==sum((Counter({'raw':x['raw'],'canonical':x['canonical'],'excluded':x['excluded'],'fair':x['fair'],'unknown':x['unknown']}) for x in s['shards']),Counter()) and c['canonical']==c['excluded']+c['fair']+c['unknown'];total.update(c)
  for line,x in enumerate((d/'unknown.jsonl').read_text().splitlines(),1):actual.append((d.name,line))
  for f in ('summary.json','verification.json','job.json','unknown.jsonl'):small.append({'path':str(d/f),'sha256':sha(d/f)})
 assert set(actual)==set(keys) and len(actual)==46 and total==Counter(old['counts'])
 out={'status':'nine34_small_receipts_and_source_bound_to_prior_full_hash_receipt','prior_full_hash_receipt_sha256':sha(a.old_receipt),'source_hashes':source,'counts':dict(total),'unknown_mapping_unique_and_bound':True,'small_file_whitelist':small,'scope':'No ledger decompression here; binds the prior full-byte receipt rather than replacing it.'};a.record.write_text(json.dumps(out,indent=2)+'\n');print(json.dumps({'status':out['status'],'small_files':len(small)}))
if __name__=='__main__':main()
