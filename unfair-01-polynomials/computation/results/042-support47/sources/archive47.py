#!/usr/bin/env python3
"""Compact archive for a passed per-total independent verification report."""
import argparse,gzip,hashlib,json
from collections import Counter,defaultdict
from pathlib import Path
def sha(p):
 h=hashlib.sha256()
 with open(p,'rb') as f:
  for x in iter(lambda:f.read(1<<20),b''):h.update(x)
 return h.hexdigest()
def gz(p,x):
 raw=(json.dumps(x,sort_keys=True,separators=(',',':'))+'\n').encode()
 with gzip.GzipFile(p,'wb',compresslevel=9,mtime=0) as f:f.write(raw)
 return sha(p),hashlib.sha256(raw).hexdigest()
def main():
 ap=argparse.ArgumentParser();ap.add_argument('report',type=Path);ap.add_argument('--out',type=Path,required=True);a=ap.parse_args()
 if a.out.exists():raise RuntimeError('archive exists')
 r=json.loads(a.report.read_text())
 if r.get('status')!='passed_complete_total_47':raise RuntimeError('not passed 47 report')
 a.out.mkdir(parents=True);pairs=defaultdict(Counter);bind=[]
 for x in r['chunk_results']:
  m,n=map(int,(x['job'][1:3],x['job'][4:6]));pairs[m,n]['chunks']+=1;pairs[m,n].update(x['counts']);bind.append({'job':x['job'],'gzip_sha256':x['gzip_sha256'],'raw_sha256':x['raw_sha256'],'counts':x['counts'],'verified_rows':x['verified_rows']})
 full,fullraw=gz(a.out/'total47-full-report.json.gz',r);bgz,braw=gz(a.out/'total47-chunk-bindings.json.gz',bind)
 degree=[{'m':m,'n':n,**dict(pairs[m,n])} for m,n in sorted(pairs)];(a.out/'total47-degree-pairs.json').write_text(json.dumps(degree,indent=2,sort_keys=True)+'\n')
 receipt={'format':'a30-total47-independent-archive-v1','status':r['status'],'report_sha256':sha(a.report),'full_gzip_sha256':full,'full_uncompressed_sha256':fullraw,'bindings_gzip_sha256':bgz,'bindings_uncompressed_sha256':braw,'degree_pairs_sha256':sha(a.out/'total47-degree-pairs.json'),'chunks':len(bind),'counts':r['counts'],'scope':'total47 complete; total48 is not represented or certified here'}
 (a.out/'total47-exit-receipt.json').write_text(json.dumps(receipt,indent=2,sort_keys=True)+'\n');print(json.dumps(receipt,sort_keys=True))
if __name__=='__main__':main()
