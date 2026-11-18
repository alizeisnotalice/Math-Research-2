#!/usr/bin/env python3
"""Extract code-255 rows from a fresh S={2} root-screen output."""
from __future__ import annotations
import argparse,gzip,hashlib,json
from pathlib import Path
def main():
 ap=argparse.ArgumentParser();ap.add_argument('--root',type=Path,required=True);ap.add_argument('--out',type=Path,required=True);a=ap.parse_args()
 s=json.loads((a.root/'summary.json').read_text());assert s['status']=='complete'
 fs=sorted((a.root/'cores').glob('*.json'));assert len(fs)==11760;rows=[];h=hashlib.sha256()
 for p in fs:
  raw=p.read_bytes();h.update(p.name.encode()+b'\0'+raw);x=json.loads(raw);assert x['status']=='complete';rows.extend(r['row'] for r in x['records'] if r['code']==255)
 assert len(rows)==len(set(map(tuple,rows)))==s['counts']['passed'];rows.sort();payload=b''.join(json.dumps(x,separators=(',',':')).encode()+b'\n' for x in rows)
 a.out.parent.mkdir(parents=True,exist_ok=True)
 with a.out.open('wb') as fp:
  with gzip.GzipFile(filename='',mode='wb',fileobj=fp,mtime=0) as g:g.write(payload)
 meta={'status':'complete','mothers':len(rows),'rows_payload_sha256':hashlib.sha256(payload).hexdigest(),'source_summary_sha256':hashlib.sha256((a.root/'summary.json').read_bytes()).hexdigest(),'source_core_ledger_sha256':h.hexdigest(),'source_code_sha256':s['code_sha256'],'gzip_sha256':hashlib.sha256(a.out.read_bytes()).hexdigest()}
 a.out.with_suffix('').with_suffix('.json').write_text(json.dumps(meta,sort_keys=True,separators=(',',':'))+'\n');print(json.dumps(meta,sort_keys=True))
if __name__=='__main__':main()
