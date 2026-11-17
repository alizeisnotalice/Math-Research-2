#!/usr/bin/env python3
"""Build and independently re-read a compact public staging package for total51."""
from __future__ import annotations
import gzip, hashlib, io, json, os, sys
from pathlib import Path

ROOT=Path('/Volumes/data/CodexTemp/T2/unfair-01-polynomials/agents/semantic51_audit47')
OUT=ROOT/'package48'
CACHE=ROOT/'cache/t51'
INPUT=ROOT/'input'
MAX=10*1024*1024
EXPECTED_REPORT='274cd08c68a9f89258cc520dced90c0597e91ff95448e5c570abb3043c3370ae'
EXPECTED_RAW='f20c2b98dd856a6fa952f5d57feb10d11f409471071f27fa2cde35e1a679560f'

def sha(p:Path)->str:
 h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(1<<20),b''):h.update(b)
 return h.hexdigest()
def atomic_bytes(p:Path,b:bytes):
 p.parent.mkdir(parents=True,exist_ok=True);t=p.with_name(p.name+'.partial')
 with t.open('wb') as f:f.write(b);f.flush();os.fsync(f.fileno())
 os.replace(t,p)
def gzip_bytes(target:Path,data:bytes):
 b=io.BytesIO()
 with gzip.GzipFile(fileobj=b,mode='wb',compresslevel=9,mtime=0) as z:z.write(data)
 atomic_bytes(target,b.getvalue())

def cache_degree(mdir:Path):
 files=sorted(mdir.glob('i*.json'))
 records=[]; lines=[]
 for p in files:
  text=p.read_text(encoding='utf-8')
  rel='t51/'+p.relative_to(CACHE).as_posix()
  records.append((rel,sha(p),text))
  lines.append(json.dumps({'relativepath':rel,'sha256':records[-1][1],'utf8':text},ensure_ascii=False,separators=(',',':'))+'\n')
 target=OUT/'cache'/f'cache-{mdir.name}.jsonl.gz'
 gzip_bytes(target,''.join(lines).encode('utf-8'))
 # Independent round-trip: every UTF-8 payload and declared hash is checked
 # against the source cache before admitting this archive.
 seen=[]
 with gzip.open(target,'rt',encoding='utf-8') as z:
  for line in z:
   r=json.loads(line); seen.append((r['relativepath'],r['sha256'],r['utf8']))
 if seen!=records:raise RuntimeError(f'roundtrip mismatch {mdir.name}')
 if target.stat().st_size>=MAX:raise RuntimeError(f'archive reaches 10MiB limit: {target}')
 return {'file':target.relative_to(OUT).as_posix(),'records':len(records),'sha256':sha(target),'bytes':target.stat().st_size}

def main():
 if OUT.exists():raise SystemExit(f'refuse existing staging {OUT}')
 if sha(ROOT/'report51.json')!=EXPECTED_REPORT or sha(ROOT/'raw51-bytes-report.json')!=EXPECTED_RAW:raise SystemExit('audits are not the expected completed reports')
 expected_files=[INPUT/x for x in sorted(p.name for p in INPUT.iterdir() if p.is_file())]
 if len(expected_files)!=11:raise SystemExit(f'input file count {len(expected_files)}')
 cache_dirs=sorted(p for p in CACHE.iterdir() if p.is_dir())
 if [p.name for p in cache_dirs]!=[f'm{x:02}' for x in range(1,26)]:raise SystemExit('unexpected cache degree directories')
 OUT.mkdir(); (OUT/'cache').mkdir(); (OUT/'input').mkdir(); (OUT/'code').mkdir()
 cache_entries=[cache_degree(d) for d in cache_dirs]
 if sum(x['records'] for x in cache_entries)!=68405:raise SystemExit('cache coverage count')
 input_entries=[]
 for p in expected_files:
  target=OUT/'input'/(p.name+'.gz');gzip_bytes(target,p.read_bytes())
  if gzip.decompress(target.read_bytes())!=p.read_bytes():raise RuntimeError(f'input roundtrip {p.name}')
  if target.stat().st_size>=MAX:raise RuntimeError(f'input archive reaches 10MiB: {p.name}')
  input_entries.append({'file':target.relative_to(OUT).as_posix(),'source_name':p.name,'source_sha256':sha(p),'sha256':sha(target),'bytes':target.stat().st_size})
 for name in ('audit51.py','raw51_bytes_audit.py','build_package48.py','report51.json','raw51-bytes-report.json','REPRODUCE.md','transfer-integrity.json'):
  src=ROOT/name;dst=OUT/'code'/name;atomic_bytes(dst,src.read_bytes())
 # This concrete public guide is deliberately staged, not a pointer to raw corpus.
 readme='''# Total-51 semantic cache evidence (round 48 staging)\n\nThis package covers the finite cap=6 total-degree 51 search only.  It does not\nclaim a general theorem.  `cache/cache-mXX.jsonl.gz` stores one lossless record\nper original cache JSON: `relativepath`, its SHA-256, and its original UTF-8\nbytes.  Recreate any cache file by writing `utf8.encode("utf-8")` to\n`relativepath`; verify its SHA-256 first.\n\nAll listed archive files are below 10 MiB.  `input/` contains gzip-compressed\nexact input bytes.  `code/` contains the independent coverage/cache audit and\nthe full raw-gzip byte audit reports.  The raw `t51` gzip corpus is deliberately\nnot included because it is reproducible and 3.7 GiB; `raw51-bytes-report.json`\nrecords that it was fully read locally.\n\nTo validate staged cache archives, decompress each JSONL, verify every record\nSHA after UTF-8 reconstruction, and compare all files to `package-manifest.json`.\nThe existing `code/REPRODUCE.md` gives the original audit commands and scope.\n'''
 atomic_bytes(OUT/'README.md',readme.encode())
 files=[p for p in sorted(OUT.rglob('*')) if p.is_file()]
 manifest={'schema':1,'scope':'Finite cap-six total51 public staging package; no general conclusion and no raw corpus.',
  'python':{'executable':sys.executable,'version':sys.version},
  'coverage':{'total':51,'cap':6,'semantic_cache_records':68405,'cache_archives':cache_entries},
  'bindings':{'semantic_runner_sha256':'c9e94bfc715784238309612a98a63bb4619428eac03df41dcb3e1b1a3fdb3719','frozen_support_sha256':'d68f1010c61ad1d911749742fef8c790c0a8064bf34512b402e7e65e4703e710','audited_verifier_sha256':'a557bf120c6bed03adbd8b2ce153a88bf76e678ed94c6bc8815d14dc33edbe4b','exit51_sha256':'ebba4dc872f3d65ecbec7a0d5f78ec4295f9d72077f938d52c887ebf523dd09b','progress51_sha256':'361cc5de2537d035ee79c41f9a18bb7e073ecdd4eb764aa5c60e0b114a20bdf1'},
  'input_archives':input_entries,
  'raw_corpus_included':False,'raw_corpus_local_byte_audit_report_sha256':EXPECTED_RAW,
  'files':[{'path':p.relative_to(OUT).as_posix(),'sha256':sha(p),'bytes':p.stat().st_size} for p in files]}
 atomic_bytes(OUT/'package-manifest.json',(json.dumps(manifest,ensure_ascii=False,indent=2,sort_keys=True)+'\n').encode())
 atomic_bytes(OUT/'package-manifest.sha256',(sha(OUT/'package-manifest.json')+'  package-manifest.json\n').encode())
 print(json.dumps({'status':'passed','cache_archives':len(cache_entries),'cache_records':sum(x['records'] for x in cache_entries),'input_archives':len(input_entries),'manifest_sha256':sha(OUT/'package-manifest.json')},sort_keys=True))
if __name__=='__main__':main()
