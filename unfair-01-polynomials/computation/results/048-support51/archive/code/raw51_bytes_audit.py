#!/usr/bin/env python3
"""Single-core raw gzip/hash/row/cache audit for total51, without semantic replay."""
from __future__ import annotations
import gzip, hashlib, json, os, time
from pathlib import Path

ROOT=Path('/Volumes/data/CodexTemp/T2/unfair-01-polynomials/agents/semantic51_audit47')
RAW=Path('/Volumes/data/CodexTemp/T2/unfair-01-polynomials/round042-macbook-support52/run')
EVENT_SHA='87f2152ebfc98df3b0f19daccfc6a23d724a967d218dd6ee6316861617315391'
BINDING={'cap':6,'events':EVENT_SHA,'launch':'7d9e010f5be285ca4c1a19a3386cf2edf224f77918aa49ee4af431625cf74402','payload':'8bcc26b8898f18cb939abb7dc06236686e1847a23b65db711a4cdbcecfc1ef72','plan':'79ba39503b2ef8a10c4a8538a9a77fa48de84bbcf46e154c8d22e4dbd542a7fb','runner':'c9e94bfc715784238309612a98a63bb4619428eac03df41dcb3e1b1a3fdb3719','source':'d68f1010c61ad1d911749742fef8c790c0a8064bf34512b402e7e65e4703e710','terminal':'264bedff460d0e9d745939f2ca796ed71920eaa6cd7942460a7842618c38513c','total':51,'verifier':'a557bf120c6bed03adbd8b2ce153a88bf76e678ed94c6bc8815d14dc33edbe4b'}

def sha_bytes(b:bytes)->str:return hashlib.sha256(b).hexdigest()
def sha_file(p:Path)->str:
 h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(1<<20),b''):h.update(b)
 return h.hexdigest()
def esha(e:dict)->str:return hashlib.sha256(json.dumps(e,sort_keys=True).encode()).hexdigest()
def dump(name:str,o:dict):
 p=ROOT/name;t=p.with_name(p.name+'.partial')
 with t.open('w') as f:json.dump(o,f,sort_keys=True);f.write('\n');f.flush();os.fsync(f.fileno())
 os.replace(t,p)

def main():
 start=time.monotonic(); errors=[]
 input_hash={'raw_events':sha_file(RAW/'events.jsonl'),'raw_payload':sha_file(RAW/'payload.py'),'raw_launch':sha_file(RAW/'launch-manifest.json'),'raw_terminal':sha_file(RAW/'exit-receipt.json'),'copied_events':sha_file(ROOT/'input/run__events.jsonl')}
 expected={'raw_events':EVENT_SHA,'raw_payload':BINDING['payload'],'raw_launch':BINDING['launch'],'raw_terminal':BINDING['terminal'],'copied_events':EVENT_SHA}
 for k in expected:
  if input_hash[k]!=expected[k]: errors.append({'kind':'raw_input_hash_mismatch','key':k,'got':input_hash[k],'expected':expected[k]})
 events=[]; seen=set()
 for no,line in enumerate((ROOT/'input/run__events.jsonl').read_text().splitlines(),1):
  e=json.loads(line)
  if e.get('event')!='chunk' or e.get('total')!=51:continue
  k=(e['m'],e['n'],e['index'])
  if k in seen: errors.append({'kind':'duplicate_event','line':no});continue
  seen.add(k);events.append(e)
 if len(events)!=68405:errors.append({'kind':'event_count','got':len(events)})
 events.sort(key=lambda e:(e['m'],e['index']))
 packed_bytes=raw_bytes=rows=ok=0
 for i,e in enumerate(events,1):
  m,n,index=e['m'],e['n'],e['index']; label=f'm{m:02}n{n:02}/i{index:06d}'
  p=RAW/f't51/m{m:02}n{n:02}/i{index:06d}.jsonl.gz'; c=ROOT/f'cache/t51/m{m:02}/i{index:06}.json'
  try:
   b=p.read_bytes(); packed_bytes+=len(b)
   if sha_bytes(b)!=e.get('gzip_sha256'):raise ValueError('gzip_hash')
   raw=gzip.decompress(b); raw_bytes+=len(raw)
   if sha_bytes(raw)!=e.get('raw_sha256'):raise ValueError('raw_hash')
   line_count=raw.count(b'\n')
   if line_count!=e['counts']['canonical']:raise ValueError(f'line_count:{line_count}')
   r=json.loads(c.read_text())
   if r.get('status')!='passed' or r.get('errors')!=[] or r.get('binding')!=BINDING or r.get('event_sha256')!=esha(e):raise ValueError('cache_binding')
   z=r.get('result',{})
   if z.get('counts')!=e['counts'] or z.get('verified_rows')!=line_count or z.get('raw_support_count')!=e['counts']['raw']:raise ValueError('cache_count')
   ok+=1;rows+=line_count
  except Exception as exc:
   errors.append({'kind':'raw_chunk_or_cache_failure','job':label,'error':repr(exc)})
   break
  if i%1000==0:
   dump('raw51-bytes-checkpoint.json',{'schema':1,'status':'running','checked':i,'expected':len(events),'packed_bytes':packed_bytes,'raw_bytes':raw_bytes,'rows':rows,'errors':errors,'input_hash':input_hash,'elapsed':time.monotonic()-start})
 status='passed' if not errors and ok==68405 else 'rejected'
 report={'schema':1,'scope':'Full byte read of local total51 gzip corpus: gzip SHA, decompressed raw SHA, newline row count, and hash-bound cache record; no row-by-row semantic reclassification.','status':status,'errors':errors,'expected_chunks':len(events),'checked_chunks':ok,'packed_bytes':packed_bytes,'raw_bytes':raw_bytes,'line_rows':rows,'input_hash':input_hash,'elapsed_seconds':time.monotonic()-start}
 dump('raw51-bytes-report.json',report);print(json.dumps({k:report[k] for k in ('status','checked_chunks','packed_bytes','raw_bytes','line_rows','errors')},sort_keys=True))
 if status!='passed':raise SystemExit(2)
if __name__=='__main__':main()
