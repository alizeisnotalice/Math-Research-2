#!/usr/bin/env python3
"""One-core independent coverage and cache-binding audit for semantic total 51.

It checks copied bytes, regenerates the cap-six total-51 I table without
trusting the plan, binds every event to one cached audited result, and checks
the successful remote exit.  It does not replay the 68k compressed ledgers;
that stricter semantic replay requires the remote run/t51 byte corpus.
"""
from __future__ import annotations

import hashlib
import itertools
import json
from collections import Counter
from pathlib import Path

ROOT = Path('/Volumes/data/CodexTemp/T2/unfair-01-polynomials/agents/semantic51_audit47')
INPUT = ROOT / 'input'
EXPECTED = {
 'code__audited_chunk_verifier.py':'a557bf120c6bed03adbd8b2ce153a88bf76e678ed94c6bc8815d14dc33edbe4b',
 'code__runner.py':'c9e94bfc715784238309612a98a63bb4619428eac03df41dcb3e1b1a3fdb3719',
 'plan__frozen_support.py':'d68f1010c61ad1d911749742fef8c790c0a8064bf34512b402e7e65e4703e710',
 'plan__support52-plan.json':'79ba39503b2ef8a10c4a8538a9a77fa48de84bbcf46e154c8d22e4dbd542a7fb',
 'run__payload.py':'8bcc26b8898f18cb939abb7dc06236686e1847a23b65db711a4cdbcecfc1ef72',
 'run__events.jsonl':'87f2152ebfc98df3b0f19daccfc6a23d724a967d218dd6ee6316861617315391',
 'run__launch-manifest.json':'7d9e010f5be285ca4c1a19a3386cf2edf224f77918aa49ee4af431625cf74402',
 'run__exit-receipt.json':'264bedff460d0e9d745939f2ca796ed71920eaa6cd7942460a7842618c38513c',
 'launch51.json':'330bf842af54b59ba1f0c4eba1157d259554eccecd61c0dfd416342054eeeeb0',
 'progress51.json':'361cc5de2537d035ee79c41f9a18bb7e073ecdd4eb764aa5c60e0b114a20bdf1',
 'exit51.json':'ebba4dc872f3d65ecbec7a0d5f78ec4295f9d72077f938d52c887ebf523dd09b',
}
CAP=6

def sha(p:Path)->str:
 h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(1<<20),b''):h.update(b)
 return h.hexdigest()
def event_sha(e:dict)->str:return hashlib.sha256(json.dumps(e,sort_keys=True).encode()).hexdigest()
def supports(m:int):
 for size in range(2,min(m+1,CAP)+1):
  for inside in itertools.combinations(range(1,m),size-2):yield (0,*inside,m)
def key(e:dict):return (int(e['total']),int(e['m']),int(e['n']),int(e['index']),tuple(e['I']))

def main():
 errors=[]
 byte_hashes={name:sha(INPUT/name) if (INPUT/name).is_file() else None for name in EXPECTED}
 for name,want in EXPECTED.items():
  if byte_hashes[name]!=want:errors.append({'kind':'input_byte_hash_mismatch','name':name,'got':byte_hashes[name],'expected':want})
 plan=json.loads((INPUT/'plan__support52-plan.json').read_text())
 generated=[]
 for m in range(1,26):
  n=51-m
  generated += [{'total':51,'m':m,'n':n,'index':i,'I':list(I)} for i,I in enumerate(supports(m))]
 listed=[x for x in plan.get('jobs',[]) if x.get('total')==51]
 if listed!=generated:errors.append({'kind':'plan_not_independent_total51_enumeration','listed':len(listed),'generated':len(generated)})
 wanted={key(x):x for x in generated}
 by_m_index={(k[1],k[3]):k for k in wanted}
 if len(wanted)!=68405:errors.append({'kind':'wrong_independent_count','got':len(wanted)})
 events=[];seen=set()
 for no,line in enumerate((INPUT/'run__events.jsonl').read_text().splitlines(),1):
  e=json.loads(line)
  if e.get('event')!='chunk' or e.get('total')!=51:continue
  try:k=key(e)
  except Exception as exc:errors.append({'kind':'malformed_event','line':no,'error':repr(exc)});continue
  if k not in wanted:errors.append({'kind':'event_outside_independent_total51_plan','line':no});continue
  if wanted[k]!= {'total':51,'m':e['m'],'n':e['n'],'index':e['index'],'I':e['I']}:errors.append({'kind':'event_I_index_mismatch','line':no});continue
  if k in seen:errors.append({'kind':'duplicate_event_key','line':no});continue
  seen.add(k);events.append(e)
 if seen!=set(wanted):errors.append({'kind':'event_coverage_incomplete','seen':len(seen),'expected':len(wanted)})
 binding={'cap':6,'events':EXPECTED['run__events.jsonl'],'launch':EXPECTED['run__launch-manifest.json'],
          'payload':EXPECTED['run__payload.py'],'plan':EXPECTED['plan__support52-plan.json'],
          'runner':EXPECTED['code__runner.py'],'source':EXPECTED['plan__frozen_support.py'],
          'terminal':EXPECTED['run__exit-receipt.json'],'total':51,'verifier':EXPECTED['code__audited_chunk_verifier.py']}
 aggregate=Counter(); cache_seen=set(); cache_files=list((ROOT/'cache/t51').rglob('*.json'))
 bykey={key(e):e for e in events}
 for p in cache_files:
  try:r=json.loads(p.read_text()); rel=p.relative_to(ROOT/'cache/t51'); m=int(rel.parts[0][1:]); i=int(rel.stem[1:])
  except Exception as exc:errors.append({'kind':'bad_cache_file','path':str(p),'error':repr(exc)});continue
  k=by_m_index.get((m,i))
  if k is None:errors.append({'kind':'cache_path_outside_plan','path':str(rel)});continue
  e=bykey.get(k)
  if k in cache_seen:errors.append({'kind':'duplicate_cache_key','path':str(rel)});continue
  cache_seen.add(k)
  if e is None:errors.append({'kind':'cache_has_no_event','path':str(rel)});continue
  if r.get('status')!='passed' or r.get('errors')!=[]:errors.append({'kind':'cache_not_passed','path':str(rel)});continue
  if r.get('binding')!=binding:errors.append({'kind':'cache_binding_mismatch','path':str(rel)});continue
  if r.get('event_sha256')!=event_sha(e):errors.append({'kind':'cache_event_hash_mismatch','path':str(rel)});continue
  result=r.get('result',{}); counts=result.get('counts')
  if counts!=e.get('counts') or result.get('raw_support_count')!=e.get('counts',{}).get('raw') or result.get('verified_rows')!=e.get('counts',{}).get('canonical'):
   errors.append({'kind':'cache_result_count_mismatch','path':str(rel)});continue
  aggregate.update(counts)
 if cache_seen!=set(wanted):errors.append({'kind':'cache_coverage_incomplete','seen':len(cache_seen),'expected':len(wanted)})
 progress=json.loads((INPUT/'progress51.json').read_text()); exit51=json.loads((INPUT/'exit51.json').read_text())
 if (progress.get('status'),progress.get('verified'),progress.get('remaining'),progress.get('failures'),progress.get('binding')) != ('complete',68405,0,[],binding):errors.append({'kind':'progress51_terminal_contract_failure'})
 if (exit51.get('returncode'),exit51.get('status'),exit51.get('verified'),exit51.get('remaining'),exit51.get('child_waited')) != (0,'complete',68405,0,True):errors.append({'kind':'exit51_terminal_contract_failure'})
 # Digest uses portable relative names and content hashes.
 rows=[f'{sha(p)}  {p.relative_to(ROOT/"cache/t51").as_posix()}\n' for p in sorted(cache_files)]
 report={'schema':1,'scope':'Byte-bound complete total51 cache/event/coverage audit; no fresh replay of remote gzip ledgers.',
  'status':'passed' if not errors else 'rejected','errors':errors,'input_sha256':byte_hashes,'plan_jobs_independently_generated':len(wanted),
  'event_keys':len(seen),'cache_files':len(cache_files),'cache_keys':len(cache_seen),'cache_relative_content_list_sha256':hashlib.sha256(''.join(rows).encode()).hexdigest(),
  'aggregate_counts_from_cache':dict(aggregate),'progress51':{k:progress.get(k) for k in ('status','verified','remaining','failures','elapsed','workers')},'exit51':exit51}
 (ROOT/'report51.json').write_text(json.dumps(report,indent=2,sort_keys=True)+'\n')
 print(json.dumps({k:report[k] for k in ('status','plan_jobs_independently_generated','event_keys','cache_files','cache_keys','cache_relative_content_list_sha256','errors')},sort_keys=True))
 if errors:raise SystemExit(2)
if __name__=='__main__':main()
