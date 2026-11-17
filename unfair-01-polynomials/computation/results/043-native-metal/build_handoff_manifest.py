#!/usr/bin/env python3
"""Build a strict, prefix-only GPU handoff manifest from durable event streams.

Non-JSON progress banners are counted and ignored.  A JSON `batch` event that
is malformed is fatal; it is never silently skipped.  The resulting manifest
binds the raw SHA-256 of every named v1/v2/v3 stream and gives one exact
contiguous next boundary for every task.
"""
import argparse, hashlib, json
from pathlib import Path

INPUTS=Path('/Volumes/data/CodexTemp/T2/unfair-01-polynomials/round038-macbook-gpu52-direct/inputs.json')
def sha(p): return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def main():
 p=argparse.ArgumentParser();p.add_argument('--layer',action='append',required=True,help='name=event.jsonl');p.add_argument('--out',type=Path,required=True);a=p.parse_args()
 d=json.loads(INPUTS.read_text()); intervals={i:[] for i in range(len(d['tasks']))};layers=[]
 for spec in a.layer:
  if spec.count('=')!=1:p.error('layer must be name=path')
  name,raw=spec.split('=');path=Path(raw);nonjson=0;events=0
  skips=[];native_receipt=None
  for line_no,line in enumerate(path.read_text(errors='replace').splitlines(),1):
   try: event=json.loads(line)
   except json.JSONDecodeError:
    stripped=line.strip()
    if not stripped or stripped.startswith(('Checking GPU compatibility','GPU compatibility','Wolfram Language','GPUArray::','Set::','MessageTemplate::')): nonjson+=1;continue
    if '"event"' in line or '"batch"' in line:raise ValueError(f'{path}:{line_no}: truncated or malformed event-looking JSON')
    raise ValueError(f'{path}:{line_no}: unrecognised non-JSON line')
   if not isinstance(event,dict):raise ValueError(f'{path}:{line_no}: JSON line not object')
   typ=event.get('event') or event.get('kind')
   if typ=='batch':
    t,lo,hi,count=(event.get(k) for k in ('task','first','last','count'))
    if type(t)is not int or type(lo)is not int or type(hi)is not int or type(count)is not int or not(0<=t<len(d['tasks'])) or lo<0 or count<=0 or hi-lo+1!=count:raise ValueError(f'{path}:{line_no}: malformed batch')
    intervals[t].append((lo,hi+1,name,line_no));events+=1
   elif typ=='streamed_batch_numeric_ranking':
    t,span=event.get('task'),event.get('range')
    count=event.get('count')
    if event.get('coverage_credit') is not True or type(t)is not int or not(isinstance(span,list) and len(span)==2 and type(span[0])is int and type(span[1])is int) or type(count)is not int:raise ValueError(f'{path}:{line_no}: malformed native range event')
    lo,hi=span
    if not(0<=t<len(d['tasks'])) or lo<0 or hi<=lo or hi>d['tasks'][t]['masks'] or count!=hi-lo or any(not(isinstance(event.get(k),str) and len(event[k])==64) for k in ('input_sha256','source_sha256','binary_sha256')):raise ValueError(f'{path}:{line_no}: native range bounds/binding')
    intervals[t].append((lo,hi,name,line_no));events+=1
    receipt_path=path.parent/'receipt.json'
    if native_receipt is None:
     native_receipt=json.loads(receipt_path.read_text())
     if native_receipt.get('status')!='complete' or native_receipt.get('coverage_credit') is not True:raise ValueError(f'{path}:{line_no}: native receipt not coverage-complete')
    if native_receipt.get('task')!=t or not (native_receipt['requested_range'][0]<=lo and hi<=native_receipt['requested_range'][1]):raise ValueError(f'{path}:{line_no}: native receipt binding')
    if any(native_receipt.get('remote_'+k)!=event.get(k) for k in ('input_sha256','source_sha256','binary_sha256')):raise ValueError(f'{path}:{line_no}: native receipt hash binding')
    if 'original_ledger_sha256' in event:
     original=path.parent/'ledger.jsonl'
     if event.get('original_ledger_sha256')!=sha(original):raise ValueError(f'{path}:{line_no}: original ledger hash binding')
     if event.get('receipt_sha256')!=sha(receipt_path):raise ValueError(f'{path}:{line_no}: receipt hash binding')
    elif native_receipt.get('ledger_sha256')!=sha(path):raise ValueError(f'{path}:{line_no}: native ledger hash binding')
   elif typ=='task_skipped_completed':
    t,start,masks=(event.get(k) for k in ('task','start','masks'))
    if type(t)is not int or type(start)is not int or type(masks)is not int or not(0<=t<len(d['tasks'])) or start!=masks or masks!=d['tasks'][t]['masks']:raise ValueError(f'{path}:{line_no}: malformed completed skip')
    skips.append((t,masks,line_no));events+=1
  layers.append({'name':name,'path':str(path),'sha256':sha(path),'recognized_events':events,'non_json_progress_lines':nonjson,'completed_skip_assertions':skips,'native_receipt_sha256':sha(receipt_path) if native_receipt is not None else None})
 nexts={}
 for t,rows in intervals.items():
  end=0
  for lo,hi,name,line in sorted(rows):
   if lo<end: raise ValueError(f'task {t}: overlap at {name}:{line}')
   if lo!=end: raise ValueError(f'task {t}: gap [{end},{lo}) before {name}:{line}')
   end=hi
  if end>d['tasks'][t]['masks']:raise ValueError(f'task {t}: exceeds mask universe')
  for layer in layers:
   for st,sm,line in layer['completed_skip_assertions']:
    if st==t and end<sm:raise ValueError(f'task {t}: skip assertion {layer["name"]}:{line} lacks actual batch coverage')
  nexts[str(t)]=end
 out={'kind':'a30-gpu52-gap-handoff-v1','source_inputs_sha256':sha(INPUTS),'layers':layers,'task_next':nexts,'rule':'future run must start exactly at task_next; this is a coverage boundary, not a mathematical exclusion'}
 a.out.parent.mkdir(parents=True,exist_ok=True);a.out.write_text(json.dumps(out,sort_keys=True,indent=2)+'\n');print(json.dumps({'out':str(a.out),'sha256':sha(a.out),'task_next':nexts},sort_keys=True))
if __name__=='__main__':main()
