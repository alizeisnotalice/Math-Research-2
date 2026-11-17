#!/usr/bin/env python3
"""Validate an already-finished remote stream after transport loss, no rerun."""
import argparse,hashlib,importlib.util,json,os
from pathlib import Path
HERE=Path(__file__).resolve().parent
spec=importlib.util.spec_from_file_location('scheduler',HERE/'native_stream_scheduler.py');m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def main():
 p=argparse.ArgumentParser();p.add_argument('--original',type=Path,required=True);p.add_argument('--journal',type=Path,required=True);p.add_argument('--final',type=Path,required=True);p.add_argument('--out',type=Path,required=True);a=p.parse_args()
 man=json.loads((a.original/'manifest.json').read_text())
 receipt_path=a.original/'receipt.json'
 if receipt_path.exists():
  old=json.loads(receipt_path.read_text());task=old['task'];requested_range=old['requested_range']
  input_sha=old['remote_input_sha256'];source_sha=old['remote_source_sha256'];binary_sha=old['remote_binary_sha256'];old_receipt_sha=sha(receipt_path);receipt_missing=False
 else:
  # A local SSH client may die after persisting its immutable manifest/input,
  # while the remote single-process stream continues.  Never infer coverage
  # from that condition: bind it anew to the locally persisted bytes and the
  # remote final/journal below.
  task=man['task'];requested_range=man['range'];input_sha=man['input_sha256'];source_sha=man['source_sha256'];binary_sha=None;old_receipt_sha=None;receipt_missing=True
  local_input=a.original/f'task{task}-input.json'
  if not local_input.exists() or sha(local_input)!=input_sha: raise ValueError('missing-receipt local input binding')
 doc=json.loads(m.INPUTS.read_text());t=doc['tasks'][task];dims=man['dims'];expected=[tuple(x) for x in man['batches']]
 final=json.loads(a.final.read_text())
 if binary_sha is None: binary_sha=final.get('binary_sha256')
 for k,v in {'command_status':'completed','task':str(task),'input_sha256':input_sha,'source_sha256':source_sha,'binary_sha256':binary_sha}.items():
  if final.get(k)!=v:raise ValueError('final binding '+k)
 events=m.validate_events(a.journal.read_text().splitlines(),expected,t,dims,input_sha,source_sha,binary_sha)
 if len(events)!=len(expected):raise ValueError('incomplete stream')
 a.out.mkdir(parents=True,exist_ok=True);ledger=a.out/'ledger.jsonl'
 with open(ledger,'w') as f:
  for e in events:f.write(json.dumps(e,sort_keys=True,separators=(',',':'))+'\n')
  f.flush();os.fsync(f.fileno())
 receipt={'kind':'recovered-native-stream-receipt','status':'complete','coverage_credit':True,'task':task,'requested_range':requested_range,'verified_batches':len(events),'ledger_sha256':sha(ledger),'remote_journal_sha256':sha(a.journal),'remote_final_sha256':sha(a.final),'original_partial_receipt_sha256':old_receipt_sha,'local_receipt_missing':receipt_missing,'manifest_sha256':sha(a.original/'manifest.json'),'remote_input_sha256':input_sha,'remote_source_sha256':source_sha,'remote_binary_sha256':binary_sha}
 with open(a.out/'receipt.json','w') as f:json.dump(receipt,f,sort_keys=True,indent=2);f.write('\n');f.flush();os.fsync(f.fileno())
 print(json.dumps(receipt,sort_keys=True))
if __name__=='__main__':main()
