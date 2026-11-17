#!/usr/bin/env python3
"""Derive immutable coverage-range events from an older native receipt/ledger.

It never edits the original ledger.  It emits only ranges whose receipt says
coverage_credit=true and whose ledger/output bindings agree.
"""
import argparse,hashlib,json
from pathlib import Path
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def main():
 p=argparse.ArgumentParser();p.add_argument('--ledger',type=Path,required=True);p.add_argument('--receipt',type=Path,required=True);p.add_argument('--out',type=Path,required=True);a=p.parse_args()
 r=json.loads(a.receipt.read_text())
 if r.get('status')!='complete' or r.get('coverage_credit') is not True:raise ValueError('receipt is not coverage-complete')
 rows=[]
 for line in a.ledger.read_text().splitlines():
  e=json.loads(line);span=e.get('range');t=e.get('task')
  if e.get('kind')!='streamed_batch_numeric_ranking' or not(isinstance(span,list) and len(span)==2 and type(t)is int):raise ValueError('unexpected original ledger event')
  lo,hi=span
  if r.get('task')!=t or not(r['requested_range'][0]<=lo<hi<=r['requested_range'][1]):raise ValueError('receipt/range disagreement')
  rows.append({'kind':'streamed_batch_numeric_ranking','coverage_credit':True,'task':t,'range':[lo,hi],'count':hi-lo,'input_sha256':e['input_sha256'],'source_sha256':e['source_sha256'],'binary_sha256':e['binary_sha256'],'original_ledger_sha256':sha(a.ledger),'receipt_sha256':sha(a.receipt)})
 a.out.parent.mkdir(parents=True,exist_ok=True);a.out.write_text(''.join(json.dumps(x,sort_keys=True,separators=(',',':'))+'\n' for x in rows));print(json.dumps({'out':str(a.out),'rows':len(rows),'sha256':sha(a.out)},sort_keys=True))
if __name__=='__main__':main()
