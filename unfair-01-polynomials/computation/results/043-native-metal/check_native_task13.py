#!/usr/bin/env python3
"""Independently replay saved masks from roots; numerical only, not exclusion."""
import gzip,hashlib,json,time
from pathlib import Path
import reconstruct_verify as v
D=Path(__file__).resolve().parent
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
 start=time.monotonic();doc=json.loads(gzip.decompress((D/'inputs.json.gz').read_bytes()));task=v.task_from_source(doc,13)
 end=11272192;rows=[];batches=0
 for name in ['round043-native-task13-first','round043-native-task13-rest']:
  d=D/name;r=json.loads((d/'receipt.json').read_text());assert r['status']=='complete' and r['coverage_credit'] is True
  assert sha(d/'ledger.jsonl')==r['ledger_sha256'];assert sha(d/'remote-stream.jsonl')==r['remote_journal_sha256'];assert sha(d/'task13-input.json')==r['remote_input_sha256'];assert sha(D/'NativeScoreProduction.swift')==r['remote_source_sha256']
  count=0
  for line in (d/'ledger.jsonl').read_text().splitlines():
   e=json.loads(line);lo,hi=e['range'];assert e['task']==13 and lo==end and hi<=task['masks'];end=hi
   for k in ['input_sha256','source_sha256','binary_sha256']:assert e[k]==r['remote_'+k]
   masks=[x['mask'] for x in e['candidates_cpu_reconstructed']];assert len(set(masks))==len(masks)==4
   for mask in masks:
    assert lo<=mask<hi;z=v.reconstruct(task,mask,1e-8,1e-8,1e-5)
    rows.append({k:z[k] for k in ['mask','classification','min_coefficient','max_product_target_error','max_imaginary_residual']})
   count+=1
  assert count==r['verified_batches'] and end==r['verified_prefix_end'];batches+=count
 assert end==33554432
 out=dict(status='completed_parent_independent_numerical_replay',range=[11272192,end],masks=end-11272192,batches=batches,saved_candidates=len(rows),classifications={c:sum(z['classification']==c for z in rows) for c in sorted({z['classification'] for z in rows})},max_product_error=max(z['max_product_target_error'] for z in rows),rows=rows,seconds=time.monotonic()-start,source_sha256=sha(Path(__file__)),checker_sha256=sha(D/'reconstruct_verify.py'),scope='Only the saved candidates are independently reconstructed; all float roots, no exact certificate and no general exclusion.')
 (D/'parent-task13-replay.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps({k:z for k,z in out.items() if k!='rows'}))
if __name__=='__main__':main()
