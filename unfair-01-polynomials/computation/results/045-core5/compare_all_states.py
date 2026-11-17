#!/usr/bin/env python3
"""Compare every accepted state and terminal rejection in the two certificates."""
import gzip,hashlib,json
from pathlib import Path
ROOT=Path(__file__).resolve().parent
def load(name):return json.loads(gzip.decompress((ROOT/name).read_bytes()))
def main():
 parent=load('direct-prefix-trees.json.gz')
 second=load('hgq-all-layer-certificates.json.gz')
 summary=load('hgq-summary.json.gz');layers=0
 assert len(parent)==len(summary['results'])==404
 for i,(tree,result) in enumerate(zip(parent,summary['results'])):
  assert tree['core']==result['core'] and result['status']=='dead_exhaustive_prefix'
  assert len(tree['tree'])==len(result['layers'])
  for st,claimed in zip(tree['tree'],result['layers']):
   k=st['k'];observed=second[f'core-{i:03d}/layer-{k:02d}.json.gz']
   expected=[{'G':p['G']+[g],'Q':p['Q']+[q]} for p in st['parents'] for g,q in p['accepted_last']]
   assert observed['accepted']==expected and claimed['accepted']==len(expected)
   assert not claimed['capped'];layers+=1
  last=tree['tree'][-1];death=second[f"core-{i:03d}/death-{last['k']:02d}.json.gz"]
  assert death['previous_states']==[{'G':p['G'],'Q':p['Q']} for p in last['parents']]
  assert death['rejected']==[[j,r['g'],r['q'],r['product_coefficient']] for j,p in enumerate(last['parents']) for r in p['rejected']]
 out={'status':'passed_full_archived_state_comparison','cores':404,'layers':layers,'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
 (ROOT/'archived-state-comparison.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out))
if __name__=='__main__':main()
