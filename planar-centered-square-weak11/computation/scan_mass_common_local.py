#!/usr/bin/env python3
"""Frozen 41x41 local rational grid in the R20 common-mass family."""
import argparse,json,os
from fractions import Fraction as F
from itertools import product
from pathlib import Path
from periodic_halo import sha
from scan_supercell_mass import digest
from scan_mass_response import run

def main():
 p=argparse.ArgumentParser();p.add_argument('--base',type=Path,required=True);p.add_argument('--target',type=Path,required=True);p.add_argument('--binary',type=Path,required=True);p.add_argument('--previous',type=Path,required=True);p.add_argument('--out',type=Path,required=True);p.add_argument('--workers',type=int,default=12);args=p.parse_args()
 if not os.path.ismount('/Volumes/data'):raise RuntimeError('data not mounted')
 base=json.loads(args.base.read_text());target=json.loads(args.target.read_text());bm=json.loads((args.binary/'manifest.json').read_text());kernel=sha(Path(__file__).with_name('periodic_halo_dominated.py'))
 if bm['base_sha256']!=sha(args.base) or bm['target_sha256']!=sha(args.target) or bm['dominated_kernel_sha256']!=kernel:raise ValueError('input/kernel mismatch')
 vectors=[[F(x) for x in v['vector']] for v in bm['selected_vectors']];orig=[F(p[2]) for p in target['pts']];upper=[F(1001,1000)*F(base['pts'][i%len(base['pts'])][2]) for i in range(len(orig))];cache={}
 for path in (args.previous/'batch0',args.previous/'batch1'):
  m=json.loads((path/'manifest.json').read_text());rr={r['id']:r for r in json.loads((path/'results.json').read_text()) if 'error' not in r}
  if m['base_sha256']!=sha(args.base) or m['target_sha256']!=sha(args.target) or m['dominated_kernel_sha256']!=kernel:raise ValueError('cache mismatch')
  for t in m['tasks']:
   r=rr.get(t['id'])
   if r and r['input_sha256']==t['input_sha256'] and digest(t['input'])==t['input_sha256']:cache[t['input_sha256']]=(t,r,sha(path/'results.json'))
 points=[];excluded=[]
 for i,j in product(range(-20,21),repeat=2):
  a=F(-3,20)+F(i,200);b=F(33,4)+F(j,40);c=(a,a,F(1),F(1),b,b,b,b);vals=[w+sum((k*v[q] for k,v in zip(c,vectors)),F(0)) for q,w in enumerate(orig)]
  if min(vals)<=0 or any(w>u for w,u in zip(vals,upper)) or sum(vals)!=sum(orig):excluded.append({'a':str(a),'b':str(b)});continue
  points.append(c)
 args.out.mkdir(parents=True,exist_ok=True);provenance={'orchestrator':'scan_mass_common_local.py','orchestrator_sha256':sha(__file__),'shared_batch_worker':'scan_mass_response.py','shared_batch_worker_sha256':sha(Path(__file__).with_name('scan_mass_response.py')),'source_binary_sha256':sha(args.binary/'manifest.json'),'previous_coarse_sha256':sha(args.previous/'batch0'/'results.json'),'previous_fine_sha256':sha(args.previous/'batch1'/'results.json')};(args.out/'orchestrator_provenance.json').write_text(json.dumps(provenance,indent=2)+'\n')
 selection={'grid':'a=-3/20+i/200,b=33/4+j/40, i,j=-20..20','center':['-3/20','33/4'],'excluded':excluded,'kappa_cap':'1001/1000','orchestrator_sha256':sha(__file__)}
 run(base,target,args.base,args.target,args.out,args.workers,0,points,vectors,cache,selection)
if __name__=='__main__':main()
