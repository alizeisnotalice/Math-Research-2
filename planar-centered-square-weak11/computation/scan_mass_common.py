#!/usr/bin/env python3
"""Two common amplitudes with exact positivity and kappa cap."""
import argparse,json,os
from fractions import Fraction as F
from itertools import product
from pathlib import Path
from periodic_halo import sha
from scan_supercell_mass import canon,digest
from scan_mass_response import run

def main():
 p=argparse.ArgumentParser();p.add_argument('--base',type=Path,required=True);p.add_argument('--target',type=Path,required=True);p.add_argument('--binary',type=Path,required=True);p.add_argument('--response',type=Path,required=True);p.add_argument('--out',type=Path,required=True);p.add_argument('--workers',type=int,default=12);args=p.parse_args()
 if not os.path.ismount('/Volumes/data'):raise RuntimeError('data not mounted')
 base=json.loads(args.base.read_text());target=json.loads(args.target.read_text());bm=json.loads((args.binary/'manifest.json').read_text());kernel=sha(Path(__file__).with_name('periodic_halo_dominated.py'))
 if bm['base_sha256']!=sha(args.base) or bm['target_sha256']!=sha(args.target) or bm['dominated_kernel_sha256']!=kernel:raise ValueError('input/kernel mismatch')
 vectors=[[F(x) for x in v['vector']] for v in bm['selected_vectors']];orig=[F(p[2]) for p in target['pts']];upper=[F(1001,1000)*F(base['pts'][i%len(base['pts'])][2]) for i in range(len(orig))];cache={}
 for path in (args.binary,args.response/'batch0',args.response/'batch1'):
  m=json.loads((path/'manifest.json').read_text());rr={r['id']:r for r in json.loads((path/'results.json').read_text()) if 'error' not in r}
  if m['base_sha256']!=sha(args.base) or m['target_sha256']!=sha(args.target) or m['dominated_kernel_sha256']!=kernel:raise ValueError('incompatible cache')
  for t in m['tasks']:
   r=rr.get(t['id'])
   if r and r['input_sha256']==t['input_sha256'] and digest(t['input'])==t['input_sha256']:cache[t['input_sha256']]=(t,r,sha(path/'results.json'))
 def convert(a,b):return (a,a,F(1),F(1),b,b,b,b)
 def viable(c):
  vals=[w+sum((k*v[i] for k,v in zip(c,vectors)),F(0)) for i,w in enumerate(orig)]
  return min(vals)>0 and all(w<=u for w,u in zip(vals,upper)) and sum(vals)==sum(orig)
 def prepare(ab):
  yes=[];no=[]
  for a,b in ab:
   c=convert(a,b)
   if viable(c):yes.append(c)
   else:no.append({'a':str(a),'b':str(b),'reason':'nonpositive mass or kappa above 1001/1000'})
  return yes,no
 center=(F(1,2),F(3));coarse=[(F(i,8),F(5*j,8)) for i,j in product(range(-8,9),range(17))];coarse.append(center);points,rejected=prepare(coarse)
 tasks,res,top=run(base,target,args.base,args.target,args.out,args.workers,0,points,vectors,cache,{'grid':'a=-1..1 step1/8, b=0..10 step5/8 plus R18 baseline','excluded':rejected,'kappa_cap':'1001/1000','source_binary_sha256':sha(args.binary/'manifest.json')});byid={t['id']:t for t in tasks};rr={r['id']:r for r in res};winner=top[0];c=tuple(map(F,winner['coefficients'][0]));print(f'coarse winner a={c[0]} b={c[4]}',flush=True)
 for t in tasks:cache[t['input_sha256']]=(t,rr[t['id']],sha(args.out/'batch0'/'results.json'))
 fine=[(F(i,20),F(j,4)) for i,j in product(range(-20,21),range(41))];points,rejected=prepare(fine)
 run(base,target,args.base,args.target,args.out,args.workers,1,points,vectors,cache,{'grid':'a=-1..1 step1/20, b=0..10 step1/4','excluded':rejected,'kappa_cap':'1001/1000','coarse_winner':[str(c[0]),str(c[4])]})
if __name__=='__main__':main()
