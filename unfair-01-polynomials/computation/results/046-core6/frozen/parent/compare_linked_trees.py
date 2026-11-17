#!/usr/bin/env python3
"""Compare every accepted state in direct-AB and independent H*G*Q trees."""
import argparse,gzip,json,hashlib,itertools
from pathlib import Path

def rows(p):
 with gzip.open(p,'rt') as f:
  for l in f: yield json.loads(l)
def main():
 ap=argparse.ArgumentParser();ap.add_argument('direct',type=Path);ap.add_argument('hgq',nargs='+',type=Path);ap.add_argument('--out',type=Path,required=True);a=ap.parse_args()
 cores=layers=0;seen=set()
 other={}
 for y in itertools.chain.from_iterable(rows(p) for p in a.hgq):
  r=y['row'];key=tuple(r[k] for k in ('b','T','N','ell','d','e'));assert key not in other;other[key]=y
 for x in rows(a.direct):
  y=other.pop(tuple(x['core']))
  r=y['row'];core=tuple(r[k] for k in ('b','T','N','ell','d','e'))
  assert tuple(x['core'])==core and core not in seen;seen.add(core)
  assert y['status']=='DEAD' and len(x['levels'])==len(y['levels'])
  for u,v in zip(x['levels'],y['levels']):
   assert u['k']==v['k'];live=[];rejected=0
   for p in u['parents']:
    live.extend((tuple(p['G']+[g]),tuple(p['Q']+[q])) for g,q in p['accepted'])
    rejected+=sum(z[2]=='C_nonbinary' for z in p['rejected'])
   assert len(live)==len(set(live))
   assert sorted(live)==sorted((tuple(g),tuple(q)) for g,q in v['live'])
   assert rejected==v['rejected']['binary_C'];layers+=1
  assert not y['levels'][-1]['live'];cores+=1
 assert not other and cores==7417 and layers==55149
 out={'status':'passed_all_accepted_states_and_binary_rejections','cores':cores,'layers':layers,'bindings':{str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in [a.direct,*a.hgq]}}
 a.out.write_text(json.dumps(out,indent=2)+'\n');print(json.dumps({'cores':cores,'layers':layers,'status':out['status']}))
if __name__=='__main__':main()
