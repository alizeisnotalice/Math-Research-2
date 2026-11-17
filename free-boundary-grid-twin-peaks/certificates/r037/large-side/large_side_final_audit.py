#!/usr/bin/env python3
"""Read-only final coverage audit: old labelled l=1..3 plus row-orbit l=4,5."""
from __future__ import annotations
import argparse, hashlib, json, os, sys, tempfile
from pathlib import Path
N=11
class Bad(RuntimeError):pass
def load(p):
 try:x=json.loads(Path(p).read_text())
 except Exception as e:raise Bad(f'cannot parse {p}: {e}')
 if not isinstance(x,dict):raise Bad(f'not object {p}')
 return x
def sha(p):
 h=hashlib.sha256()
 with open(p,'rb') as f:
  for b in iter(lambda:f.read(1<<20),b''):h.update(b)
 return h.hexdigest()
def comb(n,k):
 import math
 return math.comb(n,k)
def orbit_total(l):return comb((1<<(N-l))+l-1,l)
def check_receipts(prod,sid):
 a=prod/'logs'/(sid+'.spawn.json');b=prod/'logs'/(sid+'.wait.json');c=prod/'logs'/(sid+'.reap.json')
 if not(a.is_file() and b.is_file() and c.is_file()):raise Bad(f'{sid}: missing callback receipt')
 x,y,z=load(a),load(b),load(c);pid=x.get('pid')
 if not isinstance(pid,int) or y.get('pid')!=pid or z.get('pid')!=pid or y.get('returncode')!=0 or z.get('returncode')!=0:raise Bad(f'{sid}: bad callback receipt')
def intervals(items,total):
 cur=0;gap=[];over=[]
 for a,b in sorted(items):
  if a<cur:over.append([a,min(b,cur)])
  elif a>cur:gap.append([cur,a])
  cur=max(cur,b)
 if cur<total:gap.append([cur,total])
 return {'total':total,'covered':sum(b-a for a,b in items),'intervals':len(items),'gaps':gap,'overlaps':over,'complete':not gap and not over and cur==total}
def write_new(p,x):
 if p.exists():raise Bad(f'refusing overwrite {p}')
 p.parent.mkdir(parents=True,exist_ok=True);q=p.with_name('.'+p.name+'.tmp')
 with q.open('w')as f:json.dump(x,f,indent=2,sort_keys=True);f.write('\n')
 os.replace(q,p)
def main():
 ap=argparse.ArgumentParser();ap.add_argument('--old-manifest',type=Path,required=True);ap.add_argument('--row-production',type=Path,required=True);ap.add_argument('--output',type=Path,required=True);ap.add_argument('--require-complete',action='store_true');z=ap.parse_args()
 try:
  old=load(z.old_manifest)
  if old.get('status')!='complete_l123_authenticated' or old.get('N')!=N:raise Bad('old manifest is not complete authenticated l123')
  oldcov={l:[] for l in range(1,4)}; olddigest=hashlib.sha256()
  for r in old.get('records',[]):
   raw=r.get('interval');p=Path(r.get('path',''))
   if not(isinstance(raw,list) and len(raw)==3 and raw[0] in (1,2,3) and p.is_file() and r.get('sha256')==sha(p)):raise Bad('bad old manifest record/hash')
   l,a,b=raw
   if not(0<=a<b<=1<<(l*(N-l))):raise Bad('old interval bounds')
   oldcov[l].append((a,b));olddigest.update(('%s\t%s\n'%(raw,r['sha256'])).encode())
  if old.get('records_digest')!=olddigest.hexdigest():raise Bad('old record digest mismatch')
  meta=load(z.row_production/'started.json'); source=meta.get('source_sha256')
  if meta.get('N')!=N or meta.get('sorted_transposed_left_sizes')!=[7,6] or meta.get('original_small_left_by_sorted_transposed_left')!={'7':4,'6':5} or not isinstance(source,str):raise Bad('bad large-side production metadata')
  rowcov={4:[],5:[]}; weights={4:0,5:0}; noncomplete={}
  for p in sorted((z.row_production/'results').glob('*.json')):
   d=load(p);l=d.get('left_size');s=d.get('status')
   if l not in (7,6):raise Bad(f'wrong transposed-row result family {p}')
   if s!='completed':noncomplete[str(s)]=noncomplete.get(str(s),0)+1;continue
   sid=d.get('shard_id');a=d.get('orbit_start');b=d.get('orbit_stop_exclusive')
   want={'mode':'scan-roworbit-j7','N':N,'left_size':l,'next_orbit_rank':b,'orbits_scanned':b-a,'max_abs_Rj7':0,'source_sha256':source}
   if not(isinstance(sid,str) and sid.startswith('rowprod-') and all(d.get(k)==v for k,v in want.items())):raise Bad(f'bad row result {p}')
   if not(isinstance(a,int) and isinstance(b,int) and 0<=a<b<=orbit_total(l) and isinstance(d.get('weighted_labeled_graphs'),int) and d['weighted_labeled_graphs']>=0):raise Bad(f'bad transposed-row interval/weight {p}')
   original=11-l;check_receipts(z.row_production,sid);rowcov[original].append((a,b));weights[original]+=d['weighted_labeled_graphs']
  fam={str(l):intervals(oldcov[l],1<<(l*(N-l))) for l in range(1,4)}
  for l in (4,5):
   sorted_transposed=11-l;fam[str(l)]=intervals(rowcov[l],orbit_total(sorted_transposed));fam[str(l)]['sorted_transposed_left']=sorted_transposed;fam[str(l)]['weighted_labeled_graphs']=weights[l];fam[str(l)]['expected_weighted_labeled_graphs']=1<<(l*(N-l));fam[str(l)]['weight_complete']=weights[l]==(1<<(l*(N-l)))
  complete=all(fam[str(l)]['complete'] for l in range(1,6)) and all(fam[str(l)].get('weight_complete',True) for l in (4,5)) and not noncomplete
  report={'schema':1,'N':N,'identity':'Rj7=h2-7h3+21h4-35h5+35h6-21h7+7h8-2h9','old_manifest_sha256':sha(z.old_manifest),'old_source_sha256':old.get('source_sha256'),'row_source_sha256':source,'families':fam,'noncompleted_row_results_excluded':noncomplete,'total_labelled_graphs':sum(1<<(l*(N-l)) for l in range(1,6)),'status':'complete_authenticated' if complete else 'partial_authenticated','strict_complete_coverage':complete}
  if report['total_labelled_graphs']!=1359217664:raise Bad('internal total mismatch')
  write_new(z.output,report);print(report['status'],z.output);return 0 if complete or not z.require_complete else 2
 except Bad as e:print('AUDIT ERROR:',e,file=sys.stderr);return 1
if __name__=='__main__':raise SystemExit(main())
