#!/usr/bin/env python3
"""Independent exact verifier for bounded P31 Q-factor allocation ledgers.

This intentionally does not import the searcher.  It independently rebuilds
the support DFS, factors every saved H over QQ, enumerates every allocation,
and validates each signed-allocation witness with integer convolutions.
"""
from __future__ import annotations

import argparse, hashlib, itertools, json, platform, sys, time
from collections import Counter
from pathlib import Path
import sympy as sp

X=sp.symbols('x'); P,Q=3,5
def conv(a,b):
 r=[0]*(len(a)+len(b)-1)
 for i,u in enumerate(a):
  for j,v in enumerate(b): r[i+j]+=u*v
 return r
def stretch(a,k):
 r=[0]*((len(a)-1)*k+1)
 for i,u in enumerate(a):r[i*k]=u
 return r
def firstneg(a):
 for i,u in enumerate(a):
  if u<0:return {'degree':i,'coefficient':str(u)}
 return None
def nonneg(a):return all(u>=0 for u in a)
def supports(maxdegree):
 def rec(s,old):
  yield s
  for t in range(s[-1]+1,maxdegree+1):
   fresh={3*t+5*u for u in s}|{3*u+5*t for u in s}|{8*t}
   if len(fresh)==2*len(s)+1 and not(fresh&old):yield from rec(s+(t,),old|frozenset(fresh))
 yield from rec((0,),frozenset({0}))
def hfrom(s):
 h=[0]*(s[-1]+1)
 for i in s:h[i]=1
 return h
def factors(h):
 poly=sp.Poly.from_list(list(reversed(h)),X,domain=sp.ZZ);content,data=sp.factor_list(poly.as_expr(),X)
 if content!=1:raise ValueError('content')
 ps=[sp.Poly(f,X,domain=sp.QQ) for f,_ in data];ms=[int(m) for _,m in data]
 if not all(p.LC()==1 and p.is_irreducible for p in ps):raise ValueError('irreducibility backend')
 cc=[]
 for p in ps:
  a=[p.nth(i) for i in range(p.degree()+1)]
  if any(u.q!=1 for u in a):raise ValueError('noninteger QQ factor')
  cc.append([int(u) for u in a])
 return cc,ms
def product(fs,ee):
 out=[1]
 for f,e in zip(fs,ee):
  for _ in range(e):out=conv(out,f)
 return out
def loadjsonl(path):
 with path.open() as f:return [json.loads(line) for line in f if line.strip()]
def main():
 ap=argparse.ArgumentParser();ap.add_argument('--search-out',type=Path,required=True);ap.add_argument('--record',type=Path,required=True);a=ap.parse_args();root=a.search_out
 meta=json.loads((root/'metadata.json').read_text());final=json.loads((root/'final.json').read_text());lo,hi=meta['scope']['degrees_inclusive'];started=time.monotonic()
 expected=[s for s in supports(hi) if lo<=s[-1]<=hi]; coverage=loadjsonl(root/'h-coverage.jsonl'); signed=loadjsonl(root/'signed-allocations.jsonl')
 if len(expected)!=len(coverage):raise SystemExit('coverage length')
 signed_by_id={r['allocation_id']:r for r in signed}
 if len(signed_by_id)!=len(signed):raise SystemExit('duplicate signed allocation id')
 totals=Counter();degreecounts=Counter();seenids=set()
 for index,(s,row) in enumerate(zip(expected,coverage)):
  h=hfrom(s);hid=f'{index:05d}-d{s[-1]}'
  if row['h_id']!=hid or row['support']!=list(s) or row['H']!=h:raise SystemExit(f'H mismatch {index}')
  if not all(u in (0,1) for u in conv(stretch(h,P),stretch(h,Q))):raise SystemExit('binary C mismatch')
  fs,ms=factors(h)
  if row['factor_coefficients_low_to_high']!=fs or row['multiplicities']!=ms:raise SystemExit(f'factor mismatch {hid}')
  local=Counter();wanted=[]
  for ee in itertools.product(*(range(m+1) for m in ms)):
   local['all']+=1;f=product(fs,ee);comp=tuple(m-e for m,e in zip(ms,ee));g=product(fs,comp)
   if conv(f,g)!=h:raise SystemExit('factor product')
   if f[0]<=0 or g[0]<=0:local['nonpositive_constant']+=1;continue
   local['positive_constant']+=1
   if nonneg(f) and nonneg(g):local['both_factor_nonnegative_skipped']+=1;continue
   local['signed_checked']+=1;aid=hid+':'+','.join(map(str,ee));wanted.append(aid);rec=signed_by_id.get(aid)
   if rec is None:raise SystemExit('missing signed record '+aid)
   if rec['F']!=f or rec['G']!=g:raise SystemExit('signed factor mismatch '+aid)
   aa=conv(stretch(f,P),stretch(g,Q));bb=conv(stretch(g,P),stretch(f,Q));na,nb=firstneg(aa),firstneg(bb)
   if rec['first_negative_A_equals_F_x3_G_x5']!=na or rec['first_negative_B_equals_G_x3_F_x5']!=nb:raise SystemExit('negative witness mismatch '+aid)
   if na is None and nb is None:
    local['both_cross_nonnegative']+=1
    if rec['status']!='candidate_both_cross_nonnegative' or rec.get('A')!=aa or rec.get('B')!=bb:raise SystemExit('candidate mismatch')
   else:
    local['cross_negative_witnessed']+=1
    if rec['status']!='rejected_by_exact_cross_negative_witness':raise SystemExit('reject status')
  if row['allocation_counts']!=dict(local) or row['signed_allocation_ids']!=wanted:raise SystemExit('coverage allocation counts '+hid)
  seenids.update(wanted);totals.update({'H_completed':1,'allocations_all':local['all'],'allocations_positive_constant':local['positive_constant'],'allocations_nonpositive_constant':local['nonpositive_constant'],'allocations_nonnegative_skipped':local['both_factor_nonnegative_skipped'],'signed_allocations_checked':local['signed_checked'],'cross_negative_witnessed':local['cross_negative_witnessed'],'both_cross_nonnegative_candidates':local['both_cross_nonnegative']});degreecounts[s[-1]]+=1
 if seenids!=set(signed_by_id):raise SystemExit('extraneous signed record')
 if dict(totals)!=final['counts']:raise SystemExit('final count mismatch')
 out={'status':'independent_exact_verification_passed','scope':meta['scope'],'actual_count_by_degree':{str(k):degreecounts[k] for k in range(lo,hi+1)},'counts':dict(totals),'signed_records_verified':len(signed),'search_final_sha256':hashlib.sha256((root/'final.json').read_bytes()).hexdigest(),'software':{'python':sys.version,'sympy':sp.__version__,'platform':platform.platform()},'elapsed_seconds':time.monotonic()-started,'claim_boundary':'Independent verification of the bounded rational allocation ledger only; no real-factor or all-degree conclusion.'}
 a.record.parent.mkdir(parents=True,exist_ok=True);a.record.write_text(json.dumps(out,indent=2,sort_keys=True)+'\n');print(json.dumps(out,sort_keys=True))
if __name__=='__main__':main()
