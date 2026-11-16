#!/usr/bin/env python3
"""Bounded parallel discovery of modular factor-degree obstruction certificates."""
import argparse,concurrent.futures,gzip,hashlib,json,os,sys,time
from pathlib import Path
PRIMES=(2,3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71,73,79,83,89,97)
def solve(args):
 row,deadline=args
 import sympy as sp
 x=sp.Symbol('x');f=row['coefficients_low_to_high'];n=len(f)-1;remaining=set(range(n+1));evidence=[]
 start=time.monotonic()
 for p in PRIMES:
  if time.time()>=deadline:break
  c,factors=sp.factor_list(sp.Poly.from_list(f[::-1],gens=x,modulus=p))
  if int(c)%p!=1:raise ValueError('nonmonic factorization')
  fac=[{'coefficients':[int(v.nth(i))%p for i in range(v.degree()+1)],'multiplicity':int(m)} for v,m in factors]
  subset={0}
  for v in fac:
   for _ in range(v['multiplicity']):subset|={i+len(v['coefficients'])-1 for i in list(subset)}
  new=remaining&subset
  if new!=remaining:
   evidence.append({'prime':p,'factors':fac});remaining=new
  if remaining=={0,n}:break
 return {'factor_sha256':row['factor_sha256'],'coefficients_low_to_high':f,'degree':n,'status':'candidate_degree_sieve_certificate' if remaining=={0,n} else ('unknown_time_limit' if time.time()>=deadline else 'unknown_degree_sieve'),'remaining_degrees':sorted(remaining),'evidence':evidence,'seconds':time.monotonic()-start}
def main():
 ap=argparse.ArgumentParser();ap.add_argument('--input',type=Path,required=True);ap.add_argument('--out',type=Path,required=True);ap.add_argument('--workers',type=int,default=9);ap.add_argument('--seconds',type=int,default=14400);ap.add_argument('--limit',type=int);a=ap.parse_args();assert os.path.ismount('/Volumes/data') and not a.out.exists();a.out.mkdir(parents=True)
 raw=a.input.read_bytes();doc=json.loads(gzip.decompress(raw));rows=[r for r in doc['certificates'] if r['certificate']['kind'].startswith('unknown')];rows=rows[:a.limit] if a.limit else rows
 meta={'status':'running','pid':os.getpid(),'workers':a.workers,'source_sha256':hashlib.sha256(raw).hexdigest(),'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'input_count':len(rows),'primes':PRIMES,'seconds_limit':a.seconds,'scope':'QQ irreducibility discovery for previously unknown finite-ledger factors; not real-factor enumeration'};(a.out/'input.json').write_text(json.dumps(meta,indent=2)+'\n');started=time.monotonic();deadline=time.time()+a.seconds;counts={}
 with (a.out/'rows.jsonl').open('w') as log,concurrent.futures.ProcessPoolExecutor(max_workers=a.workers) as pool:
  futures=[pool.submit(solve,(row,deadline)) for row in rows]
  for i,f in enumerate(concurrent.futures.as_completed(futures),1):
   try:r=f.result()
   except Exception as e:r={'status':'failure','error':repr(e)}
   log.write(json.dumps(r,separators=(',',':'))+'\n');log.flush();counts[r['status']]=counts.get(r['status'],0)+1
   if i%25==0 or i==len(rows):(a.out/'checkpoint.json').write_text(json.dumps({'done':i,'total':len(rows),'counts':counts,'elapsed_seconds':time.monotonic()-started})+'\n')
 meta.update(status='complete',counts=counts,elapsed_seconds=time.monotonic()-started,workers_joined=True,python=sys.version);(a.out/'result.json').write_text(json.dumps(meta,indent=2)+'\n');print(json.dumps(meta))
if __name__=='__main__':main()
