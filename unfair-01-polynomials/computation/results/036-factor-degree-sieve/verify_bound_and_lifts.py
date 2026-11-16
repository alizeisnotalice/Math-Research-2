#!/usr/bin/env python3
import argparse,copy,gzip,hashlib,json,time
from pathlib import Path
import integer_arithmetic as z
from verify_degree_sieve import verify

def remz(a,b):
 a=a[:]
 while len(a)>=len(b) and a!=[0]:
  c=a[-1];d=len(a)-len(b)
  for i,v in enumerate(b):a[d+i]-=c*v
  z.trim(a)
 return a

def lift(row):
 f=z.valid_poly(row['coefficients_low_to_high'],'factor');c=row['certificate']
 if z.factor_key(f)!=row['factor_sha256']:raise ValueError('lift hash')
 p=c['prime'];phi=z.valid_poly(c['phi'],'phi');e=c['exponent'];d=len(phi)-1
 if not z.prime(p) or type(e)is not int or e<2 or e*d!=len(f)-1 or any(v<0 or v>=p for v in phi):raise ValueError('lift parameters')
 if d>1:
  if z.powmod([0,1],p**d,phi,p)!=[0,1]:raise ValueError('phi Frobenius')
  for v in z.distinct_prime_divisors(d):
   if z.gcdp(z.sub(z.powmod([0,1],p**(d//v),phi,p),[0,1],p),phi,p)!=[1]:raise ValueError('phi reducible')
 product=[1]
 for _ in range(e):product=[v%p for v in z.mul(product,phi)]
 if product!=[v%p for v in f]:raise ValueError('lift mod product')
 r=remz(f,phi);padded=r+[0]*(d-len(r))
 if padded!=c['integer_remainder'] or any(v%p for v in r) or not any(v%(p*p) for v in r):raise ValueError('lift valuation')
 return True

def reject(check,row,mutate):
 v=copy.deepcopy(row);mutate(v)
 try:check(v)
 except (ValueError,SystemExit):return True
 raise AssertionError('mutation accepted')

def main():
 ap=argparse.ArgumentParser();ap.add_argument('--source',type=Path,required=True);ap.add_argument('--rows',type=Path,required=True);ap.add_argument('--lifts',type=Path,required=True);ap.add_argument('--out',type=Path,required=True);a=ap.parse_args();t=time.monotonic();raw=a.source.read_bytes();source=json.loads(gzip.decompress(raw));expected={r['factor_sha256']:r for r in source['certificates'] if r['certificate']['kind'].startswith('unknown')};rows=[json.loads(l) for l in a.rows.read_text().splitlines()];seen=set();ok=set();unknown={}
 for r in rows:
  k=r['factor_sha256']
  if k not in expected or k in seen or r['coefficients_low_to_high']!=expected[k]['coefficients_low_to_high']:raise ValueError('full source binding')
  seen.add(k)
  if verify(r):ok.add(k)
  else:unknown[k]=r
 if seen!=set(expected):raise ValueError('missing source rows')
 lifts=json.loads(a.lifts.read_text());ls=set();new=[]
 for r in lifts:
  k=r['factor_sha256']
  if k not in unknown or k in ls or r['coefficients_low_to_high']!=expected[k]['coefficients_low_to_high']:raise ValueError('lift source binding')
  ls.add(k)
  if r['certificate'] is not None:lift(r);ok.add(k);new.append(k)
 if ls!=set(unknown):raise ValueError('lift coverage')
 good=next(r for r in rows if r['factor_sha256'] in ok and r['status']=='candidate_degree_sieve_certificate');lc=next(r for r in lifts if r['certificate'])
 negatives={
 'composite_prime':reject(verify,good,lambda r:r['evidence'][0].__setitem__('prime',4)),
 'zero_multiplicity':reject(verify,good,lambda r:r['evidence'][0]['factors'][0].__setitem__('multiplicity',0)),
 'wrong_factor_product':reject(verify,good,lambda r:r['evidence'][0]['factors'][0]['coefficients'].__setitem__(0,(r['evidence'][0]['factors'][0]['coefficients'][0]+1)%r['evidence'][0]['prime'])),
 'false_degree_set':reject(verify,good,lambda r:r.__setitem__('remaining_degrees',[0,1,r['degree']])),
 'false_lift_remainder':reject(lift,lc,lambda r:r['certificate']['integer_remainder'].__setitem__(0,r['certificate']['integer_remainder'][0]+1)),
 'false_lift_exponent':reject(lift,lc,lambda r:r['certificate'].__setitem__('exponent',1)),
 }
 remaining=[expected[k] for k in sorted(set(expected)-ok)]
 out={'status':'passed','source_sha256':hashlib.sha256(raw).hexdigest(),'input_count':len(expected),'degree_sieve_certified':len(ok)-len(new),'lift_certified':len(new),'remaining_unknown':len(remaining),'remaining':remaining,'negative_tests':negatives,'elapsed_seconds':time.monotonic()-t,'verifier_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()};a.out.write_text(json.dumps(out,indent=2)+'\n');print(json.dumps({k:v for k,v in out.items() if k!='remaining'}))
if __name__=='__main__':main()
