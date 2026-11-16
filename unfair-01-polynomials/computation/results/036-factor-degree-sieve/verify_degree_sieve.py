#!/usr/bin/env python3
"""Integer-only verification, independent of the SymPy discovery factorizer."""
import argparse,hashlib,json,time
from pathlib import Path
import integer_arithmetic as z

def verify(row):
 f=z.valid_poly(row['coefficients_low_to_high'],'factor');n=len(f)-1
 if row['degree']!=n or z.factor_key(f)!=row['factor_sha256']:raise ValueError('factor binding')
 possible=set(range(n+1));primes=set()
 for e in row['evidence']:
  p=e['prime']
  if not z.prime(p) or p in primes:raise ValueError('prime')
  primes.add(p);product=[1];degrees=[]
  for item in e['factors']:
   g=z.valid_poly(item['coefficients'],'modular factor');m=item['multiplicity'];d=len(g)-1
   if type(m)is not int or m<1 or any(v<0 or v>=p for v in g):raise ValueError('modular encoding')
   if d>1:
    if z.powmod([0,1],p**d,g,p)!=[0,1]:raise ValueError('Rabin power')
    for r in z.distinct_prime_divisors(d):
     if z.gcdp(z.sub(z.powmod([0,1],p**(d//r),g,p),[0,1],p),g,p)!=[1]:raise ValueError('Rabin gcd')
   for _ in range(m):product=[c%p for c in z.mul(product,g)];degrees.append(d)
  if product!=[c%p for c in f]:raise ValueError('factor product')
  subset={0}
  for d in degrees:subset.update({k+d for k in tuple(subset)})
  possible.intersection_update(subset)
 if sorted(possible)!=row['remaining_degrees']:raise ValueError('degree sieve mismatch')
 if row['status']=='candidate_degree_sieve_certificate' and possible!={0,n}:raise ValueError('false certificate')
 return possible=={0,n}
def main():
 ap=argparse.ArgumentParser();ap.add_argument('--rows',type=Path,required=True);ap.add_argument('--record',type=Path,required=True);a=ap.parse_args();t=time.monotonic();counts={'certified':0,'unknown':0,'failed_rows':0};seen=set()
 for line in a.rows.read_text().splitlines():
  row=json.loads(line)
  if row['status']=='failure':counts['failed_rows']+=1;continue
  if row['factor_sha256'] in seen:raise ValueError('duplicate')
  seen.add(row['factor_sha256']);ok=verify(row);counts['certified' if ok else 'unknown']+=1
 record={'status':'passed','counts':counts,'seconds':time.monotonic()-t,'scope':'Only supplied distinct factor rows; full coverage requires separate input binding','verifier_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()};a.record.write_text(json.dumps(record,indent=2)+'\n');print(json.dumps(record))
if __name__=='__main__':main()
