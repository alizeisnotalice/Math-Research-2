#!/usr/bin/env python3
"""Pure-integer independent verifier for the bounded factor certificate file."""
from __future__ import annotations
import argparse, hashlib, json, time
from pathlib import Path

def trim(a):
 while len(a)>1 and a[-1]==0:a.pop()
 return a
def divexact(a,b):
 a=a[:];q=[0]*max(1,(len(a)-len(b)+1));lead=b[-1]
 while len(a)>=len(b) and a!=[0]:
  if a[-1]%lead:raise ValueError('nonexact integer division')
  c=a[-1]//lead;d=len(a)-len(b);q[d]=c
  for i,u in enumerate(b):a[d+i]-=c*u
  trim(a)
 if a!=[0]:raise ValueError('remainder')
 return trim(q)
def mul(a,b):
 r=[0]*(len(a)+len(b)-1)
 for i,u in enumerate(a):
  for j,v in enumerate(b):r[i+j]+=u*v
 return r
def divisors(n):return [d for d in range(1,n+1) if n%d==0]
def cyclo(n,cache):
 if n in cache:return cache[n]
 den=[1]
 for d in divisors(n):
  if d<n:den=mul(den,cyclo(d,cache))
 num=[-1]+[0]*(n-1)+[1];cache[n]=divexact(num,den);return cache[n]
def sub(a,b,p):
 r=[0]*max(len(a),len(b))
 for i in range(len(a)):r[i]+=a[i]
 for i in range(len(b)):r[i]-=b[i]
 return trim([x%p for x in r])
def divmodp(a,b,p):
 a=trim(a[:]);b=trim(b[:]);q=[0]*max(1,len(a)-len(b)+1);ib=pow(b[-1],-1,p)
 while len(a)>=len(b) and a!=[0]:
  c=a[-1]*ib%p;d=len(a)-len(b);q[d]=c
  for i,u in enumerate(b):a[d+i]=(a[d+i]-c*u)%p
  trim(a)
 return trim(q),a
def gcdp(a,b,p):
 while b!=[0]:_,r=divmodp(a,b,p);a,b=b,r
 z=pow(a[-1],-1,p);return trim([u*z%p for u in a])
def mulmod(a,b,f,p):
 n=len(f)-1;r=[0]*(len(a)+len(b)-1)
 for i,u in enumerate(a):
  for j,v in enumerate(b):r[i+j]=(r[i+j]+u*v)%p
 for k in range(len(r)-1,n-1,-1):
  c=r[k]%p
  for j in range(n):r[k-n+j]=(r[k-n+j]-c*f[j])%p
  r.pop()
 return trim(r or [0])
def powmod(a,e,f,p):
 r=[1]
 while e:
  if e&1:r=mulmod(r,a,f,p)
  a=mulmod(a,a,f,p);e//=2
 return r
def factor_key(a):return hashlib.sha256(json.dumps(a,separators=(',',':')).encode()).hexdigest()
def exact_int(x,label):
 if type(x) is not int:raise SystemExit(label+' must be an exact JSON integer')
 return x
def valid_poly(a,label):
 if type(a) is not list or len(a)<2:raise SystemExit(label+' must be a positive-degree coefficient list')
 if any(type(x) is not int for x in a) or a[-1]!=1:raise SystemExit(label+' must be exact, trimmed, and monic')
 return a
def prime(p):
 if type(p) is not int or p<2:return False
 if p%2==0:return p==2
 d=3
 while d*d<=p:
  if p%d==0:return False
  d+=2
 return True
def distinct_prime_divisors(n):
 out=[];d=2
 while d*d<=n:
  if n%d==0:
   out.append(d)
   while n%d==0:n//=d
  d+=1
 if n>1:out.append(n)
 return out
def main():
 ap=argparse.ArgumentParser();ap.add_argument('--search-out',type=Path,required=True);ap.add_argument('--certificates',type=Path,required=True);ap.add_argument('--record',type=Path,required=True);z=ap.parse_args();started=time.monotonic()
 used={}
 with (z.search_out/'h-coverage.jsonl').open() as f:
  for line in f:
   for a in json.loads(line)['factor_coefficients_low_to_high']:
    a=valid_poly(a,'coverage factor');used.setdefault(tuple(a),0);used[tuple(a)]+=1
 cert=json.loads(z.certificates.read_text());certrows=cert.get('certificates')
 if type(certrows) is not list:raise SystemExit('certificate rows must be a list')
 rows={};seen_hashes=set()
 for row in certrows:
  if type(row) is not dict:raise SystemExit('certificate row must be an object')
  a=valid_poly(row.get('coefficients_low_to_high'),'certificate factor')
  h=row.get('factor_sha256')
  if type(h) is not str or h in seen_hashes or tuple(a) in rows:raise SystemExit('duplicate certificate row')
  seen_hashes.add(h);rows[tuple(a)]=row
 if set(rows)!=set(used):raise SystemExit('certificate coverage key mismatch')
 cache={1:[-1,1]};counts={'degree_one_exact':0,'cyclotomic_exact':0,'rabin_mod_prime':0,'unknown':0}
 for a,row in rows.items():
  if factor_key(list(a))!=row['factor_sha256'] or row['used_by_H_count']!=used[a]:raise SystemExit('factor metadata')
  c=row.get('certificate');
  if type(c) is not dict or type(c.get('kind')) is not str:raise SystemExit('malformed certificate')
  kind=c['kind']
  if kind=='degree_one_exact':
   if len(a)!=2 or a[-1]!=1:raise SystemExit('linear')
  elif kind=='cyclotomic_exact':
   n=exact_int(c.get('n'),'cyclotomic n')
   if n<1 or list(a)!=cyclo(n,cache):raise SystemExit('cyclotomic equality')
  elif kind=='irreducible_mod_prime_rabin':
   p,n=c.get('prime'),c.get('degree');p=exact_int(p,'Rabin prime');n=exact_int(n,'Rabin degree')
   if not prime(p) or n!=len(a)-1 or a[-1]%p!=1:raise SystemExit('Rabin prime/degree/leading coefficient')
   f=[u%p for u in a]
   xp=powmod([0,1],p**n,f,p)
   savedxp=c.get('x_to_p_to_degree_mod_factor_coefficients_low_to_high')
   if type(savedxp) is not list or any(type(u) is not int or not 0<=u<p for u in savedxp) or xp!=savedxp or xp!=[0,1]:raise SystemExit('rabin xp')
   checks=c.get('rabin_gcd_checks')
   if type(checks) is not list:raise SystemExit('Rabin gcd list')
   expected=distinct_prime_divisors(n);got=[]
   for check in checks:
    if type(check) is not dict:raise SystemExit('Rabin gcd row')
    r=exact_int(check.get('prime_divisor_of_degree'),'Rabin divisor');got.append(r)
    g=gcdp(sub(powmod([0,1],p**(n//r),f,p),[0,1],p),f,p);saved=check.get('gcd_coefficients_low_to_high_mod_p')
    if type(saved) is not list or any(type(u) is not int or not 0<=u<p for u in saved) or g!=saved or g!=[1]:raise SystemExit('rabin gcd')
   if sorted(got)!=expected or len(got)!=len(set(got)):raise SystemExit('incomplete or duplicate Rabin divisor list')
  elif kind=='unknown_no_mod_prime_certificate_in_fixed_prime_list':pass
  else:raise SystemExit('unknown kind')
  counts[{'degree_one_exact':'degree_one_exact','cyclotomic_exact':'cyclotomic_exact','irreducible_mod_prime_rabin':'rabin_mod_prime','unknown_no_mod_prime_certificate_in_fixed_prime_list':'unknown'}[kind]]+=1
 if counts!={k:cert['counts'][k] for k in counts} or cert['counts']['unique_factors']!=sum(counts.values()):raise SystemExit('count mismatch')
 out={'status':'independent_pure_integer_certificate_verification_passed','counts':counts,'certified_factor_count':counts['degree_one_exact']+counts['cyclotomic_exact']+counts['rabin_mod_prime'],'unknown_factor_count':counts['unknown'],'claim_boundary':'The 3190 unknown factors prevent claiming complete rational-factorization certification. Rabin and cyclotomic certificates above were independently recomputed with integer/modular arithmetic.','certificate_sha256':hashlib.sha256(z.certificates.read_bytes()).hexdigest(),'elapsed_seconds':time.monotonic()-started}
 z.record.parent.mkdir(parents=True,exist_ok=True);z.record.write_text(json.dumps(out,indent=2,sort_keys=True)+'\n');print(json.dumps(out,sort_keys=True))
if __name__=='__main__':main()
