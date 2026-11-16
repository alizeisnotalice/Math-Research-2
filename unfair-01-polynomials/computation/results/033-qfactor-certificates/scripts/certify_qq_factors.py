#!/usr/bin/env python3
"""Produce checkable irreducibility evidence for QQ factors in a ledger.

Cyclotomic factors are certified by equality with an exactly generated Phi_n.
Other factors are certified only when a displayed finite-field Rabin test at
one prime proves the reduction irreducible; remaining factors are explicitly
reported as unknown, never silently treated as certified.
"""
from __future__ import annotations
import argparse, hashlib, json, platform, sys, time
from pathlib import Path
import sympy as sp

X=sp.symbols('x'); PRIMES=(2,3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71,73,79,83,89,97)
def trim(a):
 while len(a)>1 and a[-1]==0:a.pop()
 return a
def add(a,b,p):
 r=[0]*max(len(a),len(b))
 for i in range(len(a)):r[i]+=a[i]
 for i in range(len(b)):r[i]+=b[i]
 return trim([u%p for u in r])
def sub(a,b,p):
 r=[0]*max(len(a),len(b))
 for i in range(len(a)):r[i]+=a[i]
 for i in range(len(b)):r[i]-=b[i]
 return trim([u%p for u in r])
def divmodp(a,b,p):
 a=trim(a[:]);b=trim(b[:]);q=[0]*max(1,len(a)-len(b)+1);ib=pow(b[-1],-1,p)
 while len(a)>=len(b) and a!=[0]:
  c=a[-1]*ib%p;d=len(a)-len(b);q[d]=c
  for i,u in enumerate(b):a[d+i]=(a[d+i]-c*u)%p
  trim(a)
 return trim(q),a
def gcdp(a,b,p):
 while b!=[0]:_,r=divmodp(a,b,p);a,b=b,r
 inv=pow(a[-1],-1,p);return trim([u*inv%p for u in a])
def mulmod(a,b,f,p):
 n=len(f)-1;r=[0]*(len(a)+len(b)-1)
 for i,u in enumerate(a):
  for j,v in enumerate(b):r[i+j]=(r[i+j]+u*v)%p
 for k in range(len(r)-1,n-1,-1):
  c=r[k]%p
  if c:
   for j in range(n):r[k-n+j]=(r[k-n+j]-c*f[j])%p
  r.pop()
 return trim(r or [0])
def powmod(a,e,f,p):
 out=[1]
 while e:
  if e&1:out=mulmod(out,a,f,p)
  a=mulmod(a,a,f,p);e//=2
 return out
def prime_divisors(n):
 out=[];d=2
 while d*d<=n:
  if n%d==0:
   out.append(d)
   while n%d==0:n//=d
  d+=1
 if n>1:out.append(n)
 return out
def rabin_certificate(coeffs,p):
 n=len(coeffs)-1;f=[u%p for u in coeffs];xp=powmod([0,1],p**n,f,p);checks=[]
 for r in prime_divisors(n):
  z=sub(powmod([0,1],p**(n//r),f,p),[0,1],p);g=gcdp(z,f,p);checks.append({'prime_divisor_of_degree':r,'gcd_coefficients_low_to_high_mod_p':g})
 if xp==[0,1] and all(v['gcd_coefficients_low_to_high_mod_p']==[1] for v in checks):
  return {'kind':'irreducible_mod_prime_rabin','prime':p,'degree':n,'x_to_p_to_degree_mod_factor_coefficients_low_to_high':xp,'rabin_gcd_checks':checks}
 return None
def factor_key(a):return hashlib.sha256(json.dumps(a,separators=(',',':')).encode()).hexdigest()
def cyclotomic_table(max_n=500):
 table={}
 for n in range(1,max_n+1):
  p=sp.Poly(sp.cyclotomic_poly(n,X),X,domain=sp.ZZ);a=[int(p.nth(i)) for i in range(p.degree()+1)];table.setdefault(tuple(a),n)
 return table
def main():
 ap=argparse.ArgumentParser();ap.add_argument('--search-out',type=Path,required=True);ap.add_argument('--record',type=Path,required=True);a=ap.parse_args();started=time.monotonic();unique={}
 with (a.search_out/'h-coverage.jsonl').open() as s:
  for line in s:
   row=json.loads(line)
   for f in row['factor_coefficients_low_to_high']:unique.setdefault(tuple(f),[]).append(row['h_id'])
 cy=cyclotomic_table();certs=[];counts={'unique_factors':len(unique),'cyclotomic_exact':0,'degree_one_exact':0,'rabin_mod_prime':0,'unknown':0}
 for coeff_tuple,hids in sorted(unique.items()):
  coeffs=list(coeff_tuple);key=factor_key(coeffs);degree=len(coeffs)-1;entry={'factor_sha256':key,'coefficients_low_to_high':coeffs,'degree':degree,'used_by_H_count':len(hids),'first_H_ids':hids[:10]}
  if degree==1:
   entry['certificate']={'kind':'degree_one_exact'};counts['degree_one_exact']+=1
  elif coeff_tuple in cy:
   entry['certificate']={'kind':'cyclotomic_exact','n':cy[coeff_tuple]};counts['cyclotomic_exact']+=1
  else:
   cert=None
   for p in PRIMES:
    cert=rabin_certificate(coeffs,p)
    if cert is not None:break
   if cert is None:
    entry['certificate']={'kind':'unknown_no_mod_prime_certificate_in_fixed_prime_list','primes_tried':list(PRIMES)};counts['unknown']+=1
   else:
    entry['certificate']=cert;counts['rabin_mod_prime']+=1
  certs.append(entry)
 out={'status':'factor_certificate_scan_complete','scope':'Certificates apply only to distinct factors appearing in the supplied bounded ledger. Unknown means factorization completeness remains unproved by this file.','counts':counts,'certificates':certs,'software':{'python':sys.version,'sympy':sp.__version__,'platform':platform.platform()},'elapsed_seconds':time.monotonic()-started}
 a.record.parent.mkdir(parents=True,exist_ok=True);a.record.write_text(json.dumps(out,indent=2,sort_keys=True)+'\n');print(json.dumps({'counts':counts,'sha256':hashlib.sha256(a.record.read_bytes()).hexdigest()},sort_keys=True))
if __name__=='__main__':main()
