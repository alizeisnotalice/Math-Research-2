#!/usr/bin/env python3
import json,hashlib
from pathlib import Path
import sympy as s
x=s.Symbol('x');out=[]
for line in Path('full/rows.jsonl').read_text().splitlines():
 r=json.loads(line)
 if r['status']=='candidate_degree_sieve_certificate':continue
 f=s.Poly.from_list(r['coefficients_low_to_high'][::-1],x,domain=s.ZZ);cert=None
 for p in (2,3,5,7,11,13,17,19,23,29,31):
  _,ff=s.factor_list(f,modulus=p)
  if len(ff)!=1:continue
  v,e=ff[0];a=[int(v.nth(i))%p for i in range(v.degree()+1)];phi=s.Poly.from_list(a[::-1],x,domain=s.ZZ);rem=f.rem(phi);b=[int(rem.nth(i)) for i in range(phi.degree())]
  if all(c%p==0 for c in b) and any(c%(p*p)!=0 for c in b):cert={'prime':p,'phi':a,'exponent':int(e),'integer_remainder':b};break
 out.append({'factor_sha256':r['factor_sha256'],'coefficients_low_to_high':r['coefficients_low_to_high'],'certificate':cert})
Path('lift-eisenstein.json').write_text(json.dumps(out,indent=2)+'\n');print([(len(r['coefficients_low_to_high'])-1, r['certificate']) for r in out])
