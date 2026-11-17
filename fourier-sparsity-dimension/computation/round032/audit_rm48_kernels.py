#!/usr/bin/env python3
"""Nine-process exact RREF kernel avoidance audit."""
from __future__ import annotations
import concurrent.futures as cf, hashlib, json, os, platform, sys, time
from pathlib import Path
ROOT=Path(__file__).resolve().parent
SRC=ROOT/'rm48_weight32_representatives.json'
OUT=ROOT/'rm48_weight32_kernel_audit.json'
def rank(vs):
 b={}
 for v in vs:
  while v:
   p=v.bit_length()-1
   if p in b:v^=b[p]
   else:b[p]=v;break
 return len(b)
def rref3():
 count=0
 for p0 in range(8):
  for p1 in range(p0+1,8):
   for p2 in range(p1+1,8):
    piv={p0,p1,p2};slots=[(i,q) for i,p in enumerate((p0,p1,p2)) for q in range(8) if q not in piv and q>p]
    for a in range(1<<len(slots)):
     r=[1<<p0,1<<p1,1<<p2]
     for j,(i,q) in enumerate(slots):
      if a>>j&1:r[i]|=1<<q
     K=(1<<0)|(1<<r[0])|(1<<r[1])|(1<<(r[0]^r[1]))|(1<<r[2])|(1<<(r[0]^r[2]))|(1<<(r[1]^r[2]))|(1<<(r[0]^r[1]^r[2]))
     count+=1;yield tuple(r),K
 assert count==97155
def value(rep,x):
 if x>>rep['affine_span_dimension']:return 0
 z=rep['constant_term']
 for mon in rep['monomials_zero_based']:
  t=1
  for i in mon:t&=(x>>i)&1
  z^=t
 return z
def diff(pts):
 out=0
 for x in pts:
  for y in pts:out|=1<<(x^y)
 return out
def degree(pts):
 a=[int(x in pts) for x in range(256)]
 for i in range(8):
  for x in range(256):
   if x>>i&1:a[x]^=a[x^(1<<i)]
 return max((x.bit_count() for x,v in enumerate(a) if v),default=0)
def work(rep):
 t=time.perf_counter();pts=[x for x in range(256) if value(rep,x)];D=diff(pts);n=0;first=None;alln=0
 for basis,K in rref3():
  alln+=1
  if ((K&~1)&D)==0:
   n+=1
   if first is None:first=list(basis)
 return {'id':rep['id'],'pid':os.getpid(),'kernel_subspaces_checked':alln,'avoiding_kernel_count':n,'first_rref_basis':first,'weight':len(pts),'affine_span_dimension_computed':rank([x^pts[0] for x in pts]),'anf_degree_computed':degree(set(pts)),'difference_size':D.bit_count(),'elapsed_seconds':time.perf_counter()-t}
def main():
 raw=SRC.read_bytes();reps=json.loads(raw)['representatives'];print(json.dumps({'parent_pid':os.getpid(),'workers':9,'status':'started'}),flush=True);t=time.perf_counter()
 with cf.ProcessPoolExecutor(max_workers=9) as ex: rows=list(ex.map(work,reps))
 out={'scope':'Exact nine-representative RREF audit only.','input_sha256':hashlib.sha256(raw).hexdigest(),'python':sys.version,'platform':platform.platform(),'parent_pid':os.getpid(),'max_workers':9,'representatives':rows,'elapsed_seconds':time.perf_counter()-t,'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
 OUT.write_text(json.dumps(out,ensure_ascii=False,indent=2)+'\n',encoding='utf-8');print(json.dumps(out,ensure_ascii=False,indent=2),flush=True)
 assert all(x['kernel_subspaces_checked']==97155 and x['weight']==32 for x in rows)
if __name__=='__main__':main()
