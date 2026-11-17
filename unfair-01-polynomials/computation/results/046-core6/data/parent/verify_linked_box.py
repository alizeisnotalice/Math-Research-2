#!/usr/bin/env python3
"""Independent fixed-box enumeration, avoiding signed ceil/floor range formulas."""
import math,json,hashlib,time
from pathlib import Path
HERE=Path(__file__).resolve().parent
SOURCE=Path('/Volumes/data/CodexTemp/T2/unfair-01-polynomials/agents/coprime_stretches/round046/core6_linked/core6_linked_conditional_ledger_corrected.json')
def main():
 start=time.monotonic();rows=[]
 for b in range(-7,8):
  if not b:continue
  for t in range(b*b-12,b*b+13):
   for n in range(-473,474):
    delta=t*t-4*n
    if delta<=0 or math.isqrt(delta)**2==delta:continue
    s=51-4*abs(t-b*b)
    if s<=0 or 16*delta>=s*s:continue
    s=87-2*abs(t)
    if s<=0 or 4*delta>=s*s:continue
    for ell in range(-123,124):
     if 8*abs(b**3+3*ell)>=195:continue
     s=245-4*abs(b*t+2*ell)
     if s<=0 or 16*b*b*delta>=s*s:continue
     rows.append((b,t,n,ell))
 seen=set(rows);assert len(rows)==len(seen)
 assert {(-b,t,n,-ell) for b,t,n,ell in rows}==seen
 observed=json.loads(SOURCE.read_text())
 supplied=[(r['b'],r['T'],r['N'],r['ell']) for r in observed['rows']]
 assert len(supplied)==len(set(supplied)) and set(supplied)==seen
 assert len(rows)==35488
 out={'status':'passed_independent_fixed_box_reconstruction','rows':len(rows),'d_e_expansion':len(rows)*43*15,'source_sha256':hashlib.sha256(SOURCE.read_bytes()).hexdigest(),'verifier_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'seconds':time.monotonic()-start,'scope':'Complete parameter input only; no prefix or full-model exclusion yet.'}
 (HERE/'linked-box-check.json').write_text(json.dumps(out,indent=2)+'\n')
 (HERE/'linked-input-tuples.json').write_text(json.dumps(rows,separators=(',',':'))+'\n');print(json.dumps(out))
if __name__=='__main__':main()
