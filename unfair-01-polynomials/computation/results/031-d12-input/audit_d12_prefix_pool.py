#!/usr/bin/env python3
"""CPU-only, standard-library independent d=12 prefix-pool audit."""
import argparse, hashlib, itertools, json
from pathlib import Path

PAIR=((0,0),(1,0),(0,1))
def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()
def decode(code):
    a,b=[1],[1]
    for c in code:
        x,y=PAIR[int(c)]; a.append(x); b.append(y)
    return a,b
def product(a,b,end):
    return [sum(a[i]*b[k-i] for i in range(max(0,k-len(b)+1),min(k,len(a)-1)+1)) for k in range(end+1)]
def legal(code,d):
    a,b=decode(code)
    return set(product(a,b,d-1))<={0,1} and sum(a[i]*b[d-i] for i in range(1,d))==0
def filter_reasons(code,d):
    a,b=decode(code); p=product(a,b,2*d-2); out=[]
    if any(x>1 for x in p): out.append('prefix_product_coefficient_gt_1')
    for i in range(1,d):
        baseline=p[d+i] if d+i<len(p) else 0
        if baseline>=1 and a[i]+b[i]>0:
            out.append('forced_positive_fractional_cross_at_d_plus_i'); break
    return out
def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--worker',type=Path,required=True); ap.add_argument('--driver',type=Path,required=True); ap.add_argument('--module',type=Path,required=True); ap.add_argument('--record',type=Path,required=True); x=ap.parse_args(); d=12
    raw=[''.join(c) for c in itertools.product('012',repeat=d-1)]
    valid=[c for c in raw if legal(c,d)]; survivors=[c for c in valid if not filter_reasons(c,d)]
    digest=lambda z:hashlib.sha256(json.dumps(z,separators=(',',':')).encode()).hexdigest()
    worker=x.worker.read_text(); driver=x.driver.read_text(); module=x.module.read_text()
    out={'status':'cpu_only_d12_prefix_pool_audited','d':d,'model':'i<d pairs in {(0,0),(1,0),(0,1)}; low convolution binary; d-cross zero; named filters are necessary only','raw_ternary_count':len(raw),'valid_low_prefix_count':len(valid),'survivor_count':len(survivors),'valid_codes_sha256':digest(valid),'survivor_codes_sha256':digest(survivors),'all_survivors_checked_legal_and_filter_surviving':True,'filter_rejections':{'prefix_product_coefficient_gt_1':sum('prefix_product_coefficient_gt_1' in filter_reasons(c,d) for c in valid),'forced_positive_fractional_cross_at_d_plus_i':sum('forced_positive_fractional_cross_at_d_plus_i' in filter_reasons(c,d) for c in valid)},'source_sha256':{'worker':sha(x.worker),'driver':sha(x.driver),'module':sha(x.module)},'frozen_implementation_gate':{'prefix_module_general_enumerator':'def enumerate_valid_codes(d: int)' in module,'driver_accepts_d_argument':"('d',7)" in driver,'worker_explicitly_rejects_d12':'not 1 <= d <= 11' in worker,'worker_gate_conclusion':'The explicit rejection flag records the inspected source gate; GPU execution additionally requires a matching independent CPU-only audit.'},'proposed_parameters':{'m':52,'n':92,'d':12,'starts':32768,'steps':12000,'restart_every':2000,'margin':0.001,'seed':2026111610,'mathematical_shape_gate':'m,n>d and starts>=survivor_count','shape_gate_passes':52>d and 92>d and 32768>=len(survivors),'scope_note':'margin=.001 defines only the planned numerical submodel; it is not an exclusion theorem.'},'claim_boundary':'Exact finite prefix-pool enumeration and source-gate audit only; no GPU was launched, and neither survivor status nor a numerical search settles factorization.'}
    x.record.parent.mkdir(parents=True,exist_ok=True); x.record.write_text(json.dumps(out,indent=2,sort_keys=True)+'\n'); print(json.dumps({k:out[k] for k in ('raw_ternary_count','valid_low_prefix_count','survivor_count','valid_codes_sha256','survivor_codes_sha256')},sort_keys=True))
if __name__=='__main__': main()
