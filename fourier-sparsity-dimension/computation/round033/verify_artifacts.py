"""Check final bindings and independently re-evaluate the explicit seed.

This does not prove the all-r extension: see the analytic proof.
"""
from pathlib import Path
import json,hashlib
P=Path(__file__).resolve().parent
read=lambda n:json.loads((P/n).read_text())
sha=lambda n:hashlib.sha256((P/n).read_bytes()).hexdigest()
m=read('mm_completion_result.json');a=read('independent-mm-witness-audit.json');t=read('tensor_shear_result.json')
for result,src in [(m,'search_coordinate_mm.py'),(a,'independent_mm_witness_audit.py'),(t,'check_tensor_shear.py')]:assert result['source_sha256']==sha(src)
assert a['input_result_sha256']==t['input_sha256']==sha('mm_completion_result.json')
s=m['solution'];assert sorted(s['pi'])==list(range(16))
f=int(s['f_hex'],16);g=int(s['f_plus_h_hex'],16)
def fwht(z):
 h=1
 while h<len(z):
  for b in range(0,len(z),2*h):
   for i in range(b,b+h):z[i],z[i+h]=z[i]+z[i+h],z[i]-z[i+h]
  h*=2
 return z
F=fwht([1-2*(f>>i&1) for i in range(256)]);G=fwht([1-2*(g>>i&1) for i in range(256)])
assert set(F)==set(G)=={-16,16}
S={i for i in range(256) if (f^g)>>i&1};assert len(S)==32 and {x^y for x in S for y in S}==set(range(128))
for x in range(256):
 b=[x>>i&1 for i in range(8)];h=(1^b[7])*((b[0]*b[1]*b[3])^(b[0]*b[4]*b[5])^(b[1]*b[2]*b[6]))
 assert h==((f^g)>>x&1)
D={i for i in range(256) if F[i]!=G[i]};DD={x^y for x in D for y in D};K={0}
for v in [97,98,20]:K|={x^v for x in list(K)}
assert len(D)==32 and len(DD)==182 and len(K)==8 and K&DD=={0}
assert len(t['records'])==4
for r,row in enumerate(t['records']):
 assert row['r']==r and row['tensor_both_bent'] and row['shear_both_bent'] and row['shear_difference_formula_verified']
 assert row['tensor_switch_weight']==32*(1<<r) and row['shear_switch_weight']==64*(1<<r)
remote=read('macbook_tensor_metal_result.json');receipt=read('macbook_run_receipt.json')
assert receipt['source_sha256']==sha('R33MMTensorMetal.swift') and receipt['output_sha256']==sha('macbook_tensor_metal_result.json')
assert receipt['compile_exit_code']==receipt['run_exit_code']==0
assert receipt['canonical_solution_sha256']==hashlib.sha256(json.dumps(s,sort_keys=True).encode()).hexdigest()
assert remote['result_input_sha256']==receipt['result_input_sha256'] and remote['signed_mask_input_sha256']==receipt['signed_input_sha256']
assert remote['all_checks_passed'] and set(remote['difference_support'])==S and len(remote['rows'])==5
for r,row in enumerate(remote['rows']):
 assert row['r']==r and row['dimension']==8+2*r and row['difference_weight']==32*(1<<r)
 assert row['gpu_cpu_mismatch_count']==row['cpu_bad_absq']==row['gpu_bad_absq']==0
print('PASS: explicit bent pair and nonprojection seed; dual kernel; local tensor/shear checks; MacBook CPU/Metal bindings and equality')
