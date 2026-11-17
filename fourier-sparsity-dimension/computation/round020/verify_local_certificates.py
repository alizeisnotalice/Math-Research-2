"""Independent direct-collision verification of archived R20 local certificates."""
from pathlib import Path
import hashlib,json
P=Path(__file__).resolve().parent
def read(n):return json.loads((P/n).read_text())
def sha(n):return hashlib.sha256((P/n).read_bytes()).hexdigest()
def main():
 inputs=read('mm-cases.json'); audit=read('mm-cases-constructive-kernel-audit.json')
 assert inputs['source_sha256']==sha('mm_cases.py')
 assert audit['source_sha256']==sha('mm_cases_constructive_kernel_audit.py')
 assert audit['input_sha256']==sha('mm-cases.json')
 assert inputs['count']==audit['cases']==293
 counts={}
 for a,r in zip(inputs['cases'],audit['rows'],strict=True):
  A=int(a['support_hex'],16);assert A==int(r['support_hex'],16) and A.bit_count()==32
  K={0}
  for v in r['kernel_basis']:K|={z^v for z in tuple(K)}
  assert len(K)==8
  for x in range(256):
   if A>>x&1:assert all(not(A>>(x^z)&1) for z in K-{0})
  kind=tuple(r['case_type']);counts[kind]=counts.get(kind,0)+1
 assert counts=={(4,0):288,(2,1):4,(0,2):1}
 ex=read('signed_autocorrelation_obstacle_k2.json')
 assert ex['source_sha256']==sha('verify_signed_autocorrelation_obstacle_k2.py')
 assert ex['base_and_flip_bent_verified'] and ex['signed_autocorrelation']==0 and ex['unsigned_pair_count']==8
 print(json.dumps({'status':'passed','explicit_kernels_checked':293,'method':'direct x and x+kernel point collisions; not difference masks','source_input_bindings':'passed'}))
if __name__=='__main__':main()
