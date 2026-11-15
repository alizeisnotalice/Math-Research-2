"""Independent NumPy reconstruction of retained numerical root masks.

No searcher import; never converts a floating negative result to exclusion.
"""
import json,hashlib,sys
from pathlib import Path
import numpy as np

def check(input_path,report_path):
 inp=Path(input_path);report=json.loads(Path(report_path).read_text());inputs=json.loads(inp.read_text())['candidates'];out=[];total=0
 source=Path(__file__).resolve().parent/report['source']
 for run in report['runs']:
  assert run['source_sha256']==hashlib.sha256(source.read_bytes()).hexdigest()
  if 'v1_helper_sha256' in run:
   assert run['v1_helper_sha256']==hashlib.sha256((source.parent/'gpu_root_group_search.py').read_bytes()).hexdigest()
  assert run['input_sha256']==hashlib.sha256(inp.read_bytes()).hexdigest()
  c=inputs[run['prepared_candidate']['index']];C=c['coefficients_ascending'];groups=c['groups_ascending_coefficients'];g=len(groups)
  assert c['sha256']==hashlib.sha256(json.dumps(C,separators=(',',':')).encode()).hexdigest()
  assert run['status']=='complete' and run['completed_mask_half_open']==[0,1<<(g-1)]
  assert run['counts']['codes_seen']==1<<(g-1);total+=1<<(g-1)
  tol=run['parameters']['tol'];eta=run['parameters']['eta'];near=[];rounded_fair=[];max_product_error=0.
  for item in run['candidates']:
   mask=item['mask'];assert 0<mask<(1<<(g-1));A=np.array([1.]);B=np.array([1.])
   for j,f in enumerate(groups):
    if (mask>>j)&1:A=np.convolve(A,np.array(f))
    else:B=np.convolve(B,np.array(f))
   A/=A[0];B/=B[0];old=item['cpu_float64_recheck']
   assert np.allclose(A,old['A_coefficients_ascending'],rtol=1e-13,atol=1e-12)
   assert np.allclose(B,old['B_coefficients_ascending'],rtol=1e-13,atol=1e-12)
   z=np.r_[A,B];err=float(np.max(abs(np.convolve(A,B)-C)));max_product_error=max(max_product_error,err)
   assert abs(err-old['convolution_max_abs_residual'])<1e-9
   if np.all(z>=-tol) and np.all(z<=1+tol) and np.any((z>eta)&(z<1-eta)) and abs(A[-1]-1)<=tol and abs(B[-1]-1)<=tol:near.append(mask)
   if np.max(np.minimum(abs(z),abs(z-1)))<1e-6:
    aa=[int(round(x)) for x in A];bb=[int(round(x)) for x in B]
    cc=[sum(aa[i]*bb[k-i] for i in range(len(aa)) if 0<=k-i<len(bb)) for k in range(len(C))]
    if cc==C:rounded_fair.append(mask)
  out.append({'name':c['name'],'recorded_numeric_masks':1<<(g-1),'retained_recomputed':len(run['candidates']),'numerical_fractional_survivors':near,'rounded_samples_with_exact_fair_convolution':rounded_fair,'maximum_retained_float64_product_error':max_product_error})
 return {'recorded_numeric_masks':total,'results':out,'scope':'Only retained samples independently recomputed. Full numerical mask scans are not exact exclusions.'}
if __name__=='__main__':print(json.dumps(check(sys.argv[1],sys.argv[2]),indent=2))
