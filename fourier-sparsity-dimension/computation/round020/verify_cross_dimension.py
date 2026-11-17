"""Recover MM permutations/offsets and directly verify all saved projection kernels."""
from pathlib import Path
import hashlib,json,platform,time
P=Path(__file__).resolve().parent
def sha(n):return hashlib.sha256((P/n).read_bytes()).hexdigest()
def read(n):return json.loads((P/n).read_text())
def main():
 start=time.monotonic();local=read('local-final-cases.json');remote=read('cases-gpu-final.json');receipt=read('local-final-receipt.json')
 assert receipt['source_sha256']==sha('verify_common_mm_pairs.cpp')
 assert receipt['output_sha256']==sha('local-final-cases.json') and receipt['cpu_exit_code']==0
 remote_receipt=read('cpu-final-receipt.json')
 assert remote_receipt['compile']['source_sha256']==sha('verify_common_mm_pairs.cpp')
 assert remote_receipt['final_output']['sha256']==sha('cases-gpu-final.json')
 full=[r for r in remote_receipt['runs'] if r['mode']=='full']
 assert len(full)==1 and full[0]['exit_code']==0 and full[0]['case_count']==120
 assert full[0]['output_sha256']==sha('cases-gpu-final.json')
 assert {k:v for k,v in local.items() if k!='elapsed_seconds'}=={k:v for k,v in remote.items() if k!='elapsed_seconds'}
 assert remote['case_count']==len(remote['cases'])==120 and remote['seed']==20261120
 stats={}
 for row in remote['cases']:
  k=row['k'];q=1<<k;assert row['n']==2*k and row['N']==q*q and 2<=k<=6
  F,G=[sum(int(w,16)<<(64*i) for i,w in enumerate(row[n])) for n in ['f_truth_words_le','g_truth_words_le']]
  def forms(M):
   assert M.bit_length()<=q*q
   offsets=[(M>>(y<<k))&1 for y in range(q)]
   perm=[sum((((M>>((y<<k)+(1<<i)))&1)^offsets[y])<<i for i in range(k)) for y in range(q)]
   assert sorted(perm)==list(range(q))
   assert all(((M>>(x+(y<<k)))&1)==(((x&perm[y]).bit_count()&1)^offsets[y]) for y in range(q) for x in range(q))
   return perm,offsets
  a,h=forms(F);b,j=forms(G);changed=[y for y in range(q) if a[y]!=b[y]]
  full=[y for y in range(q) if a[y]==b[y] and h[y]!=j[y]]
  assert row['type']==f'c{len(changed)}_t{len(full)}'
  A=F^G;assert A.bit_count()==row['support_size']==2*q
  K={0}
  for v in row['kernel_basis_ambient']:K|={z^v for z in tuple(K)}
  assert len(K)==1<<(k-1)
  assert all(not(A>>(x^z)&1) for x in range(q*q) if A>>x&1 for z in K-{0})
  key=f'k{k}/{row["type"]}';v=stats.setdefault(key,{'samples':0,'changed_fibers':0,'nonzero_changed_offsets':0})
  v['samples']+=1;v['changed_fibers']+=len(changed);v['nonzero_changed_offsets']+=sum(h[y]^j[y] for y in changed)
 assert len(stats)==15 and all(v['samples']==8 for v in stats.values())
 assert all(v['nonzero_changed_offsets']>0 for key,v in stats.items() if 'c0' not in key)
 out={'status':'passed','scope':'120 saved MM pairs; no new full Walsh rerun','local_remote_cases_equal':True,'all_recovered_maps_permutations':True,'all_truth_tables_reconstructed':True,'all_explicit_kernels_direct_collision_free':True,'strata':stats,'python':platform.python_version(),'elapsed_seconds':time.monotonic()-start,'source_sha256':sha('verify_cross_dimension.py'),'remote_input_sha256':sha('cases-gpu-final.json'),'local_input_sha256':sha('local-final-cases.json')}
 (P/'cross-dimension-audit.json').write_text(json.dumps(out,indent=2)+'\n')
 print(json.dumps({k:out[k] for k in ['status','local_remote_cases_equal','all_recovered_maps_permutations','all_explicit_kernels_direct_collision_free','elapsed_seconds']}))
if __name__=='__main__':main()
