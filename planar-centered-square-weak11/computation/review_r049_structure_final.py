#!/usr/bin/env python3
"""R49 final source binding and tiny boundary-ledger checks, not source rerun."""
from fractions import Fraction as F
from pathlib import Path
import argparse,hashlib,json,os,platform,sys,time
ROOT=Path(__file__).resolve().parents[1]
EXPECTED={
'proofs/sol-effective-separator-recursion.md':'a913f6eb302be63208214a92ae858be7e13a95d4561b1eaf5b228a6149bb95b5',
'computation/sol_effective_separator_check.py':'260908fc99bb3415aa0757f634907bffe3acaa5818e2eda58e0217dea1dc08c0',
'computation/results/sol-effective-separator-exact.json':'6350998ebd1c149ae3563750b5f9071f53332bf2531738cc6b705177b00168e8'}
def digest(path):return hashlib.sha256(path.read_bytes()).hexdigest()
def main():
 start=time.perf_counter();assert os.path.ismount('/Volumes/data')
 parser=argparse.ArgumentParser();parser.add_argument('--author-root',type=Path,default=ROOT);parser.add_argument('--main-topic-root',type=Path)
 args=parser.parse_args();AUTHOR=args.author_root.resolve()
 before={p:digest(AUTHOR/p) for p in EXPECTED};assert before==EXPECTED
 source=json.loads((AUTHOR/'computation/results/sol-effective-separator-exact.json').read_text())
 for rel,want in source['source_sha256'].items():
  assert digest(AUTHOR/rel)==want
  if args.main_topic_root is not None:assert digest(args.main_topic_root/rel)==want
 boundary_count=0
 for a in (F(0),F(1,4),F(1,2)):
  coords=sorted({F(0),a,F(1,2),1-a,F(1)})
  for x in coords:
   for y in coords:
    h=int(x<a)+int(x>1-a);v=int(y<a)+int(y>1-a)
    center=int(a<=x<=1-a and a<=y<=1-a)
    corner=h*v
    assert h in (0,1) and v in (0,1)
    assert h+v==1+corner-center
    boundary_count+=1
 K_checks=[]
 for K in (F(7,2),F(18,5),F(15,4),F(399,100)):
  a=(4-K)/2;R=(K-3)**2
  assert F(0)<=a<=F(1,2) and K-3>0
  assert a+(K-3)/2==F(1,2)
  assert (2-a)*2==K
  K_checks.append({'K':str(K),'a_over_L':str(a),'R_over_M':str(R),'necessary_corner_minus_center_over_M':str(4*R-1)})
 assert (F(49,50)-F(1,4)+1)*(F(49,50)+1)==F(17127,5000)
 own_prior=ROOT/'computation/results/post-r40-sorting/exact-obstructions.json'
 prior=json.loads(own_prior.read_text())
 assert prior['six_cycle_proper_rank_two_checks']==41 and prior['six_cycle_rank_root_proxy_upper']=='8/3'
 after={p:digest(AUTHOR/p) for p in EXPECTED};assert after==before
 result={'status':'PASS','python':sys.version,'platform':platform.platform(),
  'actual_scope':'Source hashes and 4 existing-main dependency hashes; 43 sampled boundary identities; 4 scalar K checks; one strip scalar identity. Author source script was read but NOT executed in this review.',
  'author_root_used':str(AUTHOR),'main_topic_dependency_root_checked':str(args.main_topic_root) if args.main_topic_root else None,
  'source_final_sha256':before,'author_result_proof_sha256_at_original_run':source['proof_sha256'],
  'boundary_point_checks':boundary_count,'K_scalar_checks':K_checks,
  'own_prior_result_reused_not_rerun':{'path':str(own_prior),'sha256':digest(own_prior)},
  'source_dependencies_already_in_main':source['source_sha256'],
  'elapsed_seconds':time.perf_counter()-start,
  'review_check_script_sha256':digest(Path(__file__)),
  'scope_limit':'Finite boundary checks support the analytic review; no general theorem is justified by these finite samples alone.'}
 out=ROOT/'computation/results/r049-structure-cross-review-check.json'
 out.write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n');print(json.dumps(result,indent=2))
if __name__=='__main__':main()
