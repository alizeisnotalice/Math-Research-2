#!/usr/bin/env python3
"""Read-only numerical ledger postprocess; no tree search or classification."""
from fractions import Fraction as F
from pathlib import Path
from collections import Counter
import gzip,hashlib,json,os,sys,time
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def absent(pid):
 try:os.kill(pid,0);return False
 except ProcessLookupError:return True
def main():
 t=time.monotonic();root=Path(__file__).resolve().parents[2];rp=root/'computation/results/round079-quad';prod=rp/'production';sp=prod/'summary.json';summary=json.loads(sp.read_text());receipt=json.loads((rp/'process-receipt.json').read_text());assert receipt['exit_code']==0
 actuals={'corrected':('computation/results/round078-actual-prefix/actual-prefix-area.json','0403e1e720e2802048bb324ccd0cdb730f63eeb55fed881d5a53fcdf2b9b961c'),'farface':('computation/results/round078-farface/farface-area.json','4d649e564bdf22de9d70bfa92cbfe0291fd8ee560d7d095261ced9843c537fca')}
 cases=[]
 for case in summary['cases']:
  cid=case['case_id'];ap,expected=actuals[cid];assert sha(root/ap)==expected;actual=json.loads((root/ap).read_text());act=actual['actual_Gamma_deletion_area'];assert actual['input_sha256']==case['input_sha256']
  cp=prod/case['certificate'];assert sha(cp)==case['certificate_sha256'];cert=json.loads(gzip.decompress(cp.read_bytes()));variants={}
  for name,v in case['variants'].items():
   lo,hi=F(v['envelope_minus_true_residual_lower']),F(v['envelope_minus_true_residual_upper']);dl=lo-F(act['upper']);du=hi-F(act['lower'])
   variants[name]={'state_interval':[str(lo),str(hi)],'state_decimal_interval':[float(lo),float(hi)],'state_gap':str(hi-lo),'state_gap_decimal':float(hi-lo),'target_met':hi-lo<=F(1,10000),'state_upper_below_half_P':hi<F(1,8),'state_upper_below_main_corner_threshold':hi<F(80003,640000),'state_minus_actual_interval':[str(dl),str(du)],'state_minus_actual_decimal_interval':[float(dl),float(du)],'actual_area_is_from_frozen_R78':True,'state_lower_is_not_actual_lower':True}
  tilebrief=[{k:z[k] for k in ['tile_id','pid','termination','created_node_count','split_count','leaf_count','maximum_depth','gaps','unresolved_leaf_area']} for z in cert['tiles']]
  cases.append({'case_id':cid,'actual_input_sha256':case['input_sha256'],'actual_R78_path':ap,'actual_R78_sha256':expected,'actual_R78_area':act,'variants':variants,'tile_termination_counts':dict(Counter(z['termination'] for z in cert['tiles'])),'tiles':tilebrief,'certificate_sha256':sha(cp),'certificate_bytes':cp.stat().st_size,'certificate_leaf_count':sum(z['leaf_count'] for z in cert['tiles']),'coverage_replay_scope':'Counts/totals are read from production certificate, not independently reclassified here.'})
 pids=[summary['pid']]+summary['worker_pids'];checks=[{'pid':p,'observed_absent':absent(p)} for p in pids];assert all(z['observed_absent'] for z in checks)
 out={'schema':'r079-production-comparison-pending-audit-v1','status':'production complete; pending independent theorem and certificate audit','pid':os.getpid(),'python':sys.version,'script_sha256':sha(__file__),'summary_sha256':sha(sp),'production_elapsed_seconds':receipt['elapsed_seconds'],'production_script_seconds':summary['seconds'],'worker_limit':summary['workers'],'distinct_task_worker_pids':summary['worker_pids'],'distinct_task_worker_count':len(summary['worker_pids']),'pool_context_joined_before_summary':True,'pool_join_fact_source':'ProcessPoolExecutor context in the pinned production source; executed-task PIDs are additionally checked absent. Idle pool worker PIDs were not logged, so no list of eight observed PIDs is claimed.','observed_exits':checks,'cases':cases,'new_parameter_search_or_tree_runs_in_postprocess':0,'seconds':time.monotonic()-t,'scope':'Two fixed continuous-envelope-state brackets compared with inherited audited actual regional areas. State lower brackets do not give actual halo lower bounds. Corrected target not met; no second batch launched.'}
 p=rp/'comparison-pending-audit.json';assert not p.exists();p.write_text(json.dumps(out,indent=2)+'\n');print(json.dumps({'path':str(p),'sha256':sha(p),'production_elapsed_seconds':receipt['elapsed_seconds'],'observed_exits':checks,'cases':[{k:z[k] for k in ['case_id','tile_termination_counts','certificate_leaf_count','variants']} for z in cases]},indent=2))
if __name__=='__main__':main()
