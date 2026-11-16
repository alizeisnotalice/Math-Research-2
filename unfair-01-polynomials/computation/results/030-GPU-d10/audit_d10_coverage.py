#!/usr/bin/env python3
"""Independent deterministic coverage ledger for the completed d=10 run."""
import argparse, hashlib, itertools, json
from pathlib import Path

PAIR=((0,0),(1,0),(0,1))
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
def reasons(code,d):
    a,b=decode(code); p=product(a,b,2*d-2); out=[]
    if any(x>1 for x in p): out.append('prefix_product_coefficient_gt_1')
    for i in range(1,d):
        baseline=p[d+i] if d+i<len(p) else 0
        if baseline>=1 and a[i]+b[i]>0: out.append('forced_positive_fractional_cross_at_d_plus_i'); break
    return out
def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--run-dir',type=Path,required=True); ap.add_argument('--worker',type=Path,required=True); ap.add_argument('--record',type=Path,required=True); x=ap.parse_args()
    cfg=json.loads((x.run_dir/'config.json').read_text()); cp=json.loads((x.run_dir/'checkpoint.json').read_text()); d=int(cfg['fractional_index']); starts=int(cfg['starts']); steps=int(cfg['steps']); every=int(cfg['restart_every'])
    valid=[''.join(c) for c in itertools.product('012',repeat=d-1) if legal(''.join(c),d)]; pool=[c for c in valid if not reasons(c,d)]
    digest=lambda z:hashlib.sha256(json.dumps(z,separators=(',',':')).encode()).hexdigest(); worker=x.worker.read_text(); worker_sha=hashlib.sha256(x.worker.read_bytes()).hexdigest()
    required='assigned = [codes[(row + offset) % len(codes)] for row in range(starts)]'
    if worker_sha!=cfg['script_sha256'] or required not in worker: raise SystemExit('worker hash or deterministic assignment implementation mismatch')
    restarts=sorted(x.run_dir.glob('restart_*.json')); expected=[x.run_dir/f'restart_{s:07d}.json' for s in range(every,steps,every)]
    if restarts!=expected or not cp.get('history') or cp['history'][-1].get('step')!=steps: raise SystemExit('incomplete restart/final ledger')
    epochs=list(range(len(restarts)+1)); epoch_rows=[]
    for epoch in epochs:
        counts={code:0 for code in pool}
        for row in range(starts): counts[pool[(row+epoch)%len(pool)]]+=1
        values=list(counts.values())
        epoch_rows.append({'epoch':epoch,'assigned_rows':sum(values),'all_prefixes_assigned':min(values)>=1,'minimum_assignments_per_prefix':min(values),'maximum_assignments_per_prefix':max(values),'prefixes_at_minimum':values.count(min(values)),'prefixes_at_maximum':values.count(max(values))})
    out={'status':'deterministic_prefix_coverage_ledger_passed','claim_boundary':'This verifies the recorded finite-prefix assignment schedule from frozen source and saved completed-cohort ledger; it is not a statement about numerical optimization success or global factorization.','parameters':{'d':d,'starts':starts,'steps':steps,'restart_every':every,'cohorts':len(epochs)},'source':{'worker_sha256':worker_sha,'prefix_assignment_expression_verified':required},'pool':{'valid_low_prefix_count':len(valid),'valid_sha256':digest(valid),'survivor_count':len(pool),'survivor_sha256':digest(pool)},'cohorts':epoch_rows,'total_assignment_rows_by_schedule':starts*len(epochs),'complete_restart_snapshots':len(restarts),'final_step_completed':cp['history'][-1]['step']==steps}
    x.record.parent.mkdir(parents=True,exist_ok=True); x.record.write_text(json.dumps(out,indent=2,sort_keys=True)+'\n'); print(json.dumps(out['cohorts'],sort_keys=True))
if __name__=='__main__': main()
