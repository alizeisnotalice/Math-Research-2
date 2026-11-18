"""Check both ledgers; replay a deterministic stratified subset with direct shoelace corners."""
from argparse import ArgumentParser
from collections import Counter
from fractions import Fraction as Q
from hashlib import sha256
import json
from math import floor,prod
from pathlib import Path

from flint import arb
from certify_three_body_refinement import corner_bounds
import certify_three_body_weighted as w


def audit(source,phase1,phase2,samples):
    raw=phase1.read_bytes();s=json.loads(raw);t=json.loads(phase2.read_bytes());threshold=Q(t['target'])
    assert s['target']==t['target'] and sha256(raw).hexdigest()==t['source_stage_sha256']
    assert Q(s['volume'])==prod(Q(b)-Q(a) for a,b in zip(*s['root']))==sum(Q(p['volume']) for p in s['parts'])
    c=s['counts'];assert c['checked']==2*(s['unresolved']+c['taylor']+c['corners'])-c['seed_tiles']
    jobs=[(p['part'],i,(tuple(map(Q,a)),tuple(map(Q,b)))) for p in s['parts'] for i,(a,b) in enumerate(p['unresolved'])]
    stats=Counter();lowest=None
    assert len(jobs)==len(t['parts'])==t['residual_input']==s['unresolved']
    for (part,i,_),row in zip(jobs,t['parts']):
        assert (part,i)==(row['part'],row['index']) and not row['unresolved']
        cc=row['counts'];assert cc['checked']==2*sum(v for k,v in cc.items() if k.startswith('weights_'))-1
        z=Q(row['min_bound']);assert z>threshold;lowest=z if lowest is None else min(lowest,z);stats.update(cc)
    assert dict(stats)==t['counts'] and str(lowest)==t['min_weighted_bound'] and t['completed'] and t['residual_output']==0
    def direct(runs,pa,h):
        row=[]
        for bound in corner_bounds(runs,pa,h):
            n=floor(float(bound.lower())*w.SCALE)-1;assert bound>arb(n)/w.SCALE;row.append(n)
        return row
    w.init_worker(source,s['target'],t['max_depth']);w.corner_row=direct
    assert 2<=samples<=len(jobs)
    indices=sorted({j*(len(jobs)-1)//(samples-1) for j in range(samples)})
    subset_counts=Counter()
    for i in indices:
        record=w.solve(jobs[i]);assert record==t['parts'][i];subset_counts.update(record['counts'])
    return {'phase1_sha256':sha256(raw).hexdigest(),'phase2_sha256':sha256(phase2.read_bytes()).hexdigest(),
            'target':str(threshold),'all_ledger_checks':True,'direct_corner_indices':indices,
            'direct_replay_counts':dict(subset_counts),'minimum_second_stage_lower':str(lowest),'passed':True}

if __name__=='__main__':
    p=ArgumentParser();p.add_argument('source');p.add_argument('phase1',type=Path);p.add_argument('phase2',type=Path);p.add_argument('--samples',type=int,default=48);a=p.parse_args();print(json.dumps(audit(a.source,a.phase1,a.phase2,a.samples),indent=2))
