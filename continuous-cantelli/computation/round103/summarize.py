"""Aggregate multiregion test receipts and verify all worker exits."""
import json
from pathlib import Path

d=Path('wells')
w=sorted([json.loads(p.read_text()) for p in d.glob('shard-*.json')],key=lambda r:r['shard'])
exits=[list(map(int,l.split())) for l in (d/'exits.tsv').read_text().splitlines()]
assert len(w)==15 and exits==[[r['shard'],r['pid'],0] for r in w]
assert len(Path('remote-final-ps.txt').read_text().splitlines())==1
rows=[x for rec in w for x in rec['results']]
top=[rec['top'] for rec in w if rec['top']]
ps=[]
for p in sorted(d.glob('ps-*.txt')):
 v=[float(l.split()[2]) for l in p.read_text().splitlines()[1:]]
 ps.append(dict(workers=len(v),aggregate_cpu_percent=sum(v)))
out=dict(utc=[(d/'start_utc').read_text().strip(),(d/'end_utc').read_text().strip()],worker_exits=exits,
        processes_absent=True,software=w[0]['software'],points=sum(x['points'] for x in rows),
        feasible_budget_count=sum(x['feasible_budget'] for x in rows),
        joint_only_screened_count=sum(x['joint_witness_count'] for x in rows),
        max_joint_ratio=max(x['maximum_joint_ratio'] for x in rows),
        largest_selected_direct_quadrature_difference=max(x['independent_quad_error'] for x in top),
        best_profile=max(top,key=lambda r:r['ratio']),
        independently_integrated_best=len(top),aggregate_cpu_seconds=sum(x['cpu'] for x in w),
        elapsed_worker_bounds=[min(x['elapsed'] for x in w),max(x['elapsed'] for x in w)],ps=ps,
        system_idle_samples=[l for l in (d/'top.txt').read_text().splitlines() if l.startswith('CPU usage:')])
Path('summary.json').write_text(json.dumps(out,indent=2)+"\n")
print({k:v for k,v in out.items() if k!='worker_exits'})
