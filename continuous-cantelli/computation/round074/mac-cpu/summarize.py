"""Aggregate MacBook Cantor coefficient diagnostics and verify exit records."""
import json
from pathlib import Path

d=Path('cantor')
w=sorted([json.loads(p.read_text()) for p in d.glob('shard-*.json')],key=lambda r:r['shard'])
r=[q for row in w for q in row['results']]
exits=[list(map(int,s.split())) for s in (d/'exits.tsv').read_text().splitlines()]
assert exits==[[x['shard'],x['pid'],0] for x in w] and len(w)==15
assert len(Path('remote-final-ps.txt').read_text().splitlines())==1
ps=[]
for p in sorted(d.glob('ps-*.txt')):
 vals=[float(l.split()[2]) for l in p.read_text().splitlines()[1:]]
 ps.append(dict(workers=len(vals),aggregate_cpu_percent=sum(vals)))
out=dict(utc=[(d/'start_utc').read_text().strip(),(d/'end_utc').read_text().strip()],node_exit=exits,
         profiles=sum(x['evaluated'] for x in r),nested_rechecks=sum(x['nested'] for x in r),
         aggregate_cpu_seconds=sum(x['cpu'] for x in w),software=w[0]['software'],
         worker_elapsed_range=[min(x['elapsed'] for x in w),max(x['elapsed'] for x in w)],
         node_ps=ps,idle_samples=[x for x in (d/'top.txt').read_text().splitlines() if x.startswith('CPU usage:')],
         max_diagnostics={k:max(x[k] for x in r) for k in r[0] if k.startswith('max_')},
         final_processes_present=0)
Path('summary.json').write_text(json.dumps(out,indent=2)+'\n')
print(out)
