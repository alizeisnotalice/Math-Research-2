"""Verify MacBook worker receipts and aggregate threshold diagnostics."""
import json
from pathlib import Path

d=Path('thresholds')
workers=sorted([json.loads(p.read_text()) for p in d.glob('shard-*.json')],key=lambda r:r['shard'])
exits=[list(map(int,x.split())) for x in (d/'exits.tsv').read_text().splitlines()]
assert exits==[[r['shard'],r['pid'],0] for r in workers] and len(workers)==15
assert len(Path('remote-final-ps.txt').read_text().splitlines())==1
ps=[]
for p in sorted(d.glob('ps-*.txt')):
    v=[float(l.split()[2]) for l in p.read_text().splitlines()[1:]]
    ps.append(dict(workers=len(v),aggregate_cpu_percent=sum(v)))
rows=[r for w in workers for r in w['results']]
out=dict(utc=[(d/'start_utc').read_text().strip(),(d/'end_utc').read_text().strip()],
         worker_exits=exits,all_processes_absent=True,profiles=sum(r['points'] for r in rows),
         nested_iterations=sum(r['nested'] for r in rows),software=workers[0]['software'],
         cpu_seconds=sum(w['cpu'] for w in workers),
         actual_worker_elapsed=[min(w['elapsed'] for w in workers),max(w['elapsed'] for w in workers)],
         ps=ps,
         system_cpu_idle_samples=[x for x in (d/'top.txt').read_text().splitlines() if x.startswith('CPU usage:')],
         greatest_diagnostics={k:max(r[k] for r in rows) for k in rows[0] if k.startswith('max_')})
Path('summary.json').write_text(json.dumps(out,indent=2)+"\n")
print(out)
