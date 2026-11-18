"""Validate MacBook exits and preserve load/precision diagnostics."""
from datetime import datetime,timedelta
import json
from pathlib import Path

def load(p): return json.loads(Path(p).read_text())

d=Path('trace')
w=sorted((load(p) for p in d.glob('shard-*.json')),key=lambda r:r['shard'])
exits=[list(map(int,l.split())) for l in (d/'exits.tsv').read_text().splitlines()]
assert len(w)==15 and exits==[[r['shard'],r['pid'],0] for r in w]
assert len(Path('remote-final-ps.txt').read_text().splitlines())==1
ps=[]
for p in sorted(d.glob('ps-*.txt')):
 v=[float(l.split()[2]) for l in p.read_text().splitlines()[1:]]
 ps.append(dict(workers=len(v),cpu_percent=sum(v)))
rows=[r for x in w for r in x['results']]
m=load('remote-worker-timing.json')
skew=load('clock-comparison.json')['post_run_remote_minus_local_ns']
a=datetime.fromisoformat(m['all_worker_timers_start_upper'])-timedelta(microseconds=skew[0]/1000)
b=datetime.fromisoformat(m['all_worker_timers_end_lower'])-timedelta(microseconds=skew[1]/1000)
local_gpu_end=datetime.fromisoformat("2026-11-18T04:03:01.002+00:00")
audit=load('refine_high.json')
out=dict(utc=[(d/'start_utc').read_text().strip(),(d/'end_utc').read_text().strip()],
 exits=exits,related_processes_absent=True,profiles=sum(r['profiles'] for r in rows),nested_checks=sum(r['independent'] for r in rows),
 software=w[0]['software'],aggregate_cpu_seconds=sum(r['cpu'] for r in w),ps_samples=ps,
 worker_elapsed_range=[min(r['elapsed'] for r in w),max(r['elapsed'] for r in w)],
 cpu_samples=[l for l in (d/'top.txt').read_text().splitlines() if l.startswith('CPU usage:')],
 conservative_mac_timer_window_local_clock=[a.isoformat(),b.isoformat()],clock_offset_postrun_ns=skew,
 four_workflow_span_if_clock_stable_seconds=(local_gpu_end-a).total_seconds(),
 production_nested_max_absolute=max(r['max_nested'] for r in rows),
 refined_topten_max_absolute=max(x['max_fine_disagreement'] for x in audit['worst']),
 production_fixed_cusp_relative_not_reliable='The s=52 cutoff truncates leading mass at the deepest fixed-p cases; use the extended adaptive cutoff audit.')
Path('summary.json').write_text(json.dumps(out,indent=2)+"\n")
print({k:v for k,v in out.items() if k not in ('exits','ps_samples')})
