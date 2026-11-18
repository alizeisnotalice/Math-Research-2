"""Aggregate exact Hermite cap profiles and verified worker exits."""
import json
from pathlib import Path
from fractions import Fraction

d=Path('hermite')
w=sorted([json.loads(p.read_text()) for p in d.glob('shard-*.json')],key=lambda x:x['shard'])
exits=[list(map(int,l.split())) for l in (d/'exits.tsv').read_text().splitlines()]
assert exits==[[x['shard'],x['pid'],0] for x in w] and len(w)==15
assert len(Path('remote-final-ps.txt').read_text().splitlines())==1
rows=[r for x in w for r in x['results']]
exact=[r for x in w for r in x['exact']]
ps=[]
for path in sorted(d.glob('ps-*.txt')):
    v=[float(l.split()[2]) for l in path.read_text().splitlines()[1:]]
    ps.append(dict(workers=len(v),cpu_percent=sum(v)))
out=dict(utc=[(d/'start_utc').read_text().strip(),(d/'end_utc').read_text().strip()],exits=exits,
    related_processes_absent=True,profiles=sum(r['count'] for r in rows),exact_profiles=len(exact),
    exact_max_energy_width=str(max(Fraction(r['upper'])-Fraction(r['lower']) for r in exact)),
    exact_max_cutoff=max(r['cutoff'] for r in exact),software=w[0]['software'],
    aggregate_cpu_seconds=sum(r['cpu'] for r in w),ps_samples=ps,
    worker_elapsed_range=[min(r['elapsed'] for r in w),max(r['elapsed'] for r in w)],
    cpu_samples=[l for l in (d/'top.txt').read_text().splitlines() if l.startswith('CPU usage:')],
    max_PMF_energy_disagreement=max(r['poly_prec'] for r in rows))
Path('summary.json').write_text(json.dumps(out,indent=2)+"\n")
print(out)
