"""Summarize authoritative MacBook worker receipts and quadrature checks."""
import json
from pathlib import Path


def process(name):
    path = Path(name)
    files = [json.loads(p.read_text()) for p in path.glob('shard-*.json')]
    exits = [list(map(int, line.split())) for line in (path/'exits.tsv').read_text().splitlines()]
    assert len(files) == len(exits) == 15
    assert {row[1] for row in exits} == {f['pid'] for f in files}
    assert all(row[2] == 0 for row in exits) and (path/'status.txt').read_text().strip() == '0'
    assert len((path/'post-exit-ps.txt').read_text().splitlines()) == 1
    assert 'ps_status=1' in Path('remote-final-ps.txt').read_text()
    samples=[]
    for block in (path/'timeline.txt').read_text().split('2026-')[1:]:
        lines=block.splitlines()
        rows=[v.split(maxsplit=6) for v in lines[2:]]
        samples.append({'utc':'2026-'+lines[0], 'runnable':sum(v[3].startswith('R') for v in rows),
                        'total_cpu_percent':sum(float(v[2]) for v in rows)})
    models={}
    for model in {r['model'] for f in files for r in f['results']}:
        all_rows=[r for f in files for r in f['results'] if r['model']==model]
        candidates=[entry for r in all_rows for entry in r['tested_negative'] if entry[4] < 0]
        models[model]={'profiles':sum(r['profiles'] for r in all_rows),
                       'y_grid':all_rows[0]['y_range'], 'y_count':all_rows[0]['y_count'],
                       'stored_negative_candidates':len(candidates),
                       'first_by_j':min(candidates,default=None),
                       'first_by_k':min(candidates,key=lambda r:(r[1],r[0]),default=None),
                       'first_by_j_plus_k':min(candidates,key=lambda r:(r[0]+r[1],r[0]),default=None),
                       'strongest':min(candidates,key=lambda r:r[4],default=None),
                       'max_relative_recheck_difference':max((abs(row[3]-row[4])/abs(row[4]) for row in candidates),default=None)}
    a=min(f['start_epoch'] for f in files);b=max(f['end_epoch'] for f in files)
    common_start=max(f['start_epoch'] for f in files);common_end=min(f['end_epoch'] for f in files)
    return {'start_utc':(path/'start_utc').read_text().strip(), 'end_utc':(path/'end_utc').read_text().strip(),
            'supervisor':int((path/'supervisor.pid').read_text()), 'worker_pids':[r[1] for r in exits],
            'all_exits_zero':True, 'versions':files[0]['versions'],'wall_seconds':b-a,
            'worker_common_span_seconds':common_end-common_start,
            'total_cpu_seconds':sum(f['cpu_seconds'] for f in files), 'samples':samples, 'models':models}


if __name__ == '__main__':
    record={s:process(s) for s in ('scan','deep-final')}
    Path('summary.json').write_text(json.dumps(record,indent=2)+'\n')
    for s, data in record.items():
        print(s, data['start_utc'], data['end_utc'],data['wall_seconds'],data['total_cpu_seconds'],data['models'])
