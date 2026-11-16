"""Collect only completed, independently replayed support-layer batches.

This packages evidence; it does not replace the record-by-record verifier.
It refuses unfinished jobs and verifies source/ledger hashes and scope again.
"""
import argparse
from collections import Counter
import hashlib,json
from pathlib import Path
from audit_layer_coverage import audit

def collect(directory):
    base=Path(__file__).resolve().parent
    batch=json.loads((directory/'batch.json').read_text())
    assert batch['all_workers_exited'] and all(j['status']=='verified' and j['process_exited'] for j in batch['jobs'])
    report={'batches':[{'directory':str(directory),'elapsed_seconds':batch['elapsed_seconds'],
                      'workers':batch['workers'],'jobs':len(batch['jobs']),
                      'driver_sha256':batch.get('recovery_source_sha256',batch.get('driver_sha256')),
                      'driver_source': 'recover_support_batches.py' if 'recovery_source_sha256' in batch else batch.get('driver_source','not_recorded'),
                      'elapsed_scope':batch.get('elapsed_scope','original batch wall time'),
                      'all_workers_exited':True}],
            'source_sha256':{f:hashlib.sha256((base/f).read_bytes()).hexdigest() for f in ('support_search_layer.py','verify_support_layer.py')},'groups':[]}
    seen=set()
    for config in sorted(batch['configs'],key=lambda c:(c['max_terms'],c['total'])):
        jobs=[];counts=Counter()
        expected=set(range(1,config['total']//2+1))
        recorded={j['m']for j in batch['jobs']if j['tag']==config['tag']}
        assert expected==recorded
        for m in sorted(expected):
            p=directory/f"{config['tag']}-m{m:02d}"
            state=json.loads((p/'job.json').read_text());assert state['status']=='verified' and state['process_exited']
            summary=json.loads((p/'summary.json').read_text());verification=json.loads((p/'verification.json').read_text())
            unknown=[json.loads(s)for s in (p/'unknown.jsonl').read_text().splitlines()if s.strip()]
            for c in unknown:
                key=(tuple(c['I']),tuple(c['J']));assert key not in seen;seen.add(key)
            assert summary['complete'] and verification['verified']
            assert summary['counts']==verification['counts']
            counts.update(summary['counts'])
            jobs.append({'directory':str(p),'summary':summary,'verification':verification,'unknown':unknown})
        report['groups'].append({k:config[k]for k in ('min_total','total','max_terms')}|{'counts':dict(counts),'jobs':jobs})
    return report

def main():
    ap=argparse.ArgumentParser();ap.add_argument('directory',type=Path);ap.add_argument('--output',type=Path,required=True);a=ap.parse_args()
    assert Path('/Volumes/data').is_mount()
    assert not a.output.exists(),'Refusing to overwrite evidence'
    data=collect(a.directory);a.output.write_text(json.dumps(data,indent=2)+'\n')
    print(json.dumps(audit(a.output),indent=2))
if __name__=='__main__':main()
