"""Recompute every pure five-body terminal and every contraction from fresh area geometry."""
import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
from hashlib import sha256
import json
import os
from pathlib import Path
import subprocess
import sys
import tarfile


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for name in ('packages','upstream','tree','cache','work','output'):
        p.add_argument(name,type=Path)
    p.add_argument('--workers',type=int,default=8)
    a=p.parse_args(); assert __debug__ and 1<=a.workers<=32
    a.work.mkdir(parents=True,exist_ok=True)
    code=a.work/'src'; code.mkdir(exist_ok=True)
    expected={}; sources={}
    for label in ('a','b'):
        with tarfile.open(a.packages/f'native45-pure-{label}.tar.gz') as tar:
            for name,size,digest in json.load(tar.extractfile('MANIFEST.json')):
                if not name.startswith('pure/'):continue
                data=tar.extractfile(name).read()
                assert len(data)==size and sha256(data).hexdigest()==digest
                if name.startswith('pure/actual/'):
                    if name.endswith('.json'):
                        row=json.loads(data);k=(row['part'],row['pos']);assert k not in expected
                        expected[k]=(row,digest)
                elif name.endswith(('.py','.cpp')):
                    short=Path(name).name;assert short not in sources
                    (code/short).write_bytes(data);sources[short]=digest
    assert len(expected)==196
    assert all(all(sources[name]==h for name,h in row['sources'].items()) for row,digest in expected.values())
    from collections import Counter
    totals=Counter()
    # Longest recorded work first; the original roots remain disjoint.
    jobs=sorted(expected,key=lambda key:-expected[key][0]['result']['vis'])
    def one(key):
        row,expected_digest=expected[key];part,pos=key
        assert row['threshold']=='0.837455' and row['budget']==32768 and not row['result']['remaining']
        dest=a.work/f'part{part}-pos{pos}'; dest.mkdir(exist_ok=True)
        record=dest/f'part{part}-pos{pos}.json'; trace=record.with_suffix('.bin')
        if not record.exists():
            cmd=[sys.executable,'-I','-B',code/'replay_five_affine_selected.py',a.upstream,a.tree,dest,a.cache,'--part',str(part),'--positions',str(pos),'--target',row['threshold'],'--budget',str(row['budget'])]
            with (dest/'stdout.txt').open('w') as stdout:
                subprocess.run(cmd,stdout=stdout,check=True)
        fresh=json.loads(record.read_bytes());assert fresh['root']==row['root'] and not fresh['result']['remaining']
        assert sha256(trace.read_bytes()).hexdigest()==fresh['trace_sha256']
        from verify_five_selected_083740_traces import audit
        z=audit(fresh['root'],trace.read_bytes()); assert z=={k:fresh['result'][k] for k in z}
        return (part,pos,z['vis'],sum(z['closed'].values()),sha256(record.read_bytes()).hexdigest()==expected_digest,fresh['trace_sha256']==row['trace_sha256'])
    sys.path.insert(0,str(code.resolve()))
    rows=[]
    with ThreadPoolExecutor(max_workers=a.workers) as pool:
        for future in as_completed([pool.submit(one,k) for k in jobs]):
            rows.append(future.result())
            if len(rows)%16==0: print(len(rows),sum(r[2] for r in rows),'nodes',flush=True)
    out=dict(complete=True,threshold='0.837455', original_leaves=len(rows),replayed_geometric_visits=sum(r[2] for r in rows),replayed_terminals=sum(r[3] for r in rows),original_record_matches=sum(r[4] for r in rows),original_trace_matches=sum(r[5] for r in rows),sources=sources,rows=sorted(rows))
    a.output.write_text(json.dumps(out,indent=2)+'\n');print('complete',out['replayed_geometric_visits'],flush=True)

if __name__=='__main__':main()
