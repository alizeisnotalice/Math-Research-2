#!/usr/bin/env python3
"""Independent received-range audit; skips assert prior work, never add work."""
import argparse, hashlib, json
from pathlib import Path

def main():
    p=argparse.ArgumentParser();p.add_argument('manifest',type=Path);p.add_argument('--inputs',type=Path,required=True);p.add_argument('--out',type=Path,required=True);a=p.parse_args()
    raw=a.inputs.read_bytes();data=json.loads(raw);m=json.loads(a.manifest.read_text())
    assert hashlib.sha256(raw).hexdigest()==m['source_inputs_sha256']
    assert [v['name'] for v in m['layers']]==['v1','v2','v3']
    spans={i:[] for i in range(len(data['tasks']))};skips=[];stats=[]
    for layer in m['layers']:
        b=Path(layer['path']).read_bytes();assert hashlib.sha256(b).hexdigest()==layer['sha256']
        count=mask_count=0
        for line in b.decode('utf8',errors='strict').splitlines():
            if not line.strip().startswith('{'):
                assert '"event"' not in line and '"batch"' not in line;continue
            e=json.loads(line)
            if e.get('event')=='batch':
                t,lo,hi,n=(e[k] for k in ['task','first','last','count'])
                assert all(type(v) is int for v in [t,lo,hi,n])
                assert t in spans and 0<=lo<=hi<data['tasks'][t]['masks'] and hi-lo+1==n
                spans[t].append((lo,hi+1));count+=1;mask_count+=n
            elif e.get('event')=='task_skipped_completed':
                t=e['task'];assert e['start']==e['masks']==data['tasks'][t]['masks'];skips.append(t)
        stats.append(dict(name=layer['name'],batches=count,masks=mask_count,sha256=layer['sha256']))
    nexts={}
    for t,rows in spans.items():
        end=0
        for lo,hi in sorted(rows): assert lo==end,(t,lo,end);end=hi
        nexts[str(t)]=end
    for t in skips:assert nexts[str(t)]==data['tasks'][t]['masks']
    assert nexts==m['task_next']
    assert stats[0]['masks']==157384704 and stats[1]['masks']==44957696
    out=dict(status='passed_independent_received_range_union',layers=stats,task_next=nexts,total_received_masks=sum(s['masks'] for s in stats),manifest_sha256=hashlib.sha256(a.manifest.read_bytes()).hexdigest(),source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),scope='Received numerical masks only; not correctness of every score or mathematical exclusion.')
    a.out.write_text(json.dumps(out,indent=2)+'\n'); print(json.dumps(out,sort_keys=True))
if __name__=='__main__':main()
