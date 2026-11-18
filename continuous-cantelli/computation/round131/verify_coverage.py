"""Recover exact interval unions and verify that every pair on j+k=948 is covered."""
from datetime import datetime,UTC
from decimal import Decimal as D
import json
from pathlib import Path

LEFT=D('2.291'); RIGHT=D('3.5'); MICRO=D('0.000001')


def rows(path):
    p=Path(path)
    xs=[json.loads(f.read_text()) for f in p.glob('shard-*.json')]
    assert len(xs)==15
    exits=[list(map(int,x.split())) for x in (p/'exits.tsv').read_text().splitlines()]
    assert len(exits)==15 and all(x[2]==0 for x in exits)
    assert {x['pid'] for x in xs}=={x[1] for x in exits}
    assert (p/'status.txt').read_text().strip()=='0'
    assert len((p/'post-exit-ps.txt').read_text().splitlines())==1
    return xs


def union(data):
    out=[]
    for a,b in sorted(data):
        if out and a<=out[-1][1]: out[-1][1]=max(b,out[-1][1]);continue
        out.append([a,b])
    return out


def main():
    coarse=[t for r in rows('broad') for t in r['results']]
    middle=[t for r in rows('simpson') for t in r['results']]
    dangerous=rows('cover')
    assert len(coarse)==949 and len({t['j'] for t in coarse})==949
    bboxes={t['j']:t['passed'] for t in coarse}
    mboxes={t['j']:t['covered'] for t in middle}
    dboxes={t['j']:t['boxes'] for t in dangerous}
    counts=[0,0,0]; full=[]; incomplete={}; broad_simpson_full=[]
    for j in range(949):
        spans=[]
        for a,b,q,_ in bboxes[j]:
            assert D(q)>0
            spans.append((D(a),D(b)));counts[0]+=1
        for a,b,q in mboxes.get(j,[]):
            assert D(q)>0
            spans.append((D(a)*MICRO,D(b)*MICRO));counts[1]+=1
        independent=union(spans)
        if any(a<=LEFT and b>=RIGHT for a,b in independent):
            broad_simpson_full.append(j)
        for a,b,q,valid in dboxes.get(j,[]):
            assert valid and D(q)>0
            spans.append((D(a),D(b)));counts[2]+=1
        spans=union(spans)
        current=LEFT;missing=[]
        for a,b in spans:
            if a>current:missing.append([str(current),str(a)])
            current=max(current,b)
        if current<RIGHT:missing.append([str(current),str(RIGHT)])
        if missing:incomplete[j]=missing
        else:full.append(j)
    performance={}
    for name,a,b in [('cover','started','ended'),('broad','start','end'),('simpson','start','end')]:
        r=rows(name);lo=max(x[a] for x in r);hi=min(x[b] for x in r)
        performance[name]={'common_start_utc':datetime.fromtimestamp(lo,UTC).isoformat(),
          'common_end_utc':datetime.fromtimestamp(hi,UTC).isoformat(),'common_seconds':hi-lo,
          'total_cpu_seconds':sum(x['cpu'] for x in r),'worker_pids':[x['pid'] for x in r]}
    out=dict(total=948,full_pairs=len(full),full_pairs_using_broad_and_simpson=len(broad_simpson_full),counts=counts,missing=incomplete,performance=performance)
    Path('summary.json').write_text(json.dumps(out,indent=2)+'\n')
    print(out['full_pairs'],out['full_pairs_using_broad_and_simpson'],out['counts'],len(incomplete),performance)

if __name__=='__main__':main()
