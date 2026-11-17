#!/usr/bin/env python3
"""Independent total47 coverage/content binding; not a second semantic replay."""
import argparse, gzip, hashlib, itertools, json, math, time
from pathlib import Path
from collections import Counter

def sha(p):
    h=hashlib.sha256()
    with p.open('rb') as f:
        for b in iter(lambda:f.read(1<<20),b''): h.update(b)
    return h.hexdigest()

def main():
    ap=argparse.ArgumentParser()
    for key in ['run','plan','semantic-report','frozen','remote','audited','out']:
        ap.add_argument('--'+key,type=Path,required=True)
    a=ap.parse_args(); start=time.monotonic()
    r=json.loads(a.semantic_report.read_text()); manifest=r['manifest']
    assert r['status']=='passed_complete_total_47' and not manifest['errors']
    for key,p in [('events',a.run/'events.jsonl'),('plan',a.plan),('frozen',a.frozen),('remote',a.remote),('payload',a.run/'payload.py'),('receipt',a.run/'exit-receipt.json'),('audited',a.audited)]:
        assert sha(p)==manifest[key+'_sha256'],key
    wanted={}
    for m in range(1,24):
        index=0
        for count in range(5):
            for inner in itertools.combinations(range(1,m),count):
                wanted[(m,47-m,index)]=(0,*inner,m); index+=1
    plan=json.loads(a.plan.read_text()); planned={}
    for x in plan['jobs']:
        if x['total']!=47:continue
        k=x['m'],x['n'],x['index'];assert k not in planned;planned[k]=tuple(x['I'])
    assert planned==wanted and len(wanted)==44551
    results={v['job']:v for v in r['chunk_results']}
    assert len(results)==len(r['chunk_results'])==44551
    seen=set(); counts=Counter(); compressed=unpacked=lines=0; bind=hashlib.sha256()
    for text in (a.run/'events.jsonl').read_text().splitlines():
        x=json.loads(text)
        if x.get('event')!='chunk' or x['total']!=47:continue
        k=x['m'],x['n'],x['index'];assert k in wanted and k not in seen
        seen.add(k);assert tuple(x['I'])==wanted[k]
        job=f"m{x['m']:02d}n{x['n']:02d}/i{x['index']:06d}";v=results[job]
        assert not v['errors'] and v['counts']==x['counts']
        assert x['counts']['raw']==sum(math.comb(x['n']-1,j) for j in range(min(4,x['n']-1)+1))
        p=a.run/'t47'/(job+'.jsonl.gz'); assert sha(p)==x['gzip_sha256']==v['gzip_sha256']
        h=hashlib.sha256(); n=sz=0
        with gzip.open(p,'rb') as f:
            for b in iter(lambda:f.read(1<<20),b''):h.update(b);n+=b.count(b'\n');sz+=len(b)
        assert h.hexdigest()==x['raw_sha256']==v['raw_sha256']
        assert n==v['verified_rows']==x['counts']['canonical']
        assert n==sum(x['counts'][key] for key in ['excluded','fair','unknown'])
        counts.update(x['counts']);compressed+=p.stat().st_size;unpacked+=sz;lines+=n
        bind.update(json.dumps([job,x['gzip_sha256'],x['raw_sha256'],x['counts']],sort_keys=True).encode()+b'\n')
    assert seen==set(wanted) and dict(counts)==r['counts'] and counts['unknown']==0
    out=dict(status='passed_parent_total47_coverage_and_content_binding',chunks=len(seen),counts=dict(counts),packed_bytes=compressed,raw_bytes=unpacked,lines=lines,binding_digest=bind.hexdigest(),semantic_report_sha256=sha(a.semantic_report),source_sha256=sha(Path(__file__)),elapsed_seconds=time.monotonic()-start,scope='Total47 cap6 only; whole byte and coverage verification bound to separate completed independent semantic replay. No total48 or arbitrary-degree claim.')
    a.out.write_text(json.dumps(out,indent=2,sort_keys=True)+'\n'); print(json.dumps(out,sort_keys=True))

if __name__=='__main__':main()
