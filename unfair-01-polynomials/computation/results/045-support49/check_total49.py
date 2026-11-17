#!/usr/bin/env python3
"""Independent union, exact coverage, and complete byte binding for total 49.

The separately frozen verifier performs semantic replay. This program does not
import its enumeration or classifiers and does not infer arbitrary-degree claims.
"""
import argparse, gzip, hashlib, itertools, json, math, time
from collections import Counter
from pathlib import Path

def sha(p):
    h = hashlib.sha256()
    with Path(p).open('rb') as f:
        for b in iter(lambda: f.read(1 << 20), b''): h.update(b)
    return h.hexdigest()

def key(e): return e['total'], e['m'], e['n'], e['index']

def main():
    ap = argparse.ArgumentParser()
    for n in ('binding', 'report', 'out'): ap.add_argument('--'+n, type=Path, required=True)
    a=ap.parse_args(); start=time.monotonic(); a.out.mkdir(parents=True, exist_ok=True)
    binding=json.loads(a.binding.read_text()); paths={k:Path(v) for k,v in binding['paths'].items()}
    for k,p in paths.items(): assert sha(p)==binding['sha256'][k], k
    report_hash=sha(a.report); report=json.loads(a.report.read_text())
    assert report['total']==49 and report['errors']==0
    expected={}
    for total in (49,50):
        for m in range(1,total//2+1):
            index=0
            for size in range(2,min(m+1,6)+1):
                for inside in itertools.combinations(range(1,m),size-2):
                    expected[(total,m,total-m,index)]=(0,*inside,m); index+=1
    plan=json.loads(paths['fresh_plan'].read_text())['jobs']
    assert len(plan)==len(expected)==123859
    assert {key(e):tuple(e['I']) for e in plan}==expected
    events={}; side_counts=Counter()
    for side in ('old','resume'):
        for line in paths[side+'_events'].read_text().splitlines():
            e=json.loads(line)
            if e.get('event')!='chunk': continue
            k=key(e); assert k in expected and k not in events
            assert tuple(e['I'])==expected[k]
            events[k]=(side,e); side_counts[side]+=1
    assert set(events)==set(expected)
    assert side_counts=={'old':58036,'resume':65823}
    resume_plan=json.loads(paths['resume_plan'].read_text())['jobs']
    assert len(resume_plan)==65823
    assert {key(e):tuple(e['I']) for e in resume_plan}=={k:expected[k] for k,(s,e) in events.items() if s=='resume'}
    receipt=json.loads(paths['resume_receipt'].read_text())
    assert receipt['returncode']==0 and receipt['selected_set_complete'] and receipt['all_owned_processes_reaped']
    assert receipt['chunks_persisted']==65823
    results={(v['source'],v['job']):v for v in report['results']}
    assert len(results)==len(report['results'])==report['jobs']==55454
    counts=Counter(); pairs={}; packed_bytes=raw_bytes=lines=chunks=0
    with gzip.open(a.out/'chunk-bindings.jsonl.gz','wt',encoding='utf8') as output:
        for k,(side,e) in sorted(events.items()):
            total,m,n,index=k
            if total!=49: continue
            label=f'm{m:02d}n{n:02d}/i{index:06d}'
            v=results.pop((side,label)); assert not v['errors']
            assert v['counts']==e['counts']
            assert e['counts']['raw']==sum(math.comb(n-1,j) for j in range(min(4,n-1)+1))
            assert v['verified_rows']==e['counts']['canonical']
            assert e['counts']['canonical']==sum(e['counts'][x] for x in ('fair','excluded','unknown'))
            assert e['counts']['unknown']==0 and e['unknown']==[]
            p=paths[side+'_events'].parent/'t49'/(label+'.jsonl.gz')
            assert sha(p)==e['gzip_sha256']
            h=hashlib.sha256(); nlines=nbytes=0
            with gzip.open(p,'rb') as stream:
                for b in iter(lambda:stream.read(1<<20),b''):
                    h.update(b);nlines+=b.count(b'\n');nbytes+=len(b)
            assert h.hexdigest()==e['raw_sha256']
            assert nlines==v['verified_rows']
            counts.update(e['counts']); packed_bytes+=p.stat().st_size;raw_bytes+=nbytes;lines+=nlines;chunks+=1
            pair=pairs.setdefault(f'{m}+{n}',{'chunks':0,'counts':Counter()})
            pair['chunks']+=1;pair['counts'].update(e['counts'])
            output.write(json.dumps({'source':side,'job':label,'I':e['I'],'gzip_sha256':e['gzip_sha256'],'raw_sha256':e['raw_sha256'],'counts':e['counts']},sort_keys=True)+'\n')
            if chunks%10000==0:print(json.dumps({'checked_chunks':chunks,'raw_bytes':raw_bytes}),flush=True)
    assert not results and dict(counts)==report['counts'] and chunks==55454
    for k,p in paths.items(): assert sha(p)==binding['sha256'][k], ('changed_during_check',k)
    assert sha(a.report)==report_hash
    out={'status':'passed_parent_total49_coverage_and_content_binding','chunks':chunks,'counts':dict(counts),'packed_bytes':packed_bytes,'raw_bytes':raw_bytes,'lines':lines,'binding_sha256':sha(a.binding),'semantic_report_sha256':report_hash,'source_sha256':sha(Path(__file__)),'archive_sha256':sha(a.out/'chunk-bindings.jsonl.gz'),'elapsed_seconds':time.monotonic()-start,'scope':'Total49 cap6 only; bound to complete independent semantic replay. No total50 or arbitrary-degree claim.'}
    (a.out/'degree-pairs.json').write_text(json.dumps(pairs,indent=2,sort_keys=True)+'\n')
    (a.out/'parent-summary.json').write_text(json.dumps(out,indent=2,sort_keys=True)+'\n')
    print(json.dumps(out,sort_keys=True),flush=True)

if __name__=='__main__': main()
