"""Bounded local runner for exact conditioned completion queries."""
from argparse import ArgumentParser
from concurrent.futures import ThreadPoolExecutor
import hashlib
import json
import os
from pathlib import Path
import re
import signal
import subprocess

HERE=Path(__file__).resolve().parent


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def timed(ix,work,out,timeout):
    src=out/f'part-{ix}.txt'
    raw=out/f'part-{ix}.jsonl'
    resource=out/f'part-{ix}.resource'
    src.write_text(''.join(f"{p['w']} {p['h']} {p['seen']} {p['states']}\n" for p in work))
    command=['/usr/bin/time','-l',str(HERE/'prefix_completion')]
    with src.open('r') as inp,raw.open('w') as stdout,resource.open('w') as stderr:
        proc=subprocess.Popen(command,stdin=inp,stdout=stdout,stderr=stderr,start_new_session=True)
        try:
            status=proc.wait(timeout=timeout)
        except subprocess.TimeoutExpired:
            os.killpg(proc.pid,signal.SIGTERM)
            try:
                proc.wait(timeout=3)
            except subprocess.TimeoutExpired:
                os.killpg(proc.pid,signal.SIGKILL)
                proc.wait()
            raise
    if status:
        raise RuntimeError(f'query part {ix} exited {status}; see {resource}')
    lines=[json.loads(line) for line in raw.read_text().splitlines()]
    assert len(lines)==len(work)
    for p,q in zip(work,lines):
        for field in ('w','h','seen','states','predecessor_lookups'):
            assert p[field]==q[field],(p,q,field)
        assert int(q['completion'])>0
    rss=re.search(r'^\s*(\d+)\s+maximum resident set size$',resource.read_text(),re.M)
    assert rss
    return lines,{'part':ix,'queries':len(lines),'pid':proc.pid,'exit_status':status,'max_rss_bytes':int(rss[1]),
                  'raw_sha256':digest(raw),'resource_sha256':digest(resource)}


def main():
    parser=ArgumentParser()
    parser.add_argument('plan',type=Path)
    parser.add_argument('output',type=Path)
    parser.add_argument('--workers',type=int,choices=(1,2,3),default=3)
    parser.add_argument('--timeout',type=int,default=120)
    args=parser.parse_args()
    data=json.loads(args.plan.read_text())
    work=[[] for _ in range(args.workers)]
    totals=[0]*len(work)
    for p in sorted(data['tasks'],key=lambda q:-q['predecessor_lookups']):
        i=min(range(len(work)),key=totals.__getitem__)
        work[i].append(p)
        totals[i]+=p['predecessor_lookups']+1
    args.output.mkdir(exist_ok=False)
    with ThreadPoolExecutor(max_workers=args.workers) as pool:
        results=list(pool.map(lambda p:timed(p[0],p[1],args.output,args.timeout),enumerate(work)))
    records=[r for part,_ in results for r in part]
    summary={'plan_sha256':digest(args.plan),'source_sha256':digest(HERE/'prefix_completion.cpp'),
             'queried':len(records),'independent_sets_processed':sum(q['states'] for q in records),
             'predecessor_lookups':sum(q['predecessor_lookups'] for q in records),
             'max_worker_rss':max(r['max_rss_bytes'] for _,r in results),
             'completed_workers':[r for _,r in results],'records':records}
    (args.output/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    print('solved',len(records),'prefixes, stats',summary['independent_sets_processed'],summary['max_worker_rss'])


if __name__=='__main__':main()
