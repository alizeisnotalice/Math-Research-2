"""Exact split/domino comparison and small independent subset-DP gates."""
import argparse
from concurrent.futures import ThreadPoolExecutor, FIRST_COMPLETED, wait
from fractions import Fraction
import itertools
import json
from math import comb
from pathlib import Path
import re
import subprocess

def states(n, m, h):
    last = {0:1}
    for y in range(n):
        width = m if y%2 else m-1
        masks = [x for x in range(1<<width) if not x & (x<<1)]
        last = {a:sum(v for b,v in last.items() if not a & b) for a in masks}
    flex = sum(int(h[y])-m+1 for y in range(0,n,2))
    return sum(last.values()) * (1 << flex)

def profiles(n):
    m=n//2
    for choices in itertools.product((m-1,m+1),repeat=(n+1)//2):
        choices=iter(choices)
        yield ''.join(str(m if y%2 else next(choices)) for y in range(n))

def flip(h,n):
    return ''.join(str(n-int(w)) for w in h)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('output',type=Path)
    parser.add_argument('--workers',type=int,default=8)
    parser.add_argument('--counter',type=Path,required=True)
    parser.add_argument('--direct',type=Path,required=True)
    args=parser.parse_args()
    ns=(4,5,6,7)
    families={n:tuple(profiles(n)) for n in ns}
    jobs=set()
    for n,family in families.items():
        m=n//2
        for h in (*family,str(m)*n):
            jobs.update(((n,m,h),(n,n-m,flip(h,n))))
    todo=sorted(jobs,key=lambda z: -states(*z))
    records={}
    def run(job):
        n,m,h=job
        child=subprocess.Popen(['/usr/bin/time','-l',str(args.counter),str(n),str(m),h],
                               stdout=subprocess.PIPE,stderr=subprocess.PIPE,text=True)
        stdout,stderr=child.communicate()
        assert child.returncode==0,(job,stderr)
        count=int(stdout)
        paired=(sum(map(int,h))<=25 and n<=6)
        if paired:
            assert count==int(subprocess.check_output([args.direct,str(n),str(m),h]))
        rss=re.search(r'^\s*(\d+)\s+maximum resident set size$',stderr,re.M)
        assert rss
        return dict(count=count,peak_rss=int(rss[1]),counter_pid=child.pid,waited=True,
                    independent_subset_dp_verified=paired, theory_states=states(*job))
    with ThreadPoolExecutor(max_workers=args.workers) as pool:
        pending={}
        while todo or pending:
            budget=20_000_000-sum(states(*job) for job in pending.values())
            while len(pending)<args.workers:
                choice=next((i for i,j in enumerate(todo) if states(*j)<=budget),None)
                if choice is None: break
                job=todo.pop(choice)
                pending[pool.submit(run,job)]=job
                budget-=states(*job)
            assert pending
            done,_=wait(pending,return_when=FIRST_COMPLETED)
            for f in done:
                records[pending.pop(f)]=f.result()
    rows={}
    for n,family in families.items():
        m=n//2
        def event(h):
            return (comb(n*n,sum(map(int,h)))*records[n,m,h]['count']*
                    records[n,n-m,flip(h,n)]['count'])
        counts={h:event(h) for h in family}
        total=sum(counts.values()); base=event(str(m)*n)
        assert total>=base>0
        rows[n]=dict(events=counts,union=total,split=base,gain=str(Fraction(total,base)))
        print(n,len(counts),total,float(Fraction(total,base)),flush=True)
    result=dict(families=rows,records=[dict(n=k[0],mid=k[1],profile=k[2],**v) for k,v in sorted(records.items())],
                maximum_peak_rss=max(x['peak_rss'] for x in records.values()))
    with args.output.open('x') as stream:
        json.dump(result,stream,indent=2)
        stream.write('\n')

if __name__=='__main__':
    main()
