#!/usr/bin/env python3
"""Negative and log-underflow checks for the strict r021 aggregator."""
import argparse, json, pathlib, shutil, subprocess, sys, tempfile

def invoke(aggregator, path, ok):
    cp=subprocess.run([sys.executable,str(aggregator),str(path),'--n','3'],text=True,capture_output=True)
    if bool(cp.returncode==0)!=ok:
        raise AssertionError((path.name,ok,cp.returncode,cp.stdout[-300:],cp.stderr[-300:]))
    return cp
def main():
    p=argparse.ArgumentParser();p.add_argument('jsonl');p.add_argument('--aggregator');a=p.parse_args();src=pathlib.Path(a.jsonl);root=src.parent
    aggregator=pathlib.Path(a.aggregator) if a.aggregator else root/'aggregate_random_tau.py'
    rows=[json.loads(s) for s in src.read_text().splitlines()]
    with tempfile.TemporaryDirectory(dir=root) as td:
        td=pathlib.Path(td)
        def materialize(name, mutate=None):
            rr=[dict(x) for x in rows]
            if mutate: mutate(rr)
            q=td/(name+'.jsonl');q.write_text(''.join(json.dumps(x)+'\n' for x in rr))
            shutil.copy(src.with_suffix(src.suffix+'.manifest.json'),q.with_suffix(q.suffix+'.manifest.json'))
            return q
        invoke(aggregator,materialize('valid'),True)
        invoke(aggregator,materialize('wrong_q',lambda r:r[0].update(qnum=r[0]['qnum']+1)),False)
        def seed(r):r[1]['desired_seed_hex']=r[0]['desired_seed_hex'];r[1]['seed_hex']=r[0]['seed_hex']
        invoke(aggregator,materialize('duplicate_seed',seed),False)
        invoke(aggregator,materialize('missing',lambda r:r.pop()),False)
        def dead(r):
            x=next(x for x in r if not x['dead']);x['dead']=True
        invoke(aggregator,materialize('bad_dead',dead),False)
        def underflow(r):
            x=next(x for x in r if not x['dead']);x.update(logZ='-1000',Z=0,U=0,D=0,F=0,G=0,F_minus_2G_minus_4D=0)
        invoke(aggregator,materialize('underflow',underflow),True)
    print('aggregate_guard_tests=PASS: valid, q, duplicate-seed, missing-row, dead-payload, logZ=-1000')
if __name__=='__main__':main()
