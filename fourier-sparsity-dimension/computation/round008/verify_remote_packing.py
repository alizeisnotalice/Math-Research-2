#!/usr/bin/env python3
"""Independently reconstruct every final remote affine-packing certificate."""
import argparse
import hashlib
import json
from pathlib import Path
import sys


def rank(rows):
    p={}
    for x in rows:
        while x:
            b=x.bit_length()-1
            if b in p:x^=p[b]
            else:p[b]=x;break
    return len(p)


def verify(d):
    assert d['status']=='passed' and d['n'] in (8,9) and d['codimension']==3
    n=d['n'];N=1<<n;used=0;normals=[]
    assert 0<len(d['flats'])<=8 and len(d['flats'])==d['best_flat_count']
    for f in d['flats']:
        A=f['normal_rows'];b=f['rhs'];assert len(A)==len(b)==3
        assert all(0<x<N for x in A) and all(x in (0,1) for x in b) and rank(A)==3
        words=f['point_mask_words_hex'];assert len(words)==N//64
        assert all(0<=int(x,16)<1<<64 for x in words)
        recorded=sum(int(x,16)<<(64*j) for j,x in enumerate(words))
        direct=sum(1<<x for x in range(N) if all(((row&x).bit_count()&1)==rhs for row,rhs in zip(A,b)))
        assert recorded==direct and direct.bit_count()==N//8 and not(direct&used)
        used|=direct;normals+=A
    actual=rank(normals)
    assert actual==d['best_normal_rank'] and d['target_normal_rank']==8
    assert bool(d['counterexample_found'])==(actual>=8)
    assert 0<d['iterations_executed']<=d['iterations_requested']
    if not d['counterexample_found']:assert d['iterations_executed']==d['iterations_requested']
    return {'n':n,'seed':d['seed'],'iterations':d['iterations_executed'],
            'flats':len(d['flats']),'normal_rank':actual,'counterexample_found':d['counterexample_found']}


def main():
    ap=argparse.ArgumentParser();ap.add_argument('receipt',type=Path);ap.add_argument('--out',type=Path,required=True);a=ap.parse_args()
    raw=a.receipt.read_bytes();receipt=json.loads(raw);r=receipt['result'];here=Path(__file__).resolve().parent
    assert receipt['ssh_exit_code']==0 and r['status']=='passed' and r['workers']==14
    assert receipt['runner_sha256']==hashlib.sha256((here/'run_remote_packing.py').read_bytes()).hexdigest()
    assert r['source_sha256']==hashlib.sha256((here/'packing_search.cpp').read_bytes()).hexdigest()
    assert len(r['jobs'])==14 and len({j['result']['seed'] for j in r['jobs']})==14
    rows=[]
    for j in r['jobs']:
        assert j['exit_code']==0
        d=j['result'];cmd=j['command']
        for flag,key in [('--n','n'),('--seed','seed'),('--iterations','iterations_requested')]:
            assert int(cmd[cmd.index(flag)+1])==d[key]
        rows.append(verify(d))
    pilots=[verify(p['result']) for p in r['pilot_records']]
    result={'status':'passed','scope':'finite heuristic searches; exact final witnesses, no exhaustive rank bound',
            'arithmetic':'integer/F2 certificate verification; floating annealing used only for proposal acceptance',
            'jobs':rows,'pilots':pilots,'total_search_iterations':sum(x['iterations'] for x in rows),
            'rank8_hits':sum(x['counterexample_found'] for x in rows),
            'peak_task_cpu_percent':max(x['percent_sum'] for x in r['cpu_samples']),
            'whole_machine_top':r['whole_machine_top'],'all_workers_exited_zero':True,
            'remote_elapsed_seconds':r['elapsed_seconds'],'python':sys.version,
            'input_sha256':hashlib.sha256(raw).hexdigest(),
            'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    a.out.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k not in ['jobs','pilots','python']}))


if __name__=='__main__':main()
