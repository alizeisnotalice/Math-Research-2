"""Exact complete bbox census for a two-row maximum-block family.

All coordinates use 1/(32N); masses use 1/(4N). No sampled geometry.
The continuous N-family theorem is separate from these finite checks.
"""
from concurrent.futures import ProcessPoolExecutor
from fractions import Fraction as F
from pathlib import Path
import hashlib,json,os,platform,subprocess,time

def one(N):
    t=time.perf_counter();best=0;max_count=0;total=0;eligible=0
    def add(units,dx,dy,multiplicity=1):
        nonlocal best,max_count,total,eligible
        total+=multiplicity
        if 256*N*units>max(dx,dy)**2:
            eligible+=multiplicity
            if units>best:best=units;max_count=multiplicity
            elif units==best:max_count+=multiplicity
    # Non-external closed boxes: contiguous x interval, one or both rows.
    for i in range(N):
        for j in range(i,N):
            length=j-i+1;dx=28*(j-i)
            add(2*length,dx,0,2)
            add(4*length,dx,28*N)
    # External boxes: the external source plus a suffix of one or both rows.
    for k in range(N):
        length=N-k;dx=39*N-28*k
        add(N+2*length,dx,14*N,2)
        add(N+4*length,dx,28*N)
    add(N,0,0)
    assert total==3*N*(N+1)//2+3*N+1
    assert best==4*N and max_count==1
    ks=range(N//4+1,15*N//32)
    capacities=[]
    for k in ks:
        units=N+4*(N-k);dx=39*N-28*k
        assert 3*N<units<4*N
        assert 256*N*units>max(dx,28*N)**2
        a=F(N-k,N);assert a>F(17,32);capacities.append(a)
    assert len(capacities)==7*N//32-1
    assert sum(capacities)>F(17,32)*(7*N//32-1)
    assert 16*(28*(N-1))**2>9*(32*N)**2
    assert 16*(28*N)**2>9*(32*N)**2
    assert F(39,32)**2>F(35,33)*F(5,4)
    return {'N':N,'sources':2*N+1,'pid':os.getpid(),'all_closed_bboxes':total,
            'eligible_closed_bboxes':eligible,'maximum_mass':'1','maximizers':1,
            'S_widths':[str(F(7*(N-1),8*N)),'7/8'],'both_axes_strictly_fat':True,
            'selected_crossing_high_boxes':len(capacities),
            'selected_removed_capacity_sum':str(sum(capacities)),
            'capacity_sum_lower_formula':str(F(17,32)*(7*N//32-1)),
            'connected_reason':'The full S is eligible and every selected crossing suffix connects the outside point to S.',
            'nonterminal_whole':True,'seconds':time.perf_counter()-t}

if __name__=='__main__':
    assert os.path.ismount('/Volumes/data');start=time.monotonic()
    with ProcessPoolExecutor(max_workers=12) as pool:rows=list(pool.map(one,[32*i for i in range(1,13)]))
    pids=sorted({r['pid'] for r in rows});absent={str(pid):not subprocess.run(['ps','-p',str(pid),'-o','pid='],capture_output=True,text=True).stdout.strip() for pid in pids}
    out={'status':'PASS_EXACT_FINITE_FAMILY_CHECKS','script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
         'python':platform.python_version(),'requested_workers':12,'actual_worker_pids':pids,'owned_workers_absent':absent,
         'wall_seconds':time.monotonic()-start,'cases':rows,'total_closed_bboxes':sum(r['all_closed_bboxes'] for r in rows),
         'scope':'Complete exact bbox census for 12 explicit N. This confirms repeated per-witness capacity even in a connected nonterminal unique-fat-max family; no failure of union budget or deletion selection is claimed.'}
    Path(__file__).with_suffix('.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({k:v for k,v in out.items() if k!='cases'},indent=2))
