#!/usr/bin/env python3
"""Independent streaming audit of batch coverage, input binding, and partition structure."""
from pathlib import Path
import argparse,hashlib,json,sys
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
    ap=argparse.ArgumentParser()
    for name in ['tasks','results','summary','binary','out']:ap.add_argument('--'+name,type=Path,required=True)
    a=ap.parse_args();tasks={}
    for line in a.tasks.open():
        r=json.loads(line);assert r['id'] not in tasks;tasks[r['id']]=r
    seen=set();hist={}
    for line in a.results.open():
        r=json.loads(line);i=r['id'];assert i in tasks and i not in seen;seen.add(i)
        assert r['task']==tasks[i] and r['status']=='ok' and r['returncode']==0
        assert r['D']==tasks[i]['lift_dimension'] and r['D']+1<=r['kappa']<=tasks[i]['upper']
        z=json.loads(r['stdout']);assert z['D']==r['D'] and z['kappa']==r['kappa'] and z['partition']==r['partition']
        cells=[int(x,16) for x in z['partition']];union=0
        for c in cells:assert c and not(c&union);union|=c
        assert union==(1<<64)-1
        key=f"{r['D']}:{r['kappa']}";hist[key]=hist.get(key,0)+1
    assert seen==set(tasks)
    s=json.loads(a.summary.read_text())
    assert s['tasks_total']==s['records']==s['successful']==len(tasks)
    assert s['unstarted_tasks']==0 and not s['failures'] and not s['counterexamples']
    for key,f in [('tasks_sha256',a.tasks),('results_sha256',a.results),('binary_sha256',a.binary)]:assert s[key]==sha(f)
    out={'status':'passed','source_sha256':sha(Path(__file__)),'coverage':len(seen),'task_result_equality':True,'all_partitions_disjoint_and_complete':True,'all_trivial_lower_and_candidate_upper_checked':True,'histogram_D_kappa':hist,'python':sys.version,'tasks_sha256':sha(a.tasks),'results_sha256':sha(a.results),'summary_sha256':sha(a.summary),'binary_sha256':sha(a.binary),'scope':'All result/input bindings and structural checks; subcode ranks of all witnesses are not recomputed by this program.'}
    a.out.write_text(json.dumps(out,indent=2)+'\n');print(json.dumps({'status':'passed','coverage':len(seen)}))
if __name__=='__main__':main()
