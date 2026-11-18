"""Verify raw first-contact rank counts and audit replicate weight concentration."""
from hashlib import sha256
from math import exp,isclose
from pathlib import Path
import argparse
import json


def rows(path):
    return [json.loads(s) for s in path.read_text().splitlines()]


def series(payload):
    m=payload[0]
    n,L,P,R=m['n'],m['move_steps'],m['particles'],m['replicates']
    assert len(payload)==2*R+2 and payload[-1]['kind']=='run_end'
    spectra=[z['tau2_distance_scaled_mean_L_hist'] for z in payload[1:-1:2]]
    logs=[z['logZ'] for z in payload[1:-1:2]]
    weights=[exp(x-max(logs)) for x in logs]
    masses=[w*sum(y) for w,y in zip(weights,spectra,strict=True)]
    total=sum(masses)
    return {'n':n,'L':L,'P':P,'R':R,'replicate_ess':total**2/sum(x*x for x in masses),
            'max_replicate_weight':max(masses)/total,
            'near_distance_lt_13':sum(w*sum(y[:13]) for w,y in zip(weights,spectra,strict=True))/total},weights


def verify(base,old,source):
    f=base/'run.jsonl'
    meta=json.loads((base/'finished.json').read_text())
    assert meta['returncode']==meta['caffeinate_returncode']==0 and meta['waited_reaped']
    assert sha256(f.read_bytes()).hexdigest()==meta['raw_sha256']
    assert sha256(source.read_bytes()).hexdigest()==meta['source_sha256']
    assert sha256((base/'resource.txt').read_bytes()).hexdigest()==meta['resource_sha256']
    raw=rows(f)
    new,w=series(raw)
    assert (new['n'],new['L'],new['P'],new['R'])==(24,16,65536,8)
    old_summary,_=series(rows(old))
    assert (old_summary['n'],old_summary['L'],old_summary['P'],old_summary['R'])==(24,64,65536,8)
    old_summary['same_seed_domain']=raw[0]['domain']==rows(old)[0]['domain']
    size=47
    hist=[0.0]*576
    near_count=reach=overall=0.0
    near_particles=0
    for i,(p,q,wi) in enumerate(zip(raw[1:-1:2],raw[2:-1:2],w,strict=True)):
        expected=(0xA022000000000000+int(raw[0]['domain'],16)+24_000_000+20000+i)%(1<<64)
        assert int(p['seed_hex'],16)==expected
        rec=q['contact_counts']
        ranks=rec['contact_rank_L']
        matrix=rec['contact_joint_source_max_distance_L']
        assert len(ranks)==576 and len(matrix)==size*size
        old_hist=[sum(int(g[d]) for g in p['tau2_distance_L_hist_groups']) for d in range(size)]
        assert [sum(matrix[d*size:(d+1)*size]) for d in range(size)]==old_hist
        assert sum(matrix)==sum(ranks)==sum(map(int,p['tau2_lift_L_groups']))
        assert rec['first_contact_nonmerge_L'][0]<=sum(ranks)
        overall+=wi*sum(ranks)
        near=sum(matrix[:13*size])
        near_particles+=near>0
        near_count+=wi*near
        reach+=wi*sum(matrix[d*size+far] for d in range(13) for far in range(13,size))
        for j,value in enumerate(ranks):
            hist[j]+=wi*value
    assert isclose(near_count/overall,new['near_distance_lt_13'],rel_tol=1e-12)
    new.update(gpu_seconds=sum(r['gpu_seconds'] for r in raw[1:-1:2]),
               elapsed_seconds=raw[-1]['wall_seconds'],
               diagnostic_mean_first_contact=sum(j*x for j,x in enumerate(hist))/overall,
               diagnostic_rank_le_8=sum(hist[:9])/overall,
               nonzero_near_replicates=near_particles,
               accessible_among_near_diagnostic=reach/near_count)
    return {'n24_new':new,'n24_r053_different_domain':old_summary,
            'scope':'Finite numerical diagnostics with substantial replicate concentration; no mixing error bound'}


if __name__=='__main__':
    p=argparse.ArgumentParser()
    p.add_argument('batch',type=Path)
    p.add_argument('older_r53',type=Path)
    p.add_argument('frozen_source',type=Path)
    a=p.parse_args()
    print(json.dumps(verify(a.batch,a.older_r53,a.frozen_source),indent=2))
