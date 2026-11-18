"""Check local Metal output against a true integer normalizer and histogram."""
import hashlib
import json
import math
from pathlib import Path

P=Path(__file__).resolve().parent
C=P.parents[1]

def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def main():
    receipt=json.loads((P/'receipt.json').read_text())
    raw=P/'n7.jsonl'
    resource=P/'resource.txt'
    src=C/'r026/r026-merge-rb-calibration/source/tau2_merge_rb_calibration_r026.swift'
    assert sha(src)==receipt['source_sha256']
    assert sha(raw)==receipt['raw_sha256'] and sha(resource)==receipt['resource_sha256']
    assert receipt['returncode']==receipt['caffeinate_returncode']==0
    assert receipt['waited_reaped'] and not receipt['process_presence_after_wait']
    rows=[json.loads(x) for x in raw.read_text().splitlines()]
    m=rows[0]
    assert m['implementation']=='r026-tau2-firstmerge-rb-terminal'
    assert m['tau_lo']==m['tau_hi']==2 and m['method']=='strat'
    n,p,r,l,domain=(m[k] for k in ('n','particles','replicates','move_steps','domain'))
    assert (n,p,r,l,domain)==(7,65536,16,128,'9e50128000000001')
    assert receipt['parameters']==dict(n=n,P=p,R=r,tau=2,L=l,domain=domain)
    assert rows[-1]['kind']=='run_end' and len(rows)==2*r+2==receipt['rows']
    groups=rows[1:-1:2]
    lift=rows[2:-1:2]
    scale=max(x['logZ'] for x in groups)
    z=[math.exp(x['logZ']-scale) for x in groups]
    for j,(x,y) in enumerate(zip(groups,lift,strict=True)):
        assert x['kind']=='replicate' and y['kind']=='firstmerge_rb_calibration'
        assert x['replicate']==y['replicate']==j and x['tau']==y['tau']==2
        assert not x['dead'] and x['Z']>0 and math.isfinite(x['logZ'])
        seed=(0xA022000000000000+int(domain,16)+n*1_000_000+20000+j)%(1<<64)
        assert seed==int(x['seed_hex'],16)
        assert len(y['mean_m_hist'])==2*n-1 and y['mean_m_hist'][:2]==[0,0]
        assert math.isclose(sum(y['mean_m_hist']),y['mean_m'],rel_tol=1e-12)
        assert math.isclose(y['mean_m_F0']-2*y['mean_m_G0'],4*y['mean_m_d2'],rel_tol=1e-12)
    mass=[a*b['mean_m'] for a,b in zip(z,lift,strict=True)]
    total=sum(mass)
    estimate_A2=math.exp(scale)*total/r
    n7=C/'r026/r026-merge-rb-calibration/exact-top2-d4/n7.json'
    ex=json.loads(n7.read_text())
    fact=math.factorial(n*n)
    A=int(ex['A2'])
    hist=[sum(a*b['mean_m_hist'][d] for a,b in zip(z,lift,strict=True))/total for d in range(2*n-1)]
    q=sum(a*b['mean_m_d2'] for a,b in zip(z,lift,strict=True))/(total*n*n)
    far=sum(hist[n//2+1:])
    exactq=int(ex['squared_euclidean_total'])/(A*n*n)
    exactfar=sum(int(v) for d,v in ex['distance_counts'].items() if int(d)>n/2)/A
    answer=dict(input_sha256=sha(raw),exact_sha256=sha(n7),parameters=receipt['parameters'],
                wall_seconds=rows[-1]['wall_seconds'],reported_gpu_seconds=sum(x['gpu_seconds'] for x in groups),
                ess=total**2/sum(v*v for v in mass),max_replicate_weight=max(mass)/total,
                ratio_A2_estimated_to_exact=estimate_A2/(A/fact),Q_estimate=q,Q_exact=exactq,
                far_half_estimate=far,far_half_exact=exactfar,
                diagnostic_only=True)
    print(json.dumps(answer,indent=2))

if __name__=='__main__':
    main()
