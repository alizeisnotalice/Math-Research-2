#!/usr/bin/env python3
"""CPU polynomial reconstruction of saved remote samples, not a certificate."""
import json, hashlib, sys, time
from pathlib import Path
import numpy as np

P=Path(__file__).resolve().parent
t=json.loads((P/'root-task2.json').read_text())
g=t['groups']; assert len(t['roots'])==g

def pair(mask):
    assert type(mask) is int and 0<=mask<2**(g-1)
    sides=[np.ones(1,dtype=complex),np.ones(1,dtype=complex)]
    bits=[0]+[(mask>>j)&1 for j in range(g-2,-1,-1)]
    for side,roots in zip(bits,t['roots']):
        p=np.polynomial.polynomial.polyfromroots([complex(*r) for r in roots])
        assert np.all(np.isfinite(p)) and abs(p[0])>0
        p=p/p[0]; assert np.max(np.abs(p.imag))<1e-7
        sides[side]=np.polynomial.polynomial.polymul(sides[side],p)
    a,b=sides
    assert np.all(np.isfinite(a)) and np.all(np.isfinite(b))
    assert len(a)+len(b)-1==len(t['C'])
    score=max(0.,-a.real.min(),-b.real.min(),a.real.max()-1,b.real.max()-1,np.max(np.abs(a.imag)),np.max(np.abs(b.imag)))
    residual=float(np.max(np.abs(np.convolve(a,b)-np.array(t['C']))))
    return float(score),residual

beg=time.monotonic();six=json.loads((P/'remote-six.json').read_text());one=json.loads((P/'remote-1m.json').read_text())
assert six['masks']==[0,1,17,1234,23166975,23166976]
assert six['command_status']==one['command_status']=='completed'
assert six['finite']==[1]*6 and one['finite_count']==one['range_count']==1048576
assert len(set(one['top_indices']))==4 and all(0<=m<1048576 for m in one['top_indices'])
checks=[]
for source,masks,scores in [('six',six['masks'],six['scores']),('1m_saved_top4',one['top_indices'],one['top_scores'])]:
    assert len(masks)==len(scores)
    for mask,score in zip(masks,scores):
        expected,residual=pair(mask)
        if mask==0:
            assert score is None; continue
        assert isinstance(score,(int,float)) and np.isfinite(score)
        error=abs(score-expected); assert error<1e-3
        checks.append(dict(source=source,mask=mask,gpu_score=score,cpu_score=expected,error=error,convolution_residual_float64=residual))
out=dict(status='passed_saved_remote_scores_independent_cpu_reconstruction',checks=checks,max_score_error=max(x['error'] for x in checks),seconds=time.monotonic()-beg,python=sys.version,numpy=np.__version__,source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),scope='Nine non-control saved samples only; no proof of global GPU top4 correctness, full coverage, positivity or exact algebraic factorization.')
(P/'parent-native-score-check.json').write_text(json.dumps(out,indent=2,sort_keys=True)+'\n')
print(json.dumps({k:v for k,v in out.items() if k!='checks'},sort_keys=True))
