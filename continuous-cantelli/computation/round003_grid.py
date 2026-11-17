"""Exact finite-time-grid audit for heat-delayed stopping measures.
Not a search for an original Cantelli counterexample.
"""
from fractions import Fraction as F
from collections import defaultdict
from pathlib import Path
import json, sys, platform, hashlib, os
s,v=F(1,2),F(1,8)
u=1-s+v
weights=[F(1,3),F(11,20),F(2,3)]
def merge(terms):
    d=defaultdict(F)
    for variance,weight in terms:d[variance]+=weight
    return {k:z for k,z in d.items() if z}
rows=[]
for w in weights:
    assert w*w/u<1 # positivity of late measure p_1-w*p_u
    exact_threshold=(w*w*s-v)/(1-w*w)
    for n in range(1,65):
        early_index=(n+1)//2
        rounded=F(early_index,n)
        delay=rounded-s
        assert delay>=0 and delay<=F(1,n)
        # Source early component evolves before its new stopping time.
        early_variance=v+delay
        # End-point distribution is identically p_1, not merely moments.
        terminal=merge([(early_variance+1-rounded,w),(F(1),F(1)),(u,-w)])
        assert terminal=={F(1):F(1)}
        if rounded==1:
            # All source layers coincide; sum to p_1.
            passed=True; squared_ratio=None
        else:
            # Q_rounded=p_rounded-w*p_(v+delay).
            squared_ratio=w*w*rounded/early_variance
            passed=squared_ratio<=1
            assert passed==(delay>=exact_threshold)
        rows.append({'w':str(w),'n':n,'rounded_early_time':str(rounded),
                     'delay':str(delay),'early_ratio_squared':None if squared_ratio is None else str(squared_ratio),
                     'grid_positive':passed})
    if w<=F(1,2):assert all(z['grid_positive'] for z in rows if z['w']==str(w))
# The invalid w=11/20 source is hidden on several coarse odd meshes.
passes=[z['n'] for z in rows if z['w']=='11/20' and z['grid_positive']]
assert passes==[1,3,5,7,9,11,13]
assert (F(11,20)**2*s-v)/(1-F(11,20)**2)==F(7,186)
# Nested dyadic grids detect any such w>1/2 at the first mesh containing s.
for z in rows:
    if z['w'] in ['11/20','2/3'] and z['n'] in [2,4,8,16,32,64]:assert not z['grid_positive']
summary={str(w):{'passes':[z['n'] for z in rows if z['w']==str(w) and z['grid_positive']],
                'minimum_delay_for_separated_layers':str((w*w*s-v)/(1-w*w)),
                'source_has_stopping_realization':w<=F(1,2)} for w in weights}
result={'parameters':{'s':str(s),'v':str(v),'u':str(u),'n_range':[1,64]},
        'summary':summary,'rows':rows,
        'runtime':{'python':platform.python_version(),'pid':os.getpid(),'executable':sys.executable},
        'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        'arithmetic':'exact Fraction; deterministic; no tolerances or random seed',
        'scope':'general two-time sources, not Gaussian-spatial continuous-graph Cantelli candidates'}
Path(sys.argv[1]).write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps({'summary':summary,'row_count':len(rows),'runtime':result['runtime'],'script_sha256':result['script_sha256']},indent=2))
