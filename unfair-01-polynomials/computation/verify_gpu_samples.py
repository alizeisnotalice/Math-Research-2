#!/usr/bin/env python3
"""Exact CPU recomputation of saved floating-point GPU samples.

This verifies the actual binary-rational sample, not a nearby algebraic point.
Failure of a sample never excludes its support or a degree range.
"""
from fractions import Fraction as Q
import hashlib
import json
from pathlib import Path
import sys


def check(path):
    report=json.loads(Path(path).read_text())
    base=Path(__file__).resolve().parent
    checked=[]
    for run in report['runs']:
        c=run['config']
        assert hashlib.sha256((base/run['source']).read_bytes()).hexdigest()==c['script_sha256']
        sample=run['best_saved_case']
        A=[Q.from_float(float(x)) for x in sample['A']]
        B=[Q.from_float(float(x)) for x in sample['B']]
        assert len(A)==c['m']+1 and len(B)==c['n']+1
        assert A[0]==A[-1]==B[0]==B[-1]==1
        assert all(0<=v<=1 for v in A+B)
        d=c.get('fractional_index',1)
        assert all(A[i]==B[i]==0 for i in range(1,d))
        assert 0<A[d]<1 and 0<B[d]<1
        C=[sum((A[i]*B[k-i] for i in range(len(A)) if 0<=k-i<len(B)),Q(0))
           for k in range(len(A)+len(B)-1)]
        errors=[min(abs(v),abs(v-1)) for v in C]
        error=max(errors)
        assert abs(float(error)-sample['max_residual_float64'])<1e-12
        assert error>0, 'Exact binary product found: escalate for full validation.'
        checked.append({'degree_pair':[c['m'],c['n']],'fractional_index':d,
                        'exact_max_distance_to_binary':str(error),
                        'max_distance_decimal':float(error),
                        'sample_is_exact_counterexample':False})
    return {'verified_saved_samples':len(checked),'arithmetic':'Fraction, exact binary-rational inputs',
            'results':checked,'scope':'Only the saved numerical samples; no finite exclusion follows.'}


if __name__=='__main__':print(json.dumps(check(sys.argv[1]),indent=2))
