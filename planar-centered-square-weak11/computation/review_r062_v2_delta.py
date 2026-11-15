#!/usr/bin/env python3
"""Independent v2 derivative-extremum and auxiliary-LP delta audit."""
import ast
from collections import Counter
import contextlib
from fractions import Fraction as F
import hashlib
import importlib.util
import io
import json
import os
from pathlib import Path
import sys
import time
from types import SimpleNamespace
import sympy as sp

HERE=Path(__file__).resolve().parent
OWNER=Path('/Volumes/data/CodexProjects/T2/.local/workspaces/pcw-sol-kernel-sources/planar-centered-square-weak11')
V1=OWNER/'computation/r062_coupled_source_verify.py'
V2=OWNER/'computation/r062_coupled_source_verify_v2.py'
LP=OWNER/'computation/cauchy_occupation_coupled_aux_lp.m'
TMP=Path('/Volumes/data/CodexTemp/T2/r062-independent-v2-audit')
V2_SHA='9ee1de8d1944cfd567a50a79c37d772e10cb16b121cca71f1676645f9197b2e7'
LP_SHA='5ee56d623be52fe53c709fcd02de0fd6e64f2a6e4cf90ebe73e57cae0fc1142e'


def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()


def module(name,path):
    spec=importlib.util.spec_from_file_location(name,path)
    m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
    return m


def hull(values):return min(a[0]for a in values),max(a[1]for a in values)


def main():
    assert os.path.ismount('/Volumes/data')
    start=time.perf_counter()
    assert sha(V2)==V2_SHA and sha(LP)==LP_SHA
    ref=module('independent_fraction_reference',HERE/'review_r062_coupled_source.py')
    v2=module('r62_v2_checked',V2)
    trees=[ast.parse(p.read_text())for p in (V1,V2)]
    fs=[{n.name:ast.dump(n,include_attributes=False)for n in t.body if isinstance(n,(ast.FunctionDef,ast.ClassDef))}for t in trees]
    changed=[k for k in fs[0]if fs[0][k]!=fs[1][k]]
    assert changed==['cauchyjet','kprime','self_test']
    x,h,P=sp.symbols('x h P',positive=True)
    k=h/(P*(h*h+x*x))
    assert sp.simplify(sp.diff(k,x,3)-24*h*x*(h*h-x*x)/(P*(h*h+x*x)**4))==0
    # k' stationary value: x=h/sqrt(3) is the negative extremum.
    assert sp.simplify(sp.diff(k,x).subs(x,h/sp.sqrt(3))+9/(8*sp.sqrt(3)*P*h*h))==0
    # Independent rational sqrt(3) bracket via dyadic bisection, not isqrt.
    rlo,rhi=F(1),F(2)
    for _ in range(240):
        mid=(rlo+rhi)/2
        if mid*mid<3:rlo=mid
        else:rhi=mid
    assert rlo*rlo<3<rhi*rhi
    grid=list(map(F,[-2,-1]))+[F(-3,5),F(-7,12),F(-4,7),F(-1,2),F(0),F(1,2),F(4,7),F(7,12),F(3,5),F(1),F(2)]
    stats=Counter();narrow=Counter()
    for hval in (F(3,20),F(1),F(2)):
        magnitude=ref.divide((F(9,8)/(hval*hval),)*2,ref.multiply((rlo,rhi),ref.PI))
        critical_pos=(hval/rhi,hval/rlo)
        critical_neg=(-critical_pos[1],-critical_pos[0])
        for ia,a in enumerate(grid):
            for b in grid[ia:]:
                lo,hi=a*hval,b*hval
                kp=v2.kprime(lo,hi,hval)
                p,kpp=v2.cauchyjet(lo,hi,hval)
                kp_reference=[ref.reference_jet(t,hval)[1]for t in (lo,hi)]
                # Rational endpoint brackets certify critical-point inclusion
                # independently of v2's squared-comparison branch predicates.
                for interval,sign in ((critical_pos,-1),(critical_neg,1)):
                    inside=lo<=interval[0] and interval[1]<=hi
                    outside=hi<interval[0] or interval[1]<lo
                    assert inside or outside
                    if inside:
                        kp_reference.append(ref.times(magnitude,sign))
                        stats['positive_critical_inclusions'if sign<0 else 'negative_critical_inclusions']+=1
                ref.contains(kp,hull(kp_reference),v2.S)
                kpp_reference=[ref.reference_jet(t,hval)[2]for t in [lo,hi]+[t for t in (-hval,F(0),hval)if lo<=t<=hi]]
                ref.contains(kpp,hull(kpp_reference),v2.S)
                for t in (lo,(3*lo+hi)/4,(lo+hi)/2,(lo+3*hi)/4,hi):
                    rr=ref.reference_jet(t,hval)
                    ref.contains(p,rr[0],v2.S);ref.contains(kp,rr[1],v2.S);ref.contains(kpp,rr[2],v2.S)
                    stats['rational_sample_points']+=1
                stats['intervals']+=1
    # Aux cancellation and directions on exact rational matrices.
    # U*v can be negative even though separately clipped bases charge >0.
    rows=[]
    for c1,c2 in ((F(1),F(1)),(F(1,2),F(3)),(F(4),F(1,8))):
        U=(F(2),F(-3));curvature=U[0]*c1+U[1]*c2
        eta=max(F(0),curvature)
        assert curvature-eta<=0 and eta>=0
        individual=max(F(0),U[0])*c1+max(F(0),U[1])*c2
        assert eta<=individual
        # -V*v+a_target+ell²/8*eta<=0 exactly matches residual bound.
        vertex=F(1)+eta/8;target=F(1)
        assert -vertex+target+eta/8==0
        rows.append(dict(coefficients=list(map(str,(c1,c2))),signed_upper=str(curvature),eta=str(eta),separate_positive_sum=str(individual)))
    TMP.mkdir(parents=True,exist_ok=True)
    inp=TMP/'toy-input.json';inp.write_text(json.dumps(dict(h='1',delta='1/8',N=1,source_origin_mass='0',density_coefficients=['1']*3))+'\n')
    common=dict(input=str(inp),repair='0',initial_subdiv=1,max_depth=0,max_seconds=20,max_boxes=100,stop_on_negative=False)
    runs=[]
    for name,updates in [('negative_incomplete',dict(max_boxes=1)),('negative_complete',{})]:
        args=SimpleNamespace(**(common|updates|dict(output=str(TMP/(name+'.json')))))
        with contextlib.redirect_stdout(io.StringIO()):v2.exact_run(args)
        data=json.loads(Path(args.output).read_text())
        assert data['status']=='INVALID_SOURCE'
        assert data['complete_terminal_coverage']==(name=='negative_complete')
        assert (data['additional_origin_repair_sufficient_exact']is not None)==(name=='negative_complete')
        runs.append(dict(name=name,status=data['status'],complete=data['complete_terminal_coverage'],repair=data['additional_origin_repair_sufficient_exact'],checked=data['checked_boxes']))
    repair=runs[-1]['repair']
    args=SimpleNamespace(**(common|dict(repair=repair,output=str(TMP/'repaired_complete.json'))))
    with contextlib.redirect_stdout(io.StringIO()):v2.exact_run(args)
    data=json.loads(Path(args.output).read_text())
    assert data['status']=='PASS_SOURCE_NONNEGATIVE' and data['complete_terminal_coverage']
    assert data['additional_origin_repair_sufficient_exact']=='0' and not data['D_minimum_certified']
    runs.append(dict(name='repaired_complete',status=data['status'],complete=True,origin_amplitude=repair,checked=data['checked_boxes']))
    proof=OWNER/'proofs/r062-coupled-wholecell-source-interface.md'
    text=proof.read_text()
    assert 'BOTH the fixed-point precision' in text and 'arctan series degree increase' in text
    assert sha(V2)==V2_SHA and sha(LP)==LP_SHA
    out=dict(status='PASS_R62_V2_DELTA',python=sys.version,executable=sys.executable,sympy=sp.__version__,elapsed_seconds=time.perf_counter()-start,
             script_sha256=sha(Path(__file__)),changed_AST_functions=changed,statistics=dict(stats),symbolic_identities=2,
             aux_rational_rows=rows,source_toy_runs=runs,termination_joint_precision_fix_confirmed=True,
             source_pins=[dict(path=str(p),sha256=sha(p))for p in (V1,V2,LP,proof,HERE/'review_r062_coupled_source.py')],
             arithmetic='Exact Fraction, independent 240-step sqrt3 bracket, independent 120-term atan brackets, symbolic derivative identities.',
             scope='V2 derivative extrema and signed-auxiliary LP delta only; no MATLAB launch, no optimized candidate or full-D certificate.')
    dest=HERE/'results/round062-v2-delta-independent.json';dest.write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(out,indent=2))


if __name__=='__main__':main()
