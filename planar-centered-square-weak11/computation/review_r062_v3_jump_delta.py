#!/usr/bin/env python3
"""Independent symbolic telescope, exact grid jumps and curvature sample audit."""
import ast
from collections import defaultdict
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
V2=OWNER/'computation/r062_coupled_source_verify_v2.py'
V3=OWNER/'computation/r062_coupled_source_verify_v3.py'
LP=OWNER/'computation/cauchy_occupation_jump_aux_lp.m'
EXPECTED=['58cc151a64c30920c25fbb3a24f79ef7651e5d762210532678cf12ac5a39edc0','c23c30dd82c85fd07023bd103a1b54af2eb8fde820fb463294e48ff6e5d62bc0']
TMP=Path('/Volumes/data/CodexTemp/T2/r062-independent-v3-audit')


def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()


def load(name,path):
    spec=importlib.util.spec_from_file_location(name,path);m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m);return m


def main():
    assert os.path.ismount('/Volumes/data')
    start=time.perf_counter();assert [sha(V3),sha(LP)]==EXPECTED
    ref=load('fraction_reference_v1',HERE/'review_r062_coupled_source.py')
    v3=load('v3_checked_source',V3)
    # Pure symbolic telescope with independent interface basis symbols.
    aa=sp.symbols('a0:3');kk=sp.symbols('K0:4')
    cells=sum(aa[i]*(kk[i]-kk[i+1])for i in range(3))
    edges=aa[0]*kk[0]+sum((aa[i]-aa[i-1])*kk[i]for i in range(1,3))-aa[2]*kk[3]
    assert sp.expand(cells-edges)==0
    assert sp.expand(edges.subs(dict.fromkeys(aa,sp.Integer(1)))-(kk[0]-kk[3]))==0
    # AST guard: arithmetic, scalar derivative intervals and grid orbits frozen.
    trees=[ast.parse(p.read_text())for p in (V2,V3)]
    ff=[{n.name:ast.dump(n,include_attributes=False)for n in t.body if isinstance(n,(ast.FunctionDef,ast.ClassDef))}for t in trees]
    changed=[k for k in ff[0]if ff[0][k]!=ff[1][k]]
    assert changed==['Model','exact_run','self_test']
    mm=[]
    for tree in trees:
        cl=next(n for n in tree.body if isinstance(n,ast.ClassDef)and n.name=='Model')
        mm.append({n.name:ast.dump(n,include_attributes=False)for n in cl.body if isinstance(n,ast.FunctionDef)})
    assert [k for k in mm[0]if mm[0][k]!=mm[1][k]]==['__init__','box_bound']
    cases=[('constant',1,[F(1,32)]*3),('nonmonotone',1,[F(1,64),F(3,64),F(2,64)]),
           ('zeros_and_edge_heavy',2,[F(0),F(2,64),F(0),F(1,64),F(4,64),F(1,64)])]
    case_rows=[];samples=0;boundary_checks=0
    for name,N,coeff in cases:
        delta=F(1,8);h=F(1);lam=F(1)
        data=dict(h=str(h),delta=str(delta),N=N,source_origin_mass=str(lam),density_coefficients=list(map(str,coeff)))
        model=v3.Model(data,F(0))
        # Independent cell-divergence assembly, not neighbor differences.
        grid={};jx=defaultdict(F);jy=defaultdict(F)
        for i in range(-N,N+1):
            for j in range(-N,N+1):
                outer,inner=max(abs(i),abs(j)),min(abs(i),abs(j))
                a=coeff[outer*(outer+1)//2+inner]
                grid[i,j]=a
                jx[(F(i)-F(1,2))*delta,j*delta]+=a
                jx[(F(i)+F(1,2))*delta,j*delta]-=a
                jy[(F(j)-F(1,2))*delta,i*delta]+=a
                jy[(F(j)+F(1,2))*delta,i*delta]-=a
        expected_x={k:v for k,v in jx.items()if v};expected_y={k:v for k,v in jy.items()if v}
        assert expected_x=={(e,c):q for e,c,q in model.xjumps}
        assert expected_y=={(e,c):q for e,c,q in model.yjumps}
        for j in range(-N,N+1):
            assert jx[-(F(N)+F(1,2))*delta,j*delta]==grid[-N,j]
            assert jx[(F(N)+F(1,2))*delta,j*delta]==-grid[N,j]
            boundary_checks+=2
        if name=='constant':assert len(model.xjumps)==len(model.yjumps)==2*(2*N+1)
        if name=='nonmonotone':
            assert any(q>0 for _,_,q in model.xjumps)and any(q<0 for _,_,q in model.xjumps)
            assert grid[0,0]<grid[1,0] and grid[1,0]>grid[1,1]
        # Independent full-cell derivative reference; not a jump sum.
        lower=[]
        for target,((cx,cy),_)in enumerate(model.oo):
            B=(cx-delta/2,cx+delta/2,cy-delta/2,cy+delta/2)
            lo,err,repair,neg,curvature=model.box_bound(target,B)
            assert lo>0 and not neg and repair==0
            lower.append(str(lo))
            for x in (B[0],cx,B[1]):
                for y in (B[2],cy,B[3]):
                    rx=ref.multiply(ref.reference_jet(x,h)[2],ref.reference_jet(y,h)[0])
                    ry=ref.multiply(ref.reference_jet(x,h)[0],ref.reference_jet(y,h)[2])
                    for (i,j),a in grid.items():
                        if not a:continue
                        u,v=i*delta,j*delta
                        tx=ref.reference_prob(x-u,delta,h);ty=ref.reference_prob(y-v,delta,h)
                        ddx=ref.minus(ref.reference_jet(x-u+delta/2,h)[1],ref.reference_jet(x-u-delta/2,h)[1])
                        ddy=ref.minus(ref.reference_jet(y-v+delta/2,h)[1],ref.reference_jet(y-v-delta/2,h)[1])
                        rx=ref.plus(rx,ref.times(ref.multiply(ddx,ty),a));ry=ref.plus(ry,ref.times(ref.multiply(tx,ddy),a))
                    ref.contains(curvature[0],rx,v3.S);ref.contains(curvature[1],ry,v3.S)
                    samples+=1
        case_rows.append(dict(name=name,N=N,density_coefficients=list(map(str,coeff)),vertical_jumps=len(expected_x),horizontal_jumps=len(expected_y),positive_source_boxes=len(lower),lower_fixed=lower))
    # A signed interval product upper bound remains safe with nonminimal split.
    split_rows=[]
    for L,U,jump,extra in [(F(-3),F(2),F(-2),F(0)),(F(-3),F(2),F(2),F(5)),(F(1),F(4),F(-1),F(3))]:
        dp=max(jump,F(0))+extra;dn=max(-jump,F(0))+extra
        bound=U*dp-L*dn;true_upper=max(L*jump,U*jump)
        assert dp-dn==jump and bound>=true_upper
        assert bound-true_upper==(U-L)*extra
        split_rows.append(dict(L=str(L),U=str(U),jump=str(jump),extra=str(extra),bound=str(bound),true_upper=str(true_upper)))
    # Actual complete support/repair test with nonzero outer boundary jumps.
    TMP.mkdir(parents=True,exist_ok=True)
    inp=TMP/'toy-input.json';inp.write_text(json.dumps(dict(h='1',delta='1/8',N=1,source_origin_mass='0',density_coefficients=['1']*3))+'\n')
    common=dict(input=str(inp),repair='0',initial_subdiv=1,max_depth=0,max_seconds=20,max_boxes=100,stop_on_negative=False)
    runs=[]
    for name,updates in [('negative_incomplete',dict(max_boxes=1)),('negative_complete',{})]:
        args=SimpleNamespace(**(common|updates|dict(output=str(TMP/(name+'.json')))))
        with contextlib.redirect_stdout(io.StringIO()):v3.exact_run(args)
        out=json.loads(Path(args.output).read_text());assert out['status']=='INVALID_SOURCE'
        assert out['nonzero_vertical_jumps']==out['nonzero_horizontal_jumps']==6
        assert out['complete_terminal_coverage']==(name=='negative_complete')
        assert (out['additional_origin_repair_sufficient_exact']is not None)==(name=='negative_complete')
        runs.append(dict(name=name,status=out['status'],complete=out['complete_terminal_coverage'],repair=out['additional_origin_repair_sufficient_exact']))
    args=SimpleNamespace(**(common|dict(repair=runs[-1]['repair'],output=str(TMP/'repaired_complete.json'))))
    with contextlib.redirect_stdout(io.StringIO()):v3.exact_run(args)
    out=json.loads(Path(args.output).read_text());assert out['status']=='PASS_SOURCE_NONNEGATIVE'and out['complete_terminal_coverage']
    assert not out['D_minimum_certified']and out['additional_origin_repair_sufficient_exact']=='0'
    runs.append(dict(name='repaired_complete',status=out['status'],complete=True,origin_amplitude=args.repair))
    assert [sha(V3),sha(LP)]==EXPECTED
    out=dict(status='PASS_R62_V3_JUMP_DELTA',python=sys.version,executable=sys.executable,sympy=sp.__version__,elapsed_seconds=time.perf_counter()-start,
             script_sha256=sha(Path(__file__)),AST_changed_functions=changed,AST_changed_Model_methods=['__init__','box_bound'],
             symbolic_telescope_identities=2,cases=case_rows,curvature_fullcell_reference_points=samples,outer_boundary_checks=boundary_checks,
             signed_split_rows=split_rows,source_coverage_repair_runs=runs,
             source_pins=[dict(path=str(p),sha256=sha(p))for p in (V2,V3,LP,HERE/'review_r062_coupled_source.py')],
             scope='Jump identity, outer-zero boundaries, nonmonotone D4 coefficients, signed LP bound and exact checker delta only; no MATLAB launch/optimized candidate/full-D/new-C2 certificate.')
    dest=HERE/'results/round062-v3-jump-independent.json';dest.write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))


if __name__=='__main__':main()
