#!/usr/bin/env python3
"""Final artifact/coverage binding audit; does NOT evaluate g or rerun checker."""
from collections import Counter
from fractions import Fraction as F
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
import time

HERE=Path(__file__).resolve().parent
ROOT=Path('/Volumes/data/CodexProjects/T2')
OWNER=ROOT/'.local/workspaces/pcw-sol-kernel-sources/planar-centered-square-weak11'
PROOFS={
 'proofs/r062-coupled-wholecell-source-interface.md':'ee97c8df889c5b7658f4405531d0db47d54eab48c3c88bd4bbb6a8e2be813730',
 'proofs/r062-jump-telescoped-source-curvature.md':'130a10cff83b6a15f4298b6a04a2c898335ad7eba1751c041d9f894a0d98522e',
}
CHAINS=[
 ('v1-h015-d125-s2','r062-pilot','occupation-coupled-lp.json','cauchy_occupation_coupled_lp.m','r062_coupled_source_verify.py','fb69861a37adaeba94ac57dfbf9921889970b2b53b07b2df166e50e39d29c5a7','65143813817936d4db538aa0a4b1f4f6de7d55e1070b211bf05f4211f9710827',612),
 ('v2-h015-d125-s2','r062-v2-pilot','occupation-coupled-aux-lp.json','cauchy_occupation_coupled_aux_lp.m','r062_coupled_source_verify_v2.py','5ee56d623be52fe53c709fcd02de0fd6e64f2a6e4cf90ebe73e57cae0fc1142e','9ee1de8d1944cfd567a50a79c37d772e10cb16b121cca71f1676645f9197b2e7',612),
 ('v3-h015-d125-s2','r062-v3-pilot','occupation-jump-aux-lp.json','cauchy_occupation_jump_aux_lp.m','r062_coupled_source_verify_v3.py','c23c30dd82c85fd07023bd103a1b54af2eb8fde820fb463294e48ff6e5d62bc0','58cc151a64c30920c25fbb3a24f79ef7651e5d762210532678cf12ac5a39edc0',612),
 ('v3-h015-d125-s4','r062-v3-h015-s4','occupation-jump-aux-lp.json','cauchy_occupation_jump_aux_lp.m','r062_coupled_source_verify_v3.py','c23c30dd82c85fd07023bd103a1b54af2eb8fde820fb463294e48ff6e5d62bc0','58cc151a64c30920c25fbb3a24f79ef7651e5d762210532678cf12ac5a39edc0',2448),
 ('v3-h0075-d125-s4','r062-v3-h0075-s4','occupation-jump-aux-lp.json','cauchy_occupation_jump_aux_lp.m','r062_coupled_source_verify_v3.py','c23c30dd82c85fd07023bd103a1b54af2eb8fde820fb463294e48ff6e5d62bc0','58cc151a64c30920c25fbb3a24f79ef7651e5d762210532678cf12ac5a39edc0',2448),
]


def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def pin(p):return dict(path=str(p),sha256=sha(p),bytes=p.stat().st_size)
def read(p):return json.loads(p.read_text(),parse_float=str)
def ceil(q):return -((-q.numerator)//q.denominator)


def main():
    assert os.path.ismount('/Volumes/data')
    start=time.perf_counter()
    for rel,digest in PROOFS.items():assert sha(OWNER/rel)==digest
    rows=[]
    for label,run,cname,mname,sname,msha,ssha,count in CHAINS:
        base=ROOT/'.local/compute'/run;archive=OWNER/'computation/results/r062-coupled-source'/label
        candidate,cert,receipt=base/cname,base/'source-exact.json',base/'receipt.json'
        c,s,r=read(candidate),read(cert),read(receipt)
        assert sha(base/mname)==sha(OWNER/'computation'/mname)==msha==r['source_sha256']
        assert sha(base/sname)==sha(OWNER/'computation'/sname)==ssha==s['script_sha256']
        assert sha(candidate)==r['output_sha256']==s['input_sha256']
        assert Path(s['input'])==candidate
        assert r['exit_status']==0 and r['owned_pid_absent'] is True
        assert r['command'][1:3]==['-singleCompThread','-batch']
        assert r['command'][-1].startswith(mname[:-2]+'(')
        current=subprocess.run(['ps','-p',str(r['pid']),'-o','pid=,comm='],capture_output=True,text=True)
        assert not current.stdout.strip()
        for name in (cname,'source-exact.json','receipt.json'):
            assert (base/name).read_bytes()==(archive/name).read_bytes()
        assert s['status']=='PASS_SOURCE_NONNEGATIVE' and s['complete_terminal_coverage'] is True
        assert s['source_only'] is True and s['D_minimum_certified'] is False
        assert F(s['origin_repair_added_exact'])==F(s['additional_origin_repair_sufficient_exact'])==0
        assert s['checked_boxes']==s['accepted_leaves']==count==len(s['accepted'])
        assert s['unresolved_leaves']==s['unchecked_stack']==0
        assert not any(s[k]for k in ('unresolved','unchecked','split_records','negative_corner_witnesses'))
        assert s['bits']==160 and s['atan_terms']==88
        N=c['N'];delta=F(c['delta']);h=F(c['h']);sub=c['subdiv'];lam=F(c['source_origin_mass'])
        assert N==16 and F(s['h'])==h and F(s['delta'])==delta and s['N']==N
        assert s['limits']['initial_subdiv']==sub and count==153*sub*sub
        a=list(map(F,c['density_coefficients']));assert len(a)==153 and min(a)>=0 and lam>=0
        assert F(s['source_origin_mass_exact'])==lam
        reps=[];sizes=[];expected=Counter()
        for i in range(N+1):
            for j in range(i+1):
                cx,cy=i*delta,j*delta;cell=len(reps);reps.append((cx,cy))
                physical={(sx*cx,sy*cy)for sx in (-1,1)for sy in (-1,1)}
                physical|={(y,x)for x,y in physical};sizes.append(len(physical))
                for ix in range(sub):
                    for iy in range(sub):
                        xl=cx-delta/2+delta*ix/sub;yl=cy-delta/2+delta*iy/sub
                        expected[cell,(xl,xl+delta/sub,yl,yl+delta/sub)]+=1
        assert reps==[tuple(map(F,p))for p in c['orbit_representatives']]
        assert sizes==c['orbit_sizes'] and sum(sizes)==c['physical_cells']==1089
        mass=lam+delta*delta*sum(q*n for q,n in zip(a,sizes))
        assert F(s['mu_mass_exact'])==mass and F(s['halfmass_over_pi_exact'])==h*mass
        actual=Counter();min_lower=None
        for leaf in s['accepted']:
            B=tuple(map(F,leaf['box']));actual[leaf['cell'],B]+=1
            assert leaf['depth']==0 and int(leaf['lower_fixed'])>=0 and F(leaf['repair_upper'])==0
            xx,yy=[tuple(map(int,p))for p in leaf['curvature_intervals_fixed']]
            assert xx[0]<=xx[1] and yy[0]<=yy[1]
            err=ceil(F(max(0,xx[1]))*(B[1]-B[0])**2/8)+ceil(F(max(0,yy[1]))*(B[3]-B[2])**2/8)
            assert err==int(leaf['error_fixed'])
            lo=int(leaf['lower_fixed']);min_lower=lo if min_lower is None else min(min_lower,lo)
        assert actual==expected
        artifacts=[pin(base/name)for name in (cname,'source-exact.json','receipt.json','run_pilot.py',mname,sname)]
        rows.append(dict(label=label,root_directory=str(base),archive_directory=str(archive),h=str(h),delta=str(delta),N=N,subdiv=sub,
                         MATLAB_pid=r['pid'],MATLAB_exit_status=0,receipt_pid_absent=True,current_pid_absent=True,
                         MATLAB_wall_seconds=r['wall_seconds'],source_checker_seconds=s['elapsed_seconds'],
                         source_status=s['status'],complete_leaves=count,all_leaf_bounds_nonnegative=True,
                         geometric_initial_partition_exactly_covered=True,curvature_error_arithmetic_matched=True,
                         min_stored_lower_fixed=str(min_lower),halfmass_over_pi_exact=str(h*mass),
                         sampled_D_cost_display=c['cost_from_exact_mass_formula_float'],D_minimum_certified=False,
                         root_to_archive_three_files_byte_identical=True,artifacts=artifacts))
    oldmanifests=[]
    for name in ('round062-coupled-source-independent-manifest.json','round062-v2-delta-independent-manifest.json','round062-v3-jump-independent-manifest.json'):
        path=HERE/'results'/name;m=read(path)
        for item in m['files']:assert sha(Path(item['source']))==item['sha256']
        oldmanifests.append(pin(path))
    for rel,digest in PROOFS.items():assert sha(OWNER/rel)==digest
    out=dict(status='PASS_FINAL_R62_BINDING_AND_COVERAGE_AUDIT',python=sys.version,executable=sys.executable,
             script_sha256=sha(Path(__file__)),elapsed_seconds=time.perf_counter()-start,
             final_proofs=[pin(OWNER/rel)for rel in PROOFS],frozen_independent_review_manifests=oldmanifests,
             runs=rows,total_source_leaves=sum(r['complete_leaves']for r in rows),
             scope='Artifact hashes, receipts, exact saved candidate masses, complete stored geometric leaf coverage and stored curvature-error arithmetic checked. NO full source residual reevaluation or independent rerun of optimized candidates; no MATLAB launch or full-D certificate.')
    assert out['total_source_leaves']==6732
    target=HERE/'results/round062-final-binding-independent.json';target.write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({k:v for k,v in out.items()if k not in ('runs','frozen_independent_review_manifests')},indent=2))


if __name__=='__main__':main()
