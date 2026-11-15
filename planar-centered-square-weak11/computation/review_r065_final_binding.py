#!/usr/bin/env python3
"""R65 final record/partition/mass audit; no source residual reevaluation."""
from collections import Counter,defaultdict
from fractions import Fraction as F
import hashlib,json,os
from pathlib import Path
import subprocess,sys,time

HERE=Path(__file__).resolve().parent
ROOT=Path('/Volumes/data/CodexProjects/T2')
OWNER=ROOT/'.local/workspaces/pcw-sol-kernel-sources/planar-centered-square-weak11'
GRID=OWNER/'computation/inputs/r065-tensor-occupation/central-fine-outer-coarse.json'
V3='6121279c0386bc614bec6c949a1451016954693c3486dd9563e341f377a8d376'
CHECK='5124b452b656fb093a30f40b619bfbc31d7d3bb095385aa6b10f8e0fccc23e38'
V2='3d696e00cbbd4edb0c0a4392efd820c763bb54ac3d8c73f991d3d530418bbfe2'
PROOFS={'r065-nonuniform-tensor-occupation-interface.md':'fd78d51d9e1256e23df6d3ec7d13b7167cc2c601cead6f16cbb2751f82ebfcb5',
        'r065-solver-diagnostic-scope.md':'3b3a9818eeb82715dac35f7d98faf0d28ac8919d783ee4159064368cee92d369',
        'r065-canonical-jump-group-interface.md':'36b61d980fc6f3c6e6833fdaa2b1ecc53985362282df2bd8a9842748866bacbc'}

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def pin(p):return dict(path=str(p),sha256=sha(p),bytes=p.stat().st_size)
def read(p):return json.loads(p.read_text(),parse_float=str)
def ceil(q):return -((-q.numerator)//q.denominator)
def absent(pid):return not subprocess.run(['ps','-p',str(pid),'-o','pid=,comm='],capture_output=True,text=True).stdout.strip()


def main():
    assert os.path.ismount('/Volumes/data')
    start=time.perf_counter()
    for name,h in PROOFS.items():assert sha(OWNER/'proofs'/name)==h
    tokens=read(GRID)['axis_edges'];ed=list(map(F,tokens));centers=[(a+b)/2 for a,b in zip(ed,ed[1:])];widths=[b-a for a,b in zip(ed,ed[1:])]
    groups=defaultdict(list)
    for i,x in enumerate(centers):
        for j,y in enumerate(centers):groups[tuple(sorted((abs(x),abs(y)),reverse=True))].append((i,j))
    keys=sorted(groups);assert len(keys)==276 and len(widths)==45
    width_at_center=dict(zip(centers,widths));rows=[];artifacts=[]
    for rootname,archname,hval,chash in [
        ('r065-v3-solve','v3-h0075-s2',F(3,40),'86fcd67753fbc72714870ae3f1df3766fbd1e99d5c3ef730b3e5cb5478936576'),
        ('r065-v3-h015','v3-h015-s2',F(3,20),'69c8a802d82e57eb46eb60cc8f570be6105699a385c00b4284e68e25cfa41f7f')]:
        p=ROOT/'.local/compute'/rootname;ap=OWNER/'computation/results/r065-tensor-occupation'/archname
        c=read(p/'occupation-tensor-jump-lp.json');s=read(p/'source-certificate.json');r=read(p/'receipt.json');sr=read(p/'source-receipt.json');status=read(p/'tensor-lp-solver-status.json')
        assert sha(p/'occupation-tensor-jump-lp.json')==r['output_sha256']==s['input_sha256']==chash
        assert sha(p/'cauchy_occupation_tensor_jump_lp_v3.m')==sha(OWNER/'computation/cauchy_occupation_tensor_jump_lp_v3.m')==r['source_sha256']==V3
        assert sha(p/'r065_tensor_source_verify.py')==sha(OWNER/'computation/r065_tensor_source_verify.py')==s['script_sha256']==CHECK
        assert sha(p/'grid.json')==sha(GRID)and c['axis_edges_exact']==tokens
        assert F(c['h'])==F(s['h'])==hval and c['subdiv']==s['limits']['initial_subdiv']==2
        assert s['axis_edges_exact']==list(map(str,ed))and s['axis_intervals']==45 and s['D4_orbits']==276 and s['physical_cells']==2025
        for rr in (r,sr):assert rr['exit_status']==0 and rr['owned_pid_absent']is True and absent(rr['pid'])
        assert sr['command'][2]==str(p/'r065_tensor_source_verify.py')
        assert str(p/'occupation-tensor-jump-lp.json')in sr['command']and str(p/'source-certificate.json')in sr['command']
        printed=json.loads(sr['stdout']);assert printed['status']==s['status']and printed['accepted_leaves']==s['accepted_leaves']
        assert Path(s['input'])==p/'occupation-tensor-jump-lp.json'
        assert s['status']=='PASS_SOURCE_NONNEGATIVE'and s['complete_terminal_coverage']is True
        assert s['source_only']is True and s['D_minimum_certified']is False
        assert s['checked_boxes']==s['accepted_leaves']==len(s['accepted'])==1104
        assert s['unchecked_stack']==s['unresolved_leaves']==0 and not any(s[k]for k in ('unchecked','unresolved','split_records','negative_corner_witnesses'))
        assert F(s['origin_repair_added_exact'])==F(s['additional_origin_repair_sufficient_exact'])==0
        assert s['bits']==160 and s['atan_terms']==88
        lam=F(c['source_origin_mass']);a=list(map(F,c['density_coefficients']));assert lam>=0 and len(a)==276 and min(a)>=0
        assert F(s['source_origin_mass_exact'])==lam
        assert [tuple(map(F,v))for v in c['orbit_representatives']]==keys
        repsize=[];repwidth=[];reparea=[];expected=Counter();mass=lam
        for target,key in enumerate(keys):
            wx,wy=map(width_at_center.get,key);n=len(groups[key]);repsize.append(n);repwidth.append((wx,wy));reparea.append(wx*wy)
            mass+=sum(a[target]*widths[i]*widths[j]for i,j in groups[key])
            for ix in range(2):
                for iy in range(2):
                    xl=key[0]-wx/2+wx*ix/2;yl=key[1]-wy/2+wy*iy/2
                    expected[target,(xl,xl+wx/2,yl,yl+wy/2)]+=1
        assert repsize==c['orbit_sizes']and repwidth==[tuple(map(F,v))for v in c['orbit_widths']]
        assert reparea==list(map(F,c['orbit_cell_areas']))and sum(repsize)==2025
        assert F(s['mu_mass_exact'])==mass and F(s['halfmass_over_pi_exact'])==hval*mass
        actual=Counter();min_lower=None
        for leaf in s['accepted']:
            B=tuple(map(F,leaf['box']));actual[leaf['cell'],B]+=1
            assert leaf['depth']==0 and int(leaf['lower_fixed'])>=0 and F(leaf['repair_upper'])==0
            xx,yy=[tuple(map(int,q))for q in leaf['curvature_intervals_fixed']];assert xx[0]<=xx[1]and yy[0]<=yy[1]
            error=ceil(F(max(0,xx[1]))*(B[1]-B[0])**2/8)+ceil(F(max(0,yy[1]))*(B[3]-B[2])**2/8)
            assert error==int(leaf['error_fixed'])
            value=int(leaf['lower_fixed']);min_lower=value if min_lower is None else min(min_lower,value)
        assert actual==expected
        assert status['exitflag']==c['LP_exitflag']==1 and status['LP_output']['iterations']==c['LP_output']['iterations']
        diag=read(p/'tensor-matrix-diagnostic.json')
        assert diag['physical_interfaces']==4140 and diag['canonical_jump_groups']==529
        assert diag['inequality_matrix']['rows']==6880 and diag['inequality_matrix']['columns']==3543
        assert diag['equality_matrix']['rows']==529 and diag['source_matrix']['rows']==4416
        assert c['LP_variables']==3543 and c['canonical_jump_groups']==529
        names=['receipt.json','source-receipt.json','occupation-tensor-jump-lp.json','source-certificate.json','tensor-matrix-diagnostic.json','tensor-lambda-only-diagnostic.json','tensor-lp-solver-status.json']
        for name in names:assert (p/name).read_bytes()==(ap/name).read_bytes();artifacts.append(pin(p/name))
        rows.append(dict(root_directory=str(p),archive_directory=str(ap),h=str(hval),candidate_sha256=chash,
                         source_certificate_sha256=sha(p/'source-certificate.json'),MATLAB_pid=r['pid'],MATLAB_exit=0,MATLAB_absent_now=True,
                         MATLAB_wall_seconds=r['wall_seconds'],source_pid=sr['pid'],source_exit=0,source_absent_now=True,
                         source_certificate_seconds=s['elapsed_seconds'],source_receipt_wall_seconds=sr['wall_seconds'],
                         solver_iterations=status['LP_output']['iterations'],accepted_leaves=1104,complete_reconstructed_partition=True,
                         all_stored_lower_bounds_nonnegative=True,min_stored_lower_fixed=str(min_lower),anisotropic_error_arithmetic_exact=True,
                         nonuniform_mass_exact=str(mass),halfmass_over_pi_exact=str(hval*mass),sampled_D_cost=c['cost_from_exact_mass_formula_float'],
                         source_status=s['status'],zero_repair=True,D_minimum_certified=False,seven_archive_copies_byte_identical=True))
    legacy=[]
    for rootname,mname,expected_sha,exitcode in [('r065-h0075','cauchy_occupation_tensor_jump_lp.m','5f00c7f4d4231e5353e5de189a4e116758cdfe4e0d022fd03cd1dde35433f21b',1),('r065-lambda-diagnostic','cauchy_occupation_tensor_jump_lp_v2.m',V2,0),('r065-v2-solve','cauchy_occupation_tensor_jump_lp_v2.m',V2,1)]:
        p=ROOT/'.local/compute'/rootname;r=read(p/'receipt.json')
        assert sha(p/mname)==r['source_sha256']==expected_sha
        gridfile=p/('central-fine-outer-coarse.json'if rootname=='r065-h0075'else'grid.json')
        assert sha(gridfile)==sha(GRID)
        assert r['exit_status']==exitcode and r['owned_pid_absent']is True and absent(r['pid'])
        rec=dict(root_directory=str(p),MATLAB_pid=r['pid'],exit_status=exitcode,pid_absent_now=True,source_sha256=expected_sha)
        if rootname=='r065-h0075':
            assert not (p/'tensor-lp-solver-status.json').exists();rec['solver_exitflag']='UNKNOWN_NOT_RECONSTRUCTED'
        elif rootname=='r065-lambda-diagnostic':
            d=read(p/'tensor-lambda-only-diagnostic.json');assert sha(p/'tensor-lambda-only-diagnostic.json')==r['output_sha256']
            assert d['axis_edges_exact']==tokens and F(d['h'])==F(3,40)and d['subdiv']==2
            assert d['source_rows']==4416 and F(d['source_row_min_unit_lambda'])>0 and d['numeric_full_sampled_LP_witness_exists']
            rec.update(source_rows=4416,floating_minrow=d['source_row_min_unit_lambda'],exact_source_or_fullD_certificate=False)
            artifacts.append(pin(p/'tensor-lambda-only-diagnostic.json'))
        else:
            d=read(p/'tensor-lp-solver-status.json');assert d['exitflag']==0 and d['LP_output']['iterations']==20000
            assert d['returned_vector']==[] and d['objective']==[]
            diag=read(p/'tensor-matrix-diagnostic.json')
            assert diag['LP_MaxIterations']==20000 and diag['source_matrix']['rows']==4416
            assert diag['inequality_matrix']['columns']==10765 and diag['equality_matrix']['rows']==4140
            rec.update(solver_exitflag=0,iterations=20000,classification='ITERATION_CAP_ONLY_NO_INFEASIBILITY_CERTIFICATE')
            artifacts.append(pin(p/'tensor-lp-solver-status.json'))
        legacy.append(rec);artifacts.append(pin(p/'receipt.json'))
    for name,h in PROOFS.items():assert sha(OWNER/'proofs'/name)==h
    old=[]
    for name in ('round065-tensor-interface-independent-manifest.json','round065-v3-group-independent-manifest.json'):
        p=HERE/'results'/name;mm=read(p)
        for f in mm['files']:assert sha(Path(f['source']))==f['sha256']
        old.append(pin(p))
    out=dict(status='PASS_R65_FINAL_ARTIFACT_DOMAIN_BINDING',python=sys.version,executable=sys.executable,
             script_sha256=sha(Path(__file__)),elapsed_seconds=time.perf_counter()-start,
             final_proofs=[pin(OWNER/'proofs'/n)for n in PROOFS],frozen_review_manifests=old,runs=rows,legacy_diagnostics=legacy,
             artifacts=artifacts,total_complete_source_leaves=2208,
             scope='Saved artifact/hash/solver/source-receipt checks, exact orbit partition/mass and stored anisotropic error arithmetic. No reevaluation of all source values/derivative enclosures, no MATLAB run, no fullD or newC2 certificate.')
    p=HERE/'results/round065-final-binding-independent.json';p.write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({k:v for k,v in out.items()if k not in ('artifacts','runs','legacy_diagnostics','frozen_review_manifests')},indent=2))


if __name__=='__main__':main()
