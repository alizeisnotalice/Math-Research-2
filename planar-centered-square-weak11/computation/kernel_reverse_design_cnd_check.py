#!/usr/bin/env python3
"""Exact small four-point checks for the reciprocal diamond potential."""
from fractions import Fraction as F
from pathlib import Path
import datetime
import hashlib
import json
import os
import platform
import subprocess
import sys


def fs(x):
    return str(x.numerator) + '/' + str(x.denominator)


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def candidate_value(p, x_num, y_num, denom):
    # Coordinates are 0 or +/- 2^-denom or +/-1. Powers are rational because
    # p*denom and (2-p)*denom are integers for these deliberately chosen inputs.
    ax, ay = abs(x_num), abs(y_num)
    if ax == 0 and ay == 0:
        return F(0)
    if ax == 0:
        assert ay == 1
        return F(1)
    assert ax == 1
    a = F(1, 2**denom)
    ap = F(1, 2**p.numerator)
    aa = F(1, 2**(2*denom-p.numerator))
    if ay == 0:
        return aa
    assert ay == 1
    return (1-a*a)/(1-ap)


def main():
    assert os.path.ismount('/Volumes/data'), 'data volume not mounted'
    base=Path(__file__).resolve().parent.parent
    inputs=[F(1,4),F(1,2),F(3,4),F(1),F(5,4),F(3,2),F(7,4)]
    rows=[]
    z=[(0,0),(1,0),(0,1),(1,1)]
    coeff=[1,-1,-1,1]
    for p in inputs:
        d=p.denominator
        a=F(1,2**d)
        ap=F(1,2**p.numerator)
        aalpha=F(1,2**(2*d-p.numerator))
        diagonal=candidate_value(p,1,1,d)
        defect=diagonal-1-aalpha
        formula=(ap-aalpha)/(1-ap)
        assert defect == formula
        qform=sum(F(coeff[i]*coeff[j])*candidate_value(
            p,z[i][0]-z[j][0],z[i][1]-z[j][1],d)
            for i in range(4) for j in range(4))
        assert qform==4*defect
        assert (qform>0) == (p<1)
        assert (qform==0) == (p==1)
        rows.append({'p':fs(p),'alpha':fs(2-p),'a':fs(a),
                     'phi_a_1':fs(diagonal),'phi_a_0':fs(aalpha),
                     'phi_0_1':'1/1','defect':fs(defect),'qform':fs(qform),
                     'four_point_rejects_cnd':qform>0})
    p=F(1,2); r=F(16,9); c=F(4)
    # sqrt(r)=4/3, so c=(1-r^-p)^-1=4 and r^(2-p)=64/27.
    integral_half = c*p*F(64,27)/(2-p)
    optimized_cost = F(4,3)**4
    assert integral_half==optimized_cost==F(256,81)
    p1_half=F(2)*F(1)*F(2)/(2-F(1))
    assert p1_half==4
    assert rows[1]['defect']=='3/4' and rows[1]['qform']=='3/1'
    # Exact hinge inequalities for two p with integral hinge exponent beta.
    hinge_checks=[]
    for p in [F(1,2),F(1)]:
        beta=(2-p)/p
        assert beta.denominator==1
        constant=F(256,81) if p==F(1,2) else F(4)
        for t in [F(1,4),F(1,2),F(3,4),F(1),F(5,4)]:
            j=(1/beta)*t**(-beta.numerator)
            boundary=max(F(0),1-t)
            gap=j-constant*boundary
            assert gap>=0
            assert (gap==0)==(t==1-p/2)
            if p==1 and t<1:
                assert gap==(2*t-1)**2/t
            hinge_checks.append({'p':fs(p),'t':fs(t),'hinge_half_integral':fs(j),
                                 'boundary_value':fs(boundary),'cost_floor':fs(constant),
                                 'gap':fs(gap)})
    # Local cusp integral coefficient, independently simplified with rationals.
    cusp=[]
    for alpha in [F(5,4),F(3,2),F(7,4)]:
        direct=4*(1/(alpha-1)-1/alpha)
        closed=4/(alpha*(alpha-1))
        assert direct==closed>0
        cusp.append({'alpha':fs(alpha),'full_line_abs_cusp_coefficient':fs(closed),
                     'half_line_abs_cusp_coefficient':fs(closed/2)})
    repo=base.parent
    gitbase=subprocess.run(['git','-C',str(repo),'rev-parse','HEAD'],
                           capture_output=True,text=True,check=True).stdout.strip()
    deps=['proofs/kernel-interface-audit.md','proofs/kernel-family-floor.md',
          'proofs/kernel-reverse-design-cnd-obstruction.md',
          'proofs/kernel-subcritical-cusp-stable-obstruction-review.md']
    hashes={name:sha(base/name) for name in deps if (base/name).exists()}
    output={'status':'PASS: exact finite arithmetic checks only',
            'timestamp_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),
            'python':sys.version,'executable':sys.executable,'platform':platform.platform(),
            'git_base':gitbase,'script_sha256':sha(Path(__file__).resolve()),
            'dependencies_sha256':hashes,
            'audited_author_draft':{'path':'/Volumes/data/CodexProjects/T2/.local/workspaces/pcw-astra-kernel/planar-centered-square-weak11/proofs/kernel-subcritical-cusp-obstruction.md',
                                    'sha256_at_read':'02ed862985cce2ce8bea7aac85b6c3cb396e5b9a1bad1f3f1246116f4d263e6d',
                                    'prior_full_review_sha256':'8d8b296824a2e545bde73ebb0c95b4de39f6a6a7e47392458482c3ac8b4eb3b0'},
            'arithmetic':'fractions.Fraction, exact',
            'workers':1,'random_seed':None,'floating_tolerance':None,
            'four_point_coefficients':coeff,'point_encoding':'(x*a,y), a=2^-denominator(p)',
            'four_point_points':z,'checks':rows,
            'integral_checks':{'p_half_R':'16/9','p_half_c':'4/1',
                               'p_half_J':fs(integral_half),'p_one_R':'2/1',
                               'p_one_c':'2/1','p_one_J':fs(p1_half)},
            'axis_cusp_coefficients':cusp,'hinge_checks':hinge_checks,
            'scope':'Only exact arithmetic identities and finite CND witnesses; analytic proof supplies parameter-wide claims. p>1 passing this necessary test is not a CND certificate. No C4 feasibility test or R9 matrix rerun.'}
    out=base/'computation/results/kernel-reverse-design/exact-four-point-checks.json'
    out.parent.mkdir(parents=True,exist_ok=True)
    out.write_text(json.dumps(output,indent=2)+'\n')
    print(json.dumps({'status':output['status'],'output':str(out),'p_half_qform':rows[1]['qform'],
                      'p_half_J':fs(integral_half),'p_one_J':fs(p1_half)},indent=2))

if __name__=='__main__':
    main()
