"""Small analytic boundary checks for the read-only general-dual probe.

No random search or LP optimization is run. Expected coefficients below are
obtained by elementary two-point/collinear geometry, independently of the
probe's radius-layer construction.
"""
from fractions import Fraction as F
from pathlib import Path
import argparse
from datetime import datetime,timezone
import hashlib,json,platform,runpy,sys,time


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--module',type=Path,required=True)
    ap.add_argument('--output',type=Path,required=True)
    args=ap.parse_args()
    start=time.perf_counter()
    raw=args.module.read_bytes()
    module=runpy.run_path(str(args.module),run_name='general_dual_audit_import')
    coefficients,masks=module['coefficients'],module['masks']
    cases=[
        ('one_atom_unattained_capacity',[(F(0),F(0))],[F(3,7)],{(0,):F(3,7)}),
        ('two_atoms_exact_strict_threshold',[(F(0),F(0)),(F(2),F(1))],
         [F(1),F(3)],{(0,):F(1),(1,):F(3)}),
        ('two_atoms_above_threshold',[(F(0),F(0)),(F(2),F(1))],
         [F(1),F(13,4)],{(0,):F(1),(1,):F(13,4),(0,1):F(17,4)}),
        ('bounded_interior_singleton_capacity',[(F(0),F(0)),(F(1),F(0)),(F(2),F(0))],
         [F(1,4),F(10),F(1,4)],
         {(0,):F(1,4),(1,):F(4),(2,):F(1,4),
          (0,1):F(41,4),(1,2):F(41,4),(0,1,2):F(21,2)}),
    ]
    result=[]
    for name,points,masses,expected in cases:
        actual=coefficients(points,masses)
        assert actual==expected,(name,actual,expected)
        result.append({'name':name,'points':points,'masses':masses,
                       'exact_coefficients':{','.join(map(str,k)):v for k,v in sorted(actual.items())},
                       'matches_independent_analytic_expectation':True})
    diagonal=[(F(-1),F(-1)),(F(1),F(1))]
    assert (0,1) not in masks(diagonal,F(999,1000))
    assert (0,1) in masks(diagonal,F(1))
    assert (0,1) in masks(diagonal,F(1001,1000))
    collinear=[(F(0),F(0)),(F(1),F(0)),(F(2),F(0))]
    assert (0,1,2) in masks(collinear,F(1))
    assert (0,2) not in masks(collinear,F(1))
    # Fixed-family source-load/all-space distinction; no conclusion on C2.
    N=3
    ts=[F(i,N+1) for i in range(1,N+1)]
    centers=[(t,1-t) for t in ts]
    sources=[(1+t,2-t) for t in ts]
    incidences=[tuple(j for j,p in enumerate(sources)
                      if max(abs(x-p[0]),abs(y-p[1]))<=1) for x,y in centers]
    assert incidences==[(i,) for i in range(N)]
    assert all(max(abs(x),abs(y))<=1 for x,y in centers)
    out={'python':sys.version,'platform':platform.platform(),
         'run_at_utc':datetime.now(timezone.utc).isoformat(),
         'probe_path':str(args.module),'probe_sha256':hashlib.sha256(raw).hexdigest(),
         'audit_script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
         'small_exact_cases':result,'closed_critical_radius_and_zero_area_center_cell':'passed',
         'fixed_family_example':{'N':N,'source_incidence_labels':incidences,
                                 'source_load_objective':4*N,'all_space_objective':4,
                                 'scope':'only this selected family; not a weak-constant obstruction'},
         'elapsed_seconds':time.perf_counter()-start,
         'command':[sys.executable,*sys.argv],
         'full_100_case_search_rerun':False}
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(out,default=str,indent=2)+'\n')
    print(json.dumps(out,default=str,indent=2))

if __name__=='__main__':
    main()
